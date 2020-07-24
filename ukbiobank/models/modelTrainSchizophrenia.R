library(h2o)
library(ukbtools)
library(tidyverse)
library(dplyr)
library(ggplot2)

# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed_4splits.txt", sep = " ")

# Bipolar and Schizophrenia data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")
my_data <- select(my_ukb_data,eid,
                  datereported = date_f20_first_reported_schizophrenia_f130874_0_0,
                  sourcereported = source_of_report_of_f20_schizophrenia_f130875_0_0)

# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")
all_data <- merge(all_data, my_data_age, by.x = "ids", by.y = "eid")

schiz <- all_data[!is.na(all_data[, "datereported"]),]
no_schiz_initial <- all_data[is.na(all_data[, "datereported"]),]

# Get breakdown of  patients by age
schiz_age <- table(schiz$yearBorn)

# Randomly get non disease patients for controls so that there is an equal amount based on age
# This will ensure that the controls are age-matched to the disease sample
# For example there are 5 patients born 1937 who have AD so we will randomly grab 5 other 
# patients born 1937 who do not have AD
no_schiz <- data.frame(matrix(ncol = ncol(no_schiz_initial), nrow = 0))
colnames(no_schiz) <- colnames(no_schiz_initial)
for (i in 1:length(schiz_age)) {
  temp <- schiz_age[i]
  age_check <- as.numeric(names(temp))
  number_cases <- as.numeric(unname(temp))
  possible_controls <- no_schiz_initial[no_schiz_initial$yearBorn == age_check,]
  no_schiz <- rbind(no_schiz, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}

schiz$datereported <- TRUE
no_schiz$datereported <- FALSE

ind <- sample(c(TRUE, FALSE), nrow(schiz), replace=TRUE, prob=c(0.7, 0.3)) # Random split

train <- schiz[ind, ]
validate <- schiz[!ind, ]

controls <- no_schiz #  get controls

train_controls <- controls[ind, ]
validate_controls <- controls[!ind, ]

# Combine controls with samples
train <- rbind(train, train_controls)
validate <- rbind(validate, validate_controls)

# Set response column to factor
train$datereported <- as.factor(train$datereported)
validate$datereported <- as.factor(validate$datereported)

#Remove unnecessary columns
train <- train[,!names(train) %in% c("ids", "sex", "behavior")]
validate <- validate[,!names(validate) %in% c("ids", "sex", "behavior")]

# Free up data 
rm(schiz, no_schiz, controls, train_controls, validate_controls)
rm(my_data, my_ukb_data, my_ukb_data_cancer, my_data_age)

# Load data into h2o

train.hex <- as.h2o(train, destination_frame = "train.hex")  
validate.hex <- as.h2o(validate, destination_frame = "validate.hex")

#
#  I usually stop here and goto http://localhost:54321/flow/index.html 
#  h2o runs a local webserver on port 54321, it offers a nice little interface.
#  you can run the AutoML from the web browser there and it has some nice features so that 
# you can monitor the progress of the training.

#Response column
response <- "datereported"
#Get Predictors
predictors <- colnames(train)
predictors <- predictors[! predictors %in% response] #Response cannot be a predictor
predictors <- predictors[! predictors %in% "yearBorn"] #Response cannot be a predictor
predictors <- predictors[! predictors %in% "sourcereported"] #Response cannot be a predictor
model <- h2o.automl(x = predictors,
                    y = response,
                    training_frame = train.hex,
                    validation_frame = validate.hex,
                    nfolds=5,
                    max_runtime_secs = 360)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Schizophrenia", auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)


# Graceful shutdown of cluster
h2o.shutdown(prompt = TRUE)

###################################################
# Calculate Confidence Interval for 4 Split Model #
###################################################

# Preallocate vector for aucs
rm(results, predictions, model_types)
results <- c()
predictions <- c()
model_types <- c()
numModels <- 100
maxRuntime <- 360 # This is in seconds

# Run 100 expirements or train 100 Auto ML models using randomized set of training data each time
# Each model will also have 5 fold cross-validation as a base parameter.

####
# This section is for the 4 split model
###
# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
# Get CNV Scale Data
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed_4splits.txt", sep = " ")

# Bipolar and Schizophrenia data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")
my_data <- select(my_ukb_data,eid,
                  datereported = date_f20_first_reported_schizophrenia_f130874_0_0,
                  sourcereported = source_of_report_of_f20_schizophrenia_f130875_0_0)

# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")
all_data <- merge(all_data, my_data_age, by.x = "ids", by.y = "eid")

schiz <- all_data[!is.na(all_data[, "datereported"]),]
no_schiz_initial <- all_data[is.na(all_data[, "datereported"]),]

# Get breakdown of  patients by age
schiz_age <- table(schiz$yearBorn)

for (i in 1:numModels) {
  # Randomly get non disease patients for controls so that there is an equal amount based on age
  # This will ensure that the controls are age-matched to the disease sample
  # For example there are 5 patients born 1937 who have AD so we will randomly grab 5 other 
  # patients born 1937 who do not have AD
  no_schiz <- data.frame(matrix(ncol = ncol(no_schiz_initial), nrow = 0))
  colnames(no_schiz) <- colnames(no_schiz_initial)
  for (i in 1:length(schiz_age)) {
    temp <- schiz_age[i]
    age_check <- as.numeric(names(temp))
    number_cases <- as.numeric(unname(temp))
    possible_controls <- no_schiz_initial %>% filter(yearBorn == age_check)
    no_schiz <- rbind(no_schiz, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
  }
  
  schiz$datereported <- TRUE
  no_schiz$datereported <- FALSE
  
  ind <- sample(c(TRUE, FALSE), nrow(schiz), replace=TRUE, prob=c(0.7, 0.3)) # Random split
  
  train <- schiz[ind, ]
  validate <- schiz[!ind, ]
  
  controls <- no_schiz #  get controls
  
  train_controls <- controls[ind, ]
  validate_controls <- controls[!ind, ]
  
  # Combine controls with samples
  train <- rbind(train, train_controls)
  validate <- rbind(validate, validate_controls)
  
  # Set response column to factor
  train$datereported <- as.factor(train$datereported)
  validate$datereported <- as.factor(validate$datereported)
  
  #Remove unnecessary columns
  train <- train[,!names(train) %in% c("ids", "sex", "behavior")]
  validate <- validate[,!names(validate) %in% c("ids", "sex", "behavior")]
  
  # Load data into h2o
  
  train.hex <- as.h2o(train, destination_frame = "train.hex")  
  validate.hex <- as.h2o(validate, destination_frame = "validate.hex")
  
  #
  #  I usually stop here and goto http://localhost:54321/flow/index.html 
  #  h2o runs a local webserver on port 54321, it offers a nice little interface.
  #  you can run the AutoML from the web browser there and it has some nice features so that 
  # you can monitor the progress of the training.
  
  #Response column
  response <- "datereported"
  #Get Predictors
  predictors <- colnames(train)
  predictors <- predictors[! predictors %in% response] #Response cannot be a predictor
  predictors <- predictors[! predictors %in% "yearBorn"] #Response cannot be a predictor
  predictors <- predictors[! predictors %in% "sourcereported"] #Response cannot be a predictor
  model <- h2o.automl(x = predictors,
                      y = response,
                      training_frame = train.hex,
                      validation_frame = validate.hex,
                      nfolds=5,
                      max_runtime_secs = maxRuntime)
  
  
  #record the Leading model AUC in the dataset
  leader <- model@leader
  auc=h2o.auc(leader, train=FALSE, xval=TRUE)
  results <- c(results, auc)
  model_types <- c(model_types, leader@algorithm)
  # If desired plot out the ROC.  We type out the tissue and AUC at the top of the ROC. Keep in mind this will make many plots
  # plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("COVID 4 Splits w/o yearBorn",auc))
  
  # Print performance info of leader. Note this will also do this for every model
  #leader@algorithm
  #h2o.performance(leader,train=FALSE, xval=TRUE)
  
  # Attempt predict on validation frame
  prediction <- h2o.predict(object = leader, newdata = validate.hex)
  as.data.frame(prediction)
  summary(prediction, exact_quantiles = TRUE)
  
  validation.perf <- h2o.performance(leader, train = FALSE, xval=TRUE, newdata = validate.hex)
  validation.perf.auc <- validation.perf@metrics$AUC
  
  predictions <- c(predictions, validation.perf.auc)
}

# Function to calculate confidence interval
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

confidence_interval(results, 0.90)
confidence_interval(results, 0.95)
confidence_interval(results, 0.99)

confidence_interval(predictions, 0.90)
confidence_interval(predictions, 0.95)
confidence_interval(predictions, 0.99)

# Graceful shutdown of cluster
h2o.shutdown(prompt = TRUE)


# Make plots
aucs <- data.frame(results)
models <- data.frame(model_types)
finalModelResults <- cbind(aucs, models)
ggplot(finalModelResults, aes(x=model_types, y=results)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)

# Create Confidence interval plot
df <- data.frame(x = 1:100,
                 y = results)
plot(y ~ x, data = df)

# model
mod <- lm(y ~ x, data = df)

# predicts + interval
newx <- seq(min(df$x), max(df$x), length.out=100)
preds <- predict(mod, newdata = data.frame(x=newx), 
                 interval = 'confidence')

# plot
plot(y ~ x, data = df, type = 'p')
# add fill
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
# model
abline(mod)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
title("Confidence Interval of AUCs of Schizophrenia Models")