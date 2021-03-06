library(h2o)
library(ukbtools)
library(tidyverse)
library(dplyr)

## Create training and validation frames
# Get CNV Scale Data
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed_4splits.txt", sep = " ")

# COVID data
covid_results <- read.csv("/data/ukbiobank/covid-19/covid19_result_June19-2020.txt", sep = "\t")
covid_results <- covid_results[covid_results$result == 1, ] # Get positive results
covid_results <- covid_results[,!names(covid_results) %in% c("datereported", "specdate", "spectype", "laboratory", "origin")] # Remove unnecessary columns
covid_results <- unique(covid_results) # Remove duplicate tests


# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data_age, by.x = "ids", by.y = "eid")

# Get COVID patients
covid_data <- merge(all_data, covid_results, by.x = "ids", by.y = "eid")
covid <- covid_data[covid_data$result == 1,]
covid <- unique(covid)

# Get breakdown of COVID patients' age
covid_age <- table(covid$yearBorn)

# As an age control we will only look at individuals born within the same time as the COVID patients for our non AD patients
all_data <- subset(all_data, all_data$yearBorn < max(covid$yearBorn))
all_data <- subset(all_data, all_data$yearBorn > min(covid$yearBorn))

#Get non COVID patients
no_covid_initial <- all_data[!all_data$ids %in% covid$ids,]

# Randomly get non COVID patients for controls so that there is an equal amount based on age
# This will ensure that the controls are age-matched to the COVID Disease sample
# For example if there are 5 patients born 1937 who have covid we will randomly grab 5 other 
# patients born 1937 who do not have covid
no_covid <- data.frame(matrix(ncol = ncol(no_covid_initial), nrow = 0))
colnames(no_covid) <- colnames(no_covid_initial)
for (i in 1:length(covid_age)) {
  temp <- covid_age[i]
  age_check <- as.numeric(names(temp))
  number_cases <- as.numeric(unname(temp))
  possible_controls <- no_covid_initial[no_covid_initial$yearBorn == age_check,]
  no_covid <- rbind(no_covid, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}
no_covid <- unique(no_covid) # Remove duplicates if they exist

ind <- sample(c(TRUE, FALSE), nrow(covid), replace=TRUE, prob=c(0.7, 0.3)) # Random split

# Split COVID patients into train and validate
train <- covid[ind, ]
validate <- covid[!ind, ]

#Remove unnecessary rows
train <- train[,!names(train) %in% c("datereported", "specdate", "spectype", "laboratory", "origin")]
validate <- validate[,!names(validate) %in% c( "datereported", "specdate", "spectype", "laboratory", "origin")]

controls <- no_covid # get controls
controls['result'] = 0 # add column indicating these patients do not have COVID

# Split control patients into train and validate
train_controls <- controls[ind, ]
validate_controls <- controls[!ind, ]

# Combine controls with samples

train <- rbind(train, train_controls)
validate <- rbind(validate, validate_controls)

# Set response column to boolean
train$result <- as.logical(train$result)
validate$result <- as.logical(validate$result)
train$result <- as.factor(train$result)
validate$result <- as.factor(validate$result)

#Remove unnecessary columns
train <- train[,!names(train) %in% c("ids")]
validate <- validate[,!names(validate) %in% c( "ids")]

# Free up data 
rm(no_covid, covid, controls, train_controls, validate_controls)
rm(covid_results, condensed, my_ukb_data_cancer, my_data_age, all_data)

# Load h2o
h2o.init(nthreads=15)

# Load data into h2o
train.hex <- as.h2o(train, destination_frame = "train.hex")  
validate.hex <- as.h2o(validate, destination_frame = "validate.hex")

#
#  I usually stop here and goto http://localhost:54321/flow/index.html 
#  h2o runs a local webserver on port 54321, it offers a nice little interface.
#  you can run the AutoML from the web browser there and it has some nice features so that 
# you can monitor the progress of the training.

#Response column
response <- "result"
#Get Predictors
predictors <- colnames(train)
predictors <- predictors[! predictors %in% response] #Response cannot be a predictor
predictors <- predictors[! predictors %in% "yearBorn"] #Remove yearBorn as predictor
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
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("COVID Cross-Validated 4 Splits",auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Attempt predict on validation frame
prediction <- h2o.predict(object = leader, newdata = validate.hex)
as.data.frame(prediction)
summary(prediction, exact_quantiles = TRUE)

validation.perf <- h2o.performance(leader, train = FALSE, xval=TRUE, newdata = validate.hex)
validation.perf.auc <- validation.perf@metrics$AUC
plot(validation.perf, type = 'roc', main = paste("COVID 4 Splits predictions", validation.perf.auc))


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

# COVID data
covid_results <- read.csv("/data/ukbiobank/covid-19/covid19_result.txt", sep = "\t")
covid_results <- covid_results[covid_results$result == 1, ] # Get positive results
covid_results <- covid_results[,!names(covid_results) %in% c("datereported", "specdate", "spectype", "laboratory", "origin")] # Remove unnecessary columns
covid_results <- unique(covid_results) # Remove duplicate tests

# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data_age, by.x = "ids", by.y = "eid")

# Get COVID patients
covid_data <- merge(all_data, covid_results, by.x = "ids", by.y = "eid")
covid <- covid_data[covid_data$result == 1,]
covid <- unique(covid)

# Get breakdown of COVID patients' age
covid_age <- table(covid$yearBorn)

# As an age control we will only look at individuals born within the same time as the COVID patients for our non AD patients
all_data <- subset(all_data, all_data$yearBorn < max(covid$yearBorn))
all_data <- subset(all_data, all_data$yearBorn > min(covid$yearBorn))

#Get non COVID patients
no_covid_initial <- all_data[!all_data$ids %in% covid$ids,]

for (i in 1:numModels) {
  # Randomly get non COVID patients for controls so that there is an equal amount based on age
  # This will ensure that the controls are age-matched to the COVID Disease sample
  # For example there are 5 patients born 1937 who have covid so we will randomly grab 5 other 
  # patients born 1937 who do not have covid
  no_covid <- data.frame(matrix(ncol = ncol(no_covid_initial), nrow = 0))
  colnames(no_covid) <- colnames(no_covid_initial)
  for (i in 1:length(covid_age)) {
    temp <- covid_age[i]
    age_check <- as.numeric(names(temp))
    number_cases <- as.numeric(unname(temp))
    possible_controls <- no_covid_initial[no_covid_initial$yearBorn == age_check,]
    no_covid <- rbind(no_covid, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
  }
  no_covid <- unique(no_covid)
  
  ind <- sample(c(TRUE, FALSE), nrow(covid), replace=TRUE, prob=c(0.7, 0.3)) # Random split
  
  train <- covid[ind, ]
  validate <- covid[!ind, ]
  
  #Remove unnecessary rows
  train <- train[,!names(train) %in% c("datereported", "specdate", "spectype", "laboratory", "origin")]
  validate <- validate[,!names(validate) %in% c( "datereported", "specdate", "spectype", "laboratory", "origin")]
  
  controls <- no_covid # Get controls
  controls['result'] = 0
  
  train_controls <- controls[ind, ]
  validate_controls <- controls[!ind, ]
  
  # Combine controls with samples
  
  train <- rbind(train, train_controls)
  validate <- rbind(validate, validate_controls)
  
  # Set response column to boolean
  train$result <- as.logical(train$result)
  validate$result <- as.logical(validate$result)
  train$result <- as.factor(train$result)
  validate$result <- as.factor(validate$result)
  
  #Remove unnecessary columns
  train <- train[,!names(train) %in% c("ids")]
  validate <- validate[,!names(validate) %in% c( "ids")]
  
  train <- train[!is.na(train$result),]
  validate <- validate[!is.na(validate$result),]
  
  # Load data into h2o
  
  train.hex <- as.h2o(train, destination_frame = "train.hex")  
  validate.hex <- as.h2o(validate, destination_frame = "validate.hex")
  
  #
  #  I usually stop here and goto http://localhost:54321/flow/index.html 
  #  h2o runs a local webserver on port 54321, it offers a nice little interface.
  #  you can run the AutoML from the web browser there and it has some nice features so that 
  # you can monitor the progress of the training.
  
  #Response column
  response <- "result"
  #Get Predictors
  predictors <- colnames(train)
  predictors <- predictors[! predictors %in% response] #Response cannot be a predictor
  predictors <- predictors[! predictors %in% "yearBorn"] #Remove yearBorn predictor
  model <- h2o.automl(x = predictors,
                      y = response,
                      training_frame = train.hex,
                      validation_frame = validate.hex,
                      nfolds=10,
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
