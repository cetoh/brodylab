# Must use h2o v3.32.0.2 or higher for the explainability plots
# Run this section to install latest h2o package
if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
pkgs <- c("RCurl","jsonlite")
for (pkg in pkgs) {
  if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
}
install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R")))

# Add libraries
library(h2o)
library(ukbtools)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggsci)

# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed_4splits.txt", sep = " ")

# AD data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")
my_data <- select(my_ukb_data,eid,
                  datereported = date_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0,
                  sourcereported = source_of_report_of_f00_dementia_in_alzheimers_disease_f130837_0_0)

# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")
all_data <- merge(all_data, my_data_age, by.x = "ids", by.y = "eid")

alzheimers <- all_data[!is.na(all_data[, "datereported"]),]
no_alzheimers_initial <- all_data[is.na(all_data[, "datereported"]),]

# Get breakdown of  patients by age
alzheimers_age <- table(alzheimers$yearBorn)

# Randomly get non disease patients for controls so that there is an equal amount based on age
# This will ensure that the controls are age-matched to the disease sample
# For example there are 5 patients born 1937 who have AD so we will randomly grab 5 other 
# patients born 1937 who do not have AD
no_alzheimers <- data.frame(matrix(ncol = ncol(no_alzheimers_initial), nrow = 0))
colnames(no_alzheimers) <- colnames(no_alzheimers_initial)
for (i in 1:length(alzheimers_age)) {
  temp <- alzheimers_age[i]
  age_check <- as.numeric(names(temp))
  number_cases <- as.numeric(unname(temp))
  possible_controls <- no_alzheimers_initial[no_alzheimers_initial$yearBorn == age_check,]
  no_alzheimers <- rbind(no_alzheimers, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}

alzheimers$datereported <- "Alzheimers"
no_alzheimers$datereported <- "Normal"

ind <- sample(c(TRUE, FALSE), nrow(alzheimers), replace=TRUE, prob=c(0.7, 0.3)) # Random split

train <- alzheimers[ind, ]
validate <- alzheimers[!ind, ]

controls <- no_alzheimers #  get controls

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
rm(alzheimers, no_alzheimers, controls, train_controls, validate_controls)
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
                    nfolds = 5,
                    max_runtime_secs = 300,
                    keep_cross_validation_predictions = TRUE)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# Generate explainability plots
exp_aml <- h2o.explain(model, validate.hex)
exp_aml
exp_leader <- h2o.explain(leader, validate.hex)
exp_leader
h2o.model_correlation_heatmap(model, validate.hex)
h2o.varimp_heatmap(model)
h2o.pd_multi_plot(model, validate.hex, "chrX_1.X1")
h2o.shap_summary_plot(leader, validate.hex)
h2o.ice_plot(leader, validate.hex, "chrX_1.X1")

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
mod_perf <- h2o.performance(leader,train=FALSE, xval=TRUE)


plot(mod_perf, type='roc', main=paste("Alzheimer's Disease 4-split CSLV Model Performance: AUC = ", auc))

# for example I have a list of H2OModels
list(leader) %>% 
  # map a function to each element in the list
  map(function(x) x %>% h2o.performance(valid=T) %>% 
        # from all these 'paths' in the object
        .@metrics %>% .$thresholds_and_metric_scores %>% 
        # extracting true positive rate and false positive rate
        .[c('tpr','fpr')] %>% 
        # add (0,0) and (1,1) for the start and end point of ROC curve
        add_row(tpr=0,fpr=0,.before=T) %>% 
        add_row(tpr=0,fpr=0,.before=F)) %>% 
  # add a column of model name for future grouping in ggplot2
  map2(c('4-Split'),
       function(x,y) x %>% add_column(model=y)) %>% 
  # reduce four data.frame to one
  reduce(rbind) %>% 
  # plot fpr and tpr, map model to color as grouping
  ggplot(aes(fpr,tpr,col=model))+
  geom_line()+
  geom_segment(aes(x=0,y=0,xend = 1, yend = 1),linetype = 2,col='grey')+
  xlab('False Positive Rate')+
  ylab('True Positive Rate')+
  ggtitle("ROC Curve for 4-Split Alzheimer's Disease Model")+
  annotate("text", x = .75, y = .25, label = paste("AUC =", auc)) + 
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()

# Plot Precision-Recall Graph
metrics <- as.data.frame(h2o.metric(mod_perf))
metrics %>%
  ggplot(aes(recall,precision)) + 
  geom_line() +
  theme_minimal() + ggtitle("Alzheimer's 4-split CSLV Model Precision-Recall") +
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Find Odds ratio
#cross_val_preds <- gbm@model[["cross_validation_predictions"]]

cvpreds <- as.data.frame(h2o.getFrame(model@leader@model[["cross_validation_holdout_predictions_frame_id"]][["name"]]))

#combine with the original data
predictions <- cbind(cvpreds,train)

#rank by disease predicted score ascending
predictions$rank <- rank(predictions[,"datereported"])
predictions<- select(predictions, rank, everything())

#total number of AD cases
num_patients <- nrow(predictions[predictions$datereported == "Alzheimers",])
total_samples <- nrow(predictions)

#This builds a dataframe for percentile of genetic risk score vs odds ratio
ordataframe<- predictions %>% mutate(decile = ntile(Alzheimers, 5)) %>% group_by(decile) %>% 
  summarize(Normal=summary(datereported)[["Normal"]], 
            Alzheimers=summary(datereported)[["Alzheimers"]], 
            OR=(summary(datereported)[["Alzheimers"]]/summary(datereported)[["Normal"]])/(num_patients/(total_samples-num_patients)),
            count = n())

#
# Compute 95% confidence intervals (TCI,BCI) top of ci and bottom of ci
# based on statperls https://www.ncbi.nlm.nih.gov/books/NBK431098/
ordataframe$tci<-exp(log(ordataframe$OR)+
                       1.96*sqrt(
                         (1/(ordataframe$Alzheimers+1))+
                           (1/(ordataframe$Normal+1))+
                           (1/sum(ordataframe$Alzheimers))+
                           (1/sum(ordataframe$Normal))
                       ))
ordataframe$bci<-exp(log(ordataframe$OR)-
                       1.96*sqrt(
                         (1/(ordataframe$Alzheimers+1))+
                           (1/(ordataframe$Normal+1))+
                           (1/sum(ordataframe$Alzheimers))+
                           (1/sum(ordataframe$Normal))
                       ))

ordataframe$xaxis<-sequence(5)

ggplot(ordataframe,aes(x=xaxis,y=OR)) + 
  geom_point() +
  geom_errorbar(aes(ymin=bci, ymax=tci), width=.2,position=position_dodge(0.05)) +
  ggtitle("Odds Ratio Between Quintiles of Predicted Alzheimer's Disease Patients") +
  xlab("Quintile") + ylab("Odds Ratio") +
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Graceful shutdown of cluster
h2o.removeAll()
h2o.shutdown(prompt = TRUE)


###################################################
# Create a Model using only 64 Split X Chromosome #
###################################################

# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_chrX_condensed_64splits.txt", sep = " ")

# AD data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")
my_data <- select(my_ukb_data,eid,
                  datereported = date_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0,
                  sourcereported = source_of_report_of_f00_dementia_in_alzheimers_disease_f130837_0_0)

# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")
all_data <- merge(all_data, my_data_age, by.x = "ids", by.y = "eid")

alzheimers <- all_data[!is.na(all_data[, "datereported"]),]
no_alzheimers_initial <- all_data[is.na(all_data[, "datereported"]),]

# Get breakdown of  patients by age
alzheimers_age <- table(alzheimers$yearBorn)

# Randomly get non disease patients for controls so that there is an equal amount based on age
# This will ensure that the controls are age-matched to the disease sample
# For example there are 5 patients born 1937 who have AD so we will randomly grab 5 other 
# patients born 1937 who do not have AD
no_alzheimers <- data.frame(matrix(ncol = ncol(no_alzheimers_initial), nrow = 0))
colnames(no_alzheimers) <- colnames(no_alzheimers_initial)
for (i in 1:length(alzheimers_age)) {
  temp <- alzheimers_age[i]
  age_check <- as.numeric(names(temp))
  number_cases <- as.numeric(unname(temp))
  possible_controls <- no_alzheimers_initial[no_alzheimers_initial$yearBorn == age_check,]
  no_alzheimers <- rbind(no_alzheimers, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}

alzheimers$datereported <- "Alzheimers"
no_alzheimers$datereported <- "Normal"

ind <- sample(c(TRUE, FALSE), nrow(alzheimers), replace=TRUE, prob=c(0.7, 0.3)) # Random split

train <- alzheimers[ind, ]
validate <- alzheimers[!ind, ]

controls <- no_alzheimers #  get controls

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
rm(alzheimers, no_alzheimers, controls, train_controls, validate_controls)
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
                    nfolds = 5,
                    max_runtime_secs = 300,
                    keep_cross_validation_predictions = TRUE)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# Generate explainability plots
exp_aml <- h2o.explain(model, validate.hex)
exp_aml
exp_leader <- h2o.explain(leader, validate.hex)
exp_leader
h2o.model_correlation_heatmap(model, validate.hex)
h2o.varimp_heatmap(model)
h2o.pd_multi_plot(model, validate.hex, "chrX_1.X1")
h2o.shap_summary_plot(leader, validate.hex)
h2o.ice_plot(leader, validate.hex, "chrX_1.X1")

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
mod_perf <- h2o.performance(leader,train=FALSE, xval=TRUE)


plot(mod_perf, type='roc', main=paste("Alzheimer's Disease 64-split X Chromosome CSLV Model Performance: AUC = ", auc))

# for example I have a list of H2OModels
list(leader) %>% 
  # map a function to each element in the list
  map(function(x) x %>% h2o.performance(valid=T) %>% 
        # from all these 'paths' in the object
        .@metrics %>% .$thresholds_and_metric_scores %>% 
        # extracting true positive rate and false positive rate
        .[c('tpr','fpr')] %>% 
        # add (0,0) and (1,1) for the start and end point of ROC curve
        add_row(tpr=0,fpr=0,.before=T) %>% 
        add_row(tpr=0,fpr=0,.before=F)) %>% 
  # add a column of model name for future grouping in ggplot2
  map2(c('64-Split'),
       function(x,y) x %>% add_column(model=y)) %>% 
  # reduce four data.frame to one
  reduce(rbind) %>% 
  # plot fpr and tpr, map model to color as grouping
  ggplot(aes(fpr,tpr,col=model))+
  geom_line()+
  geom_segment(aes(x=0,y=0,xend = 1, yend = 1),linetype = 2,col='grey')+
  xlab('False Positive Rate')+
  ylab('True Positive Rate')+
  ggtitle("ROC Curve for 64-Split X Chr Alzheimer's Disease Model")+
  annotate("text", x = .75, y = .25, label = paste("AUC =", auc)) + 
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()

# Plot Precision-Recall Graph
metrics <- as.data.frame(h2o.metric(mod_perf))
metrics %>%
  ggplot(aes(recall,precision)) + 
  geom_line() +
  theme_minimal() + ggtitle("Alzheimer's Disease 64-split X Chr CSLV Model Precision-Recall") +
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Find Odds ratio
#cross_val_preds <- gbm@model[["cross_validation_predictions"]]

cvpreds <- as.data.frame(h2o.getFrame(model@leader@model[["cross_validation_holdout_predictions_frame_id"]][["name"]]))

#combine with the original data
predictions <- cbind(cvpreds,train)

#rank by disease predicted score ascending
predictions$rank <- rank(predictions[,"datereported"])
predictions<- select(predictions, rank, everything())

#total number of AD cases
num_patients <- nrow(predictions[predictions$datereported == "Alzheimers",])
total_samples <- nrow(predictions)

#This builds a dataframe for percentile of genetic risk score vs odds ratio
ordataframe<- predictions %>% mutate(decile = ntile(Alzheimers, 5)) %>% group_by(decile) %>% 
  summarize(Normal=summary(datereported)[["Normal"]], 
            Alzheimers=summary(datereported)[["Alzheimers"]], 
            OR=(summary(datereported)[["Alzheimers"]]/summary(datereported)[["Normal"]])/(num_patients/(total_samples-num_patients)),
            count = n())

#
# Compute 95% confidence intervals (TCI,BCI) top of ci and bottom of ci
# based on statperls https://www.ncbi.nlm.nih.gov/books/NBK431098/
ordataframe$tci<-exp(log(ordataframe$OR)+
                       1.96*sqrt(
                         (1/(ordataframe$Alzheimers+1))+
                           (1/(ordataframe$Normal+1))+
                           (1/sum(ordataframe$Alzheimers))+
                           (1/sum(ordataframe$Normal))
                       ))
ordataframe$bci<-exp(log(ordataframe$OR)-
                       1.96*sqrt(
                         (1/(ordataframe$Alzheimers+1))+
                           (1/(ordataframe$Normal+1))+
                           (1/sum(ordataframe$Alzheimers))+
                           (1/sum(ordataframe$Normal))
                       ))

ordataframe$xaxis<-sequence(5)

ggplot(ordataframe,aes(x=xaxis,y=OR)) + 
  geom_point() +
  geom_errorbar(aes(ymin=bci, ymax=tci), width=.2,position=position_dodge(0.05)) +
  ggtitle("Odds Ratio Between Quintiles of Predicted Alzheimer's Disease Patients") +
  xlab("Quintile") + ylab("Odds Ratio") +
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Graceful shutdown of cluster
h2o.removeAll()
h2o.shutdown(prompt = TRUE)

###################################################
# Calculate Confidence Interval for 1 Split Model #
###################################################

# Preallocate vector for aucs
rm(results, predictions, model_types)
rm(results4, predictions4, model_types4)
rm(results8, predictions8, model_types8)
results <- c()
predictions <- c()
model_types <- c()
numModels <- 150
maxRuntime <- 60 # This is in seconds

# Run 100 expirements or train 100 Auto ML models using randomized set of training data each time
# Each model will also have 5 fold cross-validation as a base parameter.

####
# This section is for the 1 split model
###
# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
# Get CNV Scale Data
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed.txt", sep = " ")

# AD data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")
my_data <- select(my_ukb_data,eid,
                  datereported = date_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0,
                  sourcereported = source_of_report_of_f00_dementia_in_alzheimers_disease_f130837_0_0)

# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")
all_data <- merge(all_data, my_data_age, by.x = "ids", by.y = "eid")

alzheimers <- all_data[!is.na(all_data[, "datereported"]),]
no_alzheimers_initial <- all_data[is.na(all_data[, "datereported"]),]

# Get breakdown of  patients by age
alzheimers_age <- table(alzheimers$yearBorn)

for (i in 1:numModels) {
  # Randomly get non disease patients for controls so that there is an equal amount based on age
  # This will ensure that the controls are age-matched to the disease sample
  # For example there are 5 patients born 1937 who have AD so we will randomly grab 5 other 
  # patients born 1937 who do not have AD
  no_alzheimers <- data.frame(matrix(ncol = ncol(no_alzheimers_initial), nrow = 0))
  colnames(no_alzheimers) <- colnames(no_alzheimers_initial)
  for (i in 1:length(alzheimers_age)) {
    temp <- alzheimers_age[i]
    age_check <- as.numeric(names(temp))
    number_cases <- as.numeric(unname(temp))
    possible_controls <- no_alzheimers_initial %>% filter(yearBorn == age_check)
    no_alzheimers <- rbind(no_alzheimers, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
  }
  
  alzheimers$datereported <- TRUE
  no_alzheimers$datereported <- FALSE
  
  ind <- sample(c(TRUE, FALSE), nrow(alzheimers), replace=TRUE, prob=c(0.7, 0.3)) # Random split
  
  train <- alzheimers[ind, ]
  validate <- alzheimers[!ind, ]
  
  controls <- no_alzheimers #  get controls
  
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

mean(results)
sd(results)

confidence_interval(results, 0.90)
confidence_interval(results, 0.95)
confidence_interval(results, 0.99)

confidence_interval(predictions, 0.90)
confidence_interval(predictions, 0.95)
confidence_interval(predictions, 0.99)

# Graceful shutdown of cluster
h2o.removeAll()
h2o.shutdown(prompt = TRUE)


# Make plots
aucs <- data.frame(results)
models <- data.frame(model_types)
finalModelResults <- cbind(aucs, models)
ggplot(finalModelResults, aes(x=model_types, y=results, color=model_types)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.75) +
  scale_color_brewer(palette="Dark2") + 
  ggtitle("Comparison of Alzheimer's Disease Prediction AUCs by Model 1 Split") + 
  xlab("Model Algorithms") + ylab("AUC of Predictions") +
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()


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
plot(y ~ x, data = df, type = 'p', xlab = "Model Number", ylab = "AUC of Predictions")
# add fill
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', density=10,border = NA)
# model
abline(mod)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
labels()
title("Confidence Interval of AUCs of Alzheimer's Disease Models 1 Split")

###################################################
# Calculate Confidence Interval for 4 Split Model #
###################################################

# Preallocate vector for aucs
results4 <- c()
predictions4 <- c()
model_types4 <- c()
numModels <- 100
maxRuntime <- 60 # This is in seconds

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

# AD data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")
my_data <- select(my_ukb_data,eid,
                  datereported = date_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0,
                  sourcereported = source_of_report_of_f00_dementia_in_alzheimers_disease_f130837_0_0)

# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")
all_data <- merge(all_data, my_data_age, by.x = "ids", by.y = "eid")

alzheimers <- all_data[!is.na(all_data[, "datereported"]),]
no_alzheimers_initial <- all_data[is.na(all_data[, "datereported"]),]

# Get breakdown of  patients by age
alzheimers_age <- table(alzheimers$yearBorn)

for (i in 1:numModels) {
  # Randomly get non disease patients for controls so that there is an equal amount based on age
  # This will ensure that the controls are age-matched to the disease sample
  # For example there are 5 patients born 1937 who have AD so we will randomly grab 5 other 
  # patients born 1937 who do not have AD
  no_alzheimers <- data.frame(matrix(ncol = ncol(no_alzheimers_initial), nrow = 0))
  colnames(no_alzheimers) <- colnames(no_alzheimers_initial)
  for (i in 1:length(alzheimers_age)) {
    temp <- alzheimers_age[i]
    age_check <- as.numeric(names(temp))
    number_cases <- as.numeric(unname(temp))
    possible_controls <- no_alzheimers_initial %>% filter(yearBorn == age_check)
    no_alzheimers <- rbind(no_alzheimers, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
  }
  
  alzheimers$datereported <- TRUE
  no_alzheimers$datereported <- FALSE
  
  ind <- sample(c(TRUE, FALSE), nrow(alzheimers), replace=TRUE, prob=c(0.7, 0.3)) # Random split
  
  train <- alzheimers[ind, ]
  validate <- alzheimers[!ind, ]
  
  controls <- no_alzheimers #  get controls
  
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
  results4 <- c(results4, auc)
  model_types4 <- c(model_types4, leader@algorithm)
  
  # Attempt predict on validation frame
  prediction <- h2o.predict(object = leader, newdata = validate.hex)
  as.data.frame(prediction)
  summary(prediction, exact_quantiles = TRUE)
  
  validation.perf <- h2o.performance(leader, train = FALSE, xval=TRUE, newdata = validate.hex)
  validation.perf.auc <- validation.perf@metrics$AUC
  
  predictions4 <- c(predictions4, validation.perf.auc)
  h2o.removeAll()
  
  rm(train.hex, validate.hex, model, leader)
  
  # trigger removal of h2o back-end objects that got rm’d above, since the rm can be lazy.
  gc()
  # optional extra one to be paranoid.  this is usually very fast.
  gc()
  
  # optionally sanity check that you see only what you expect to see here, and not more.
  h2o.ls()
  
  # tell back-end cluster nodes to do three back-to-back JVM full GCs.
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
}

mean(results4)
sd(results4)

confidence_interval(results4, 0.90)
confidence_interval(results4, 0.95)
confidence_interval(results4, 0.99)

confidence_interval(predictions4, 0.90)
confidence_interval(predictions4, 0.95)
confidence_interval(predictions4, 0.99)

# Graceful shutdown of 
h2o.removeAll()
h2o.shutdown(prompt = TRUE)

# Make plots
aucs <- data.frame(results4)
models <- data.frame(model_types4)
finalModelResults <- cbind(aucs, models)
ggplot(finalModelResults, aes(x=model_types4, y=results, color=model_types4)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.75) +
  scale_color_brewer(palette="Dark2") + 
  ggtitle("Comparison of Alzheimer's Disease Prediction AUCs by Model 4 Splits") + 
  xlab("Model Algorithms") + ylab("AUC of Predictions") +
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()


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
plot(y ~ x, data = df, type = 'p', xlab = "Model Number", ylab = "AUC of Predictions")
# add fill
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', density=10,border = NA)
# model
abline(mod)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
labels()
title("Confidence Interval of AUCs of Alzheimer's Models 4 Splits")

###################################################
# Calculate Confidence Interval for 8 Split Model #
###################################################

# Preallocate vector for aucs
results8 <- c()
predictions8 <- c()
model_types8 <- c()
numModels <- 100
maxRuntime <- 60 # This is in seconds

# Run 100 expirements or train 100 Auto ML models using randomized set of training data each time
# Each model will also have 5 fold cross-validation as a base parameter.

####
# This section is for the 8 split model
###
# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
# Get CNV Scale Data
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed_8splits.txt", sep = " ")

# AD data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")
my_data <- select(my_ukb_data,eid,
                  datereported = date_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0,
                  sourcereported = source_of_report_of_f00_dementia_in_alzheimers_disease_f130837_0_0)

# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")
all_data <- merge(all_data, my_data_age, by.x = "ids", by.y = "eid")

alzheimers <- all_data[!is.na(all_data[, "datereported"]),]
no_alzheimers_initial <- all_data[is.na(all_data[, "datereported"]),]

# Get breakdown of  patients by age
alzheimers_age <- table(alzheimers$yearBorn)

for (i in 1:numModels) {
  # Randomly get non disease patients for controls so that there is an equal amount based on age
  # This will ensure that the controls are age-matched to the disease sample
  # For example there are 5 patients born 1937 who have AD so we will randomly grab 5 other 
  # patients born 1937 who do not have AD
  no_alzheimers <- data.frame(matrix(ncol = ncol(no_alzheimers_initial), nrow = 0))
  colnames(no_alzheimers) <- colnames(no_alzheimers_initial)
  for (i in 1:length(alzheimers_age)) {
    temp <- alzheimers_age[i]
    age_check <- as.numeric(names(temp))
    number_cases <- as.numeric(unname(temp))
    possible_controls <- no_alzheimers_initial %>% filter(yearBorn == age_check)
    no_alzheimers <- rbind(no_alzheimers, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
  }
  
  alzheimers$datereported <- TRUE
  no_alzheimers$datereported <- FALSE
  
  ind <- sample(c(TRUE, FALSE), nrow(alzheimers), replace=TRUE, prob=c(0.7, 0.3)) # Random split
  
  train <- alzheimers[ind, ]
  validate <- alzheimers[!ind, ]
  
  controls <- no_alzheimers #  get controls
  
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
  results8 <- c(results8, auc)
  model_types8 <- c(model_types8, leader@algorithm)
  
  
  # Attempt predict on validation frame
  prediction <- h2o.predict(object = leader, newdata = validate.hex)
  as.data.frame(prediction)
  summary(prediction, exact_quantiles = TRUE)
  
  validation.perf <- h2o.performance(leader, train = FALSE, xval=TRUE, newdata = validate.hex)
  validation.perf.auc <- validation.perf@metrics$AUC
  
  predictions8 <- c(predictions8, validation.perf.auc)
}

mean(results8)
sd(results8)

confidence_interval(results8, 0.90)
confidence_interval(results8, 0.95)
confidence_interval(results8, 0.99)

confidence_interval(predictions8, 0.90)
confidence_interval(predictions8, 0.95)
confidence_interval(predictions8, 0.99)

# Graceful shutdown of cluster
h2o.removeAll()
h2o.shutdown(prompt = TRUE)


# Make plots
aucs <- data.frame(results)
models <- data.frame(model_types8)
finalModelResults <- cbind(aucs, models)
ggplot(finalModelResults, aes(x=model_types8, y=results, color=model_types8)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.75) +
  scale_color_brewer(palette="Dark2") + 
  ggtitle("Comparison of Alzheimer's Disease Prediction AUCs by Model 8 Split") + 
  xlab("Model Algorithms") + ylab("AUC of Predictions") +
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()


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
plot(y ~ x, data = df, type = 'p', xlab = "Model Number", ylab = "AUC of Predictions")
# add fill
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', density=10,border = NA)
# model
abline(mod)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
labels()
title("Confidence Interval of AUCs of Alzheimer's Disease Models 8 Split")



###
# GRAPHS FOR 1, 4, 8 split results
###
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 12
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

split <- rep("1 split", length(results))
newRes <- round(results, digits = 3)
splitAucs <- cbind(newRes, split, model_types)

split <- rep("4 splits", length(results4))
newRes <- round(results4, digits = 3)
temp <- cbind(newRes, split, model_types4)
splitAucs <- rbind(splitAucs, temp)

split <- rep("8 splits", length(results8))
newRes <- round(results8, digits = 3)
temp <- cbind(newRes, split, model_types8)
splitAucs <- rbind(splitAucs, temp)

splitAucs <- data.frame(splitAucs)
ggplot(splitAucs,  aes(x=split, y=newRes, color=split, group = split)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.75) +
  #scale_color_brewer(palette="Dark2") + 
  scale_fill_manual(values = mycolors) +
  scale_y_discrete(breaks = c(0.50, 0.51,0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64)) +
  ggtitle("Comparison of Alzheimer's Disease Prediction AUCs by All Models Performance by Split") + 
  xlab("CSLV Splits") + ylab("AUC of Predictions") +
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()

# You can remove some model types if there were insufficient of those models made to make a relevant graph
ggplot(subset(splitAucs, split %in% c("1 split", "4 splits", "8 splits") & model_types %in% c("gbm", "glm", "stackedensemble", "deeplearning", "drf", "xgboost")),
       aes(x = split, y = newRes,  colour = interaction(model_types, split), group = split)) + facet_wrap( ~ model_types) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, outlier.alpha = 0.1) +
  #scale_color_brewer(palette="Dark2") + 
  scale_fill_manual(values = mycolors) +
  scale_y_discrete(breaks = c(0.50, 0.51,0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64)) +
  ggtitle("Comparison of Alzheimer's Disease Prediction AUCs by Model Type Performance by Split") + 
  xlab("Model Types") + ylab("AUC of Predictions") +
  theme_bw() + theme(plot.title = element_text(size = 18, face = "bold")) +
  scale_color_npg()

# Pring statistics
mean(results)
sd(results)
mean(results4)
sd(results4)
mean(results8)
sd(results8)

test <- t.test(results, results4, paired = TRUE, alternative = "two.sided")
test$p.value
test2 <- t.test(results, results8, paired = TRUE, alternative = "two.sided")
test2$p.value
test3 <- t.test(results4, results8, paired = TRUE, alternative = "two.sided")
test3$p.value


