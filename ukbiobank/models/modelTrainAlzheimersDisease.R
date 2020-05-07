library(h2o)
library(h2o4gpu)
library(ukbtools)
library(tidyverse)

# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed.txt", sep = " ")

# Dementia Alzheimer's Disease data
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

# Get Alzheimer's patients
alzheimers <- all_data[!is.na(all_data[, "sourcereported"]),]

# Get breakdown of Alzheimer's patients by age
alz_age <- table(alzheimers$yearBorn)

# As an age control we will only look at individuals born within the same time as the Alzheimer's patients for our non AD patients
all_data <- subset(all_data, all_data$yearBorn < max(alzheimers$yearBorn))
all_data <- subset(all_data, all_data$yearBorn > min(alzheimers$yearBorn))

#Get non Alzheimer's patients
no_alzheimers_initial <- all_data[is.na(all_data[, "sourcereported"]),]

# Randomly get non Alzheimer's patients for controls so that there is an equal amount based on age
# This will ensure that the controls are age-matched to the Alzheimer's Disease sample
# For example there are 5 patients born 1937 who have AD so we will randomly grab 5 other 
# patients born 1937 who do not have AD
no_alzheimers <- data.frame(matrix(ncol = ncol(no_alzheimers_initial), nrow = 0))
colnames(no_alzheimers) <- colnames(no_alzheimers_initial)
for (i in 1:length(alz_age)) {
  temp <- alz_age[i]
  age_check <- as.numeric(names(temp))
  number_cases <- as.numeric(unname(temp))
  possible_controls <- no_alzheimers_initial[no_alzheimers_initial$yearBorn == age_check,]
  no_alzheimers <- rbind(no_alzheimers, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}

ind <- sample(c(TRUE, FALSE), nrow(alzheimers), replace=TRUE, prob=c(0.7, 0.3)) # Random split

train <- alzheimers[ind, ]
validate <- alzheimers[!ind, ]

controls <- no_alzheimers[sample(nrow(no_alzheimers), nrow(alzheimers), replace = TRUE), ] # Randomly get controls

train_controls <- controls[ind, ]
validate_controls <- controls[!ind, ]

# Combine controls with samples

train <- rbind(train, train_controls)
validate <- rbind(validate, validate_controls)

# Set response column to boolean
train$sourcereported[!is.na(train$sourcereported)] <- TRUE
train$sourcereported[is.na(train$sourcereported)] <- FALSE
validate$sourcereported[!is.na(validate$sourcereported)] <- TRUE
validate$sourcereported[is.na(validate$sourcereported)] <- FALSE
train$sourcereported <- as.factor(train$sourcereported)
validate$sourcereported <- as.factor(validate$sourcereported)

#Remove unnecessary rows
train <- train[,!names(train) %in% c("ids","datereported")]
validate <- validate[,!names(train) %in% c("ids","datereported")]

# Free up data 
rm(no_alzheimers, alzheimers, controls, train_controls, validate_controls, no_alzheimers_initial)
rm(my_data, my_ukb_data, condensed, my_ukb_data_cancer, my_data_age)

# Load data into h2o

train.hex <- as.h2o(train, destination_frame = "train.hex")  
validate.hex <- as.h2o(validate, destination_frame = "validate.hex")

#
#  I usually stop here and goto http://localhost:54321/flow/index.html 
#  h2o runs a local webserver on port 54321, it offers a nice little interface.
#  you can run the AutoML from the web browser there and it has some nice features so that 
# you can monitor the progress of the training.

#Response column
response <- "sourcereported"
#Get Predictors
predictors <- colnames(train)
predictors <- predictors[! predictors %in% response] #Response cannot be a predictor
model <- h2o.automl(x = predictors,
                    y = response,
                    training_frame = train.hex,
                    validation_frame = validate.hex,
                    nfolds=5,
                    max_runtime_secs = 3600)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Alzheimer's Disease Full Condensed",auc))

# Print performanc info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Graceful shutdown of cluster
h2o.shutdown(prompt = TRUE)



####
# This section is for the 4 split model
###
# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed_4splits.txt", sep = " ")

# Dementia Alzheimer's Disease data
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

# Get Alzheimer's patients
alzheimers <- all_data[!is.na(all_data[, "sourcereported"]),]

# Get breakdown of Alzheimer's patients by age
alz_age <- table(alzheimers$yearBorn)

# As an age control we will only look at individuals born within the same time as the Alzheimer's patients for our non AD patients
all_data <- subset(all_data, all_data$yearBorn < max(alzheimers$yearBorn))
all_data <- subset(all_data, all_data$yearBorn > min(alzheimers$yearBorn))

#Get non Alzheimer's patients
no_alzheimers_initial <- all_data[is.na(all_data[, "sourcereported"]),]

# Randomly get non Alzheimer's patients for controls so that there is an equal amount based on age
# This will ensure that the controls are age-matched to the Alzheimer's Disease sample
# For example there are 5 patients born 1937 who have AD so we will randomly grab 5 other 
# patients born 1937 who do not have AD
no_alzheimers <- data.frame(matrix(ncol = ncol(no_alzheimers_initial), nrow = 0))
colnames(no_alzheimers) <- colnames(no_alzheimers_initial)
for (i in 1:length(alz_age)) {
  temp <- alz_age[i]
  age_check <- as.numeric(names(temp))
  number_cases <- as.numeric(unname(temp))
  possible_controls <- no_alzheimers_initial[no_alzheimers_initial$yearBorn == age_check,]
  no_alzheimers <- rbind(no_alzheimers, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}

combined_data <- rbind(alzheimers, no_alzheimers)

ind <- sample(c(TRUE, FALSE), nrow(alzheimers), replace=TRUE, prob=c(0.7, 0.3)) # Random split

train <- alzheimers[ind, ]
validate <- alzheimers[!ind, ]

controls <- no_alzheimers[sample(nrow(no_alzheimers), nrow(alzheimers), replace = TRUE), ] # Randomly get controls

train_controls <- controls[ind, ]
validate_controls <- controls[!ind, ]

# Combine controls with samples

train <- rbind(train, train_controls)
validate <- rbind(validate, validate_controls)

# Set response column to boolean
train$sourcereported[!is.na(train$sourcereported)] <- TRUE
train$sourcereported[is.na(train$sourcereported)] <- FALSE
validate$sourcereported[!is.na(validate$sourcereported)] <- TRUE
validate$sourcereported[is.na(validate$sourcereported)] <- FALSE
combined_data$sourcereported[!is.na(combined_data$sourcereported)] <- TRUE
combined_data$sourcereported[is.na(combined_data$sourcereported)] <- FALSE
train$sourcereported <- as.factor(train$sourcereported)
validate$sourcereported <- as.factor(validate$sourcereported)
combined_data$sourcereported <- as.factor(combined_data$sourcereported)

#Remove unnecessary rows
train <- train[,!names(train) %in% c("ids","datereported")]
validate <- validate[,!names(validate) %in% c("ids","datereported")]
combined_data <- combined_data[,!names(combined_data) %in% c("ids","datereported")]

# Free up data 
rm(no_alzheimers, alzheimers, controls, train_controls, validate_controls)
rm(my_data, my_ukb_data, condensed2, my_ukb_data_cancer, my_data_age)

# Load data into h2o

train.hex <- as.h2o(train, destination_frame = "train.hex")  
validate.hex <- as.h2o(validate, destination_frame = "validate.hex")
combined.hex <- as.h2o(combined_data, destination_frame = "combined.hex")

#
#  I usually stop here and goto http://localhost:54321/flow/index.html 
#  h2o runs a local webserver on port 54321, it offers a nice little interface.
#  you can run the AutoML from the web browser there and it has some nice features so that 
# you can monitor the progress of the training.

#Response column
response <- "sourcereported"
#Get Predictors
predictors <- colnames(train)
predictors <- predictors[! predictors %in% response] #Response cannot be a predictor

model <- h2o.automl(x = predictors,
                    y = response,
                    training_frame = combined.hex,
                    nfolds=5,
                    max_runtime_secs = 3600)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Alzheimer's Disease 4 Splits w/ yearBorn combined data",auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

#Response column
response <- "sourcereported"
#Get Predictors
predictors <- colnames(train)
predictors <- predictors[! predictors %in% response] #Response cannot be a predictor

model <- h2o.automl(x = predictors,
                    y = response,
                    training_frame = train.hex,
                    validation_frame = validate.hex,
                    nfolds=5,
                    max_runtime_secs = 3600)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Alzheimer's Disease 4 Splits w/ yearBorn",auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Attempt predict on validation frame
prediction <- h2o.predict(object = leader, newdata = validate.hex)
as.data.frame(prediction)
summary(prediction)

validation.perf <- h2o.performance(leader, train = FALSE, xval=TRUE, newdata = validate.hex)
validation.perf.auc <- validation.perf@metrics$AUC
plot(validation.perf, type = 'roc', main = paste("Alzheimer's Disease 4 Splits w/ yearBorn predictions", validation.perf.auc))


##
# Attempt without yearBorn as predictor
predictors <- predictors[! predictors %in% "yearBorn"] #Remove yearBorn predictor
model <- h2o.automl(x = predictors,
                    y = response,
                    training_frame = train.hex,
                    validation_frame = validate.hex,
                    nfolds=5,
                    max_runtime_secs = 3600)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Alzheimer's Disease 4 Splits w/o yearBorn",auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Attempt predict on validation frame
prediction <- h2o.predict(object = leader, newdata = validate.hex)
as.data.frame(prediction)
summary(prediction)

validation.perf <- h2o.performance(leader, train = FALSE, xval=TRUE, newdata = validate.hex)
validation.perf.auc <- validation.perf@metrics$AUC
plot(validation.perf, type = 'roc', main = paste("Alzheimer's Disease 4 Splits w/o yearBorn predictions", validation.perf.auc))

##
# Attempt with ONLY yearBorn as predictor
predictors <- "yearBorn"
model <- h2o.automl(x = predictors,
                    y = response,
                    training_frame = train.hex,
                    validation_frame = validate.hex,
                    nfolds=5,
                    max_runtime_secs = 3600)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Alzheimer's Disease yearBorn only",auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Graceful shutdown of cluster
h2o.shutdown(prompt = TRUE)

####
# This section is for the 8 split model
###
# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed_8splits.txt", sep = " ")

# Dementia Alzheimer's Disease data
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

all_data <- subset(all_data, all_data$yearBorn < 1945)

# Get Alzheimer's patients
alzheimers <- all_data[!is.na(all_data[, "sourcereported"]),]

# As an age control we will only look at individuals born within the same time as the Alzheimer's patients for our non AD patients
all_data <- subset(all_data, all_data$yearBorn < max(alzheimers$yearBorn))
all_data <- subset(all_data, all_data$yearBorn > min(alzheimers$yearBorn))

#Get non Alzheimer's patients
no_alzheimers_initial <- all_data[is.na(all_data[, "sourcereported"]),]

# Randomly get non Alzheimer's patients for controls so that there is an equal amount based on age
# This will ensure that the controls are age-matched to the Alzheimer's Disease sample
# For example there are 5 patients born 1937 who have AD so we will randomly grab 5 other 
# patients born 1937 who do not have AD
no_alzheimers <- data.frame(matrix(ncol = ncol(no_alzheimers_initial), nrow = 0))
colnames(no_alzheimers) <- colnames(no_alzheimers_initial)
for (i in 1:length(alz_age)) {
  temp <- alz_age[i]
  age_check <- as.numeric(names(temp))
  number_cases <- as.numeric(unname(temp))
  possible_controls <- no_alzheimers_initial[no_alzheimers_initial$yearBorn == age_check,]
  no_alzheimers <- rbind(no_alzheimers, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}

ind <- sample(c(TRUE, FALSE), nrow(alzheimers), replace=TRUE, prob=c(0.7, 0.3)) # Random split

train <- alzheimers[ind, ]
validate <- alzheimers[!ind, ]

controls <- no_alzheimers[sample(nrow(no_alzheimers), nrow(alzheimers)), ] # Randomly get controls

train_controls <- controls[ind, ]
validate_controls <- controls[!ind, ]

# Combine controls with samples

train <- rbind(train, train_controls)
validate <- rbind(validate, validate_controls)

# Set response column to boolean
train$sourcereported[!is.na(train$sourcereported)] <- TRUE
train$sourcereported[is.na(train$sourcereported)] <- FALSE
validate$sourcereported[!is.na(validate$sourcereported)] <- TRUE
validate$sourcereported[is.na(validate$sourcereported)] <- FALSE
train$sourcereported <- as.factor(train$sourcereported)
validate$sourcereported <- as.factor(validate$sourcereported)

#Remove unnecessary rows
train <- train[,!names(train) %in% c("ids","datereported")]
validate <- validate[,!names(validate) %in% c("ids","datereported")]

# Free up data 
rm(no_alzheimers, train_controls, validate_controls)
rm(my_data, my_ukb_data, condensed2, my_ukb_data_cancer, my_data_age)

# Load data into h2o

train.hex <- as.h2o(train, destination_frame = "train.hex")  
validate.hex <- as.h2o(validate, destination_frame = "validate.hex")

#
#  I usually stop here and goto http://localhost:54321/flow/index.html 
#  h2o runs a local webserver on port 54321, it offers a nice little interface.
#  you can run the AutoML from the web browser there and it has some nice features so that 
# you can monitor the progress of the training.

#Response column
response <- "sourcereported"
#Get Predictors
predictors <- colnames(train)
predictors <- predictors[! predictors %in% response] #Response cannot be a predictor

model <- h2o.automl(x = predictors,
                    y = response,
                    training_frame = train.hex,
                    validation_frame = validate.hex,
                    nfolds=5,
                    max_runtime_secs = 3600)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Alzheimer's Disease 8 Splits",auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

# Graceful shutdown of cluster
h2o.shutdown(prompt = TRUE)
