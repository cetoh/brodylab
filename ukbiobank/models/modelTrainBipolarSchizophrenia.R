library(h2o)
library(h2o4gpu)
library(ukbtools)
library(tidyverse)

# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed_4splits.txt", sep = " ")

# Bipolar and Schizophrenia data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")
my_data <- select(my_ukb_data,eid,
                  datereported = date_f31_first_reported_bipolar_affective_disorder_f130892_0_0,
                  sourcereported = source_of_report_of_f31_bipolar_affective_disorder_f130893_0_0,
                  datereported2 = date_f20_first_reported_schizophrenia_f130874_0_0,
                  sourcereported2 = source_of_report_of_f20_schizophrenia_f130875_0_0)

# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")
all_data <- merge(all_data, my_data_age, by.x = "ids", by.y = "eid")

bipolar <- all_data[!is.na(all_data[, "datereported"]),]
no_bipolar <- all_data[is.na(all_data[, "datereported"]),]

schiz <- all_data[!is.na(all_data[, "datereported2"]),]
no_schiz <- all_data[is.na(all_data[, "datereported2"]),]

#Get non Bipolar Patients
no_bipolar_initial <- all_data[is.na(all_data[, "datereported"]),]

# Randomly get non Bipolar patients for controls so that there is an equal amount based on age
# This will ensure that the controls are age-matched to the Bipolar sample
no_bipolar <- data.frame(matrix(ncol = ncol(no_bipolar_initial), nrow = 0))
colnames(no_bipolar) <- colnames(no_bipolar_initial)
for (i in 1:length(alz_age)) {
  temp <- alz_age[i]
  age_check <- as.numeric(names(temp))
  number_cases <- as.numeric(unname(temp))
  possible_controls <- no_bipolar_initial[no_bipolar_initial$yearBorn == age_check,]
  no_bipolar <- rbind(no_bipolar, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}

#Get non Bipolar Patients
no_schiz_initial <- all_data[is.na(all_data[, "datereported2"]),]

# Randomly get non Bipolar patients for controls so that there is an equal amount based on age
# This will ensure that the controls are age-matched to the Bipolar sample
no_schiz <- data.frame(matrix(ncol = ncol(no_schiz_initial), nrow = 0))
colnames(no_schiz) <- colnames(no_schiz_initial)
for (i in 1:length(alz_age)) {
  temp <- alz_age[i]
  age_check <- as.numeric(names(temp))
  number_cases <- as.numeric(unname(temp))
  possible_controls <- no_schiz_initial[no_schiz_initial$yearBorn == age_check,]
  no_schiz <- rbind(no_schiz, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}

#Set Response To Boolean
bipolar$datereported <- TRUE
no_bipolar$datereported <- FALSE

schiz$datereported2 <- TRUE
no_schiz$datereported2 <- FALSE

ind <- sample(c(TRUE, FALSE), nrow(bipolar), replace = TRUE, prob=c(0.7, 0.3)) # Random split
ind <- sample(c(TRUE, FALSE), nrow(schiz), replace = TRUE, prob=c(0.7, 0.3)) # Random split

train <- bipolar[ind, ]
validate <- bipolar[!ind, ]
train2 <- schiz[ind, ]
validate2 <- schiz[!ind, ]

controls <- no_bipolar[sample(nrow(no_bipolar), nrow(bipolar), replace = TRUE), ] # Randomly get controls
controls2 <- no_schiz[sample(nrow(no_schiz), nrow(schiz), replace = TRUE), ] # Randomly get controls

train_controls <- controls[ind, ]
validate_controls <- controls[!ind, ]
train_controls2 <- controls2[ind, ]
validate_controls2 <- controls2[!ind, ]

# Combine controls with samples

train <- rbind(train, train_controls)
validate <- rbind(validate, validate_controls)
train2 <- rbind(train2, train_controls2)
validate2 <- rbind(validate2, validate_controls2)

# Set response column to factor
train$datereported <- as.factor(train$datereported)
validate$datereported <- as.factor(validate$datereported)
train2$datereported2 <- as.factor(train2$datereported2)
validate2$datereported2 <- as.factor(validate2$datereported2)

#Remove unnecessary columns
train <- train[,!names(train) %in% c("ids", "sex", "behavior")]
validate <- validate[,!names(validate) %in% c("ids", "sex", "behavior")]
train2 <- train2[,!names(train2) %in% c("ids", "sex", "behavior")]
validate2 <- validate2[,!names(validate2) %in% c("ids", "sex", "behavior")]

# Free up data 
rm(no_bipolar, bipolar, schiz, no_schiz, controls, controls2, train_controls, validate_controls, train_controls2, validate_controls2)
rm(my_data, my_ukb_data, my_ukb_data_cancer, my_data_age, possible_controls, no_bipolar_initial, no_schiz_initial)

# Load data into h2o

train.hex <- as.h2o(train, destination_frame = "train.hex")  
validate.hex <- as.h2o(validate, destination_frame = "validate.hex")
train2.hex <- as.h2o(train2, destination_frame = "train2.hex")  
validate2.hex <- as.h2o(validate2, destination_frame = "validate2.hex")

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
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Bipolar Disease", auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

response2 <- "datereported2"

model <- h2o.automl(x = predictors,
                    y = response2,
                    training_frame = train2.hex,
                    validation_frame = validate2.hex,
                    nfolds=5,
                    max_runtime_secs = 3600)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Schizophrenia", auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)


# Graceful shutdown of cluster
h2o.shutdown(prompt = TRUE)
