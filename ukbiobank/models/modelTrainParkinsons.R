library(h2o)
library(h2o4gpu)
library(ukbtools)
library(tidyverse)

# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed_4splits.txt", sep = " ")

# Parkinson's Disease data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")
my_data <- select(my_ukb_data,eid,
                  datereported = date_g20_first_reported_parkinsons_disease_f131022_0_0)

# Get age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data_age <- select(my_ukb_data_cancer, eid, yearBorn = year_of_birth_f34_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")
all_data <- merge(all_data, my_data_age, by.x = "ids", by.y = "eid")

# Get Parkinson's Patients
parkinsons <- all_data[!is.na(all_data[, "datereported"]),]

# As an age control we will only look at individuals born within the same time as 
# the Parkinson's patients for our non Parkinson's patients
all_data <- subset(all_data, all_data$yearBorn < max(parkinsons$yearBorn))
all_data <- subset(all_data, all_data$yearBorn > min(parkinsons$yearBorn))

no_parkinsons <- all_data[is.na(all_data[, "datereported"]),]

parkinsons$datereported <- TRUE
no_parkinsons$datereported <- FALSE

ind <- sample(c(TRUE, FALSE), nrow(parkinsons), replace=TRUE, prob=c(0.7, 0.3)) # Random split

train <- parkinsons[ind, ]
validate <- parkinsons[!ind, ]

controls <- no_parkinsons[sample(nrow(no_parkinsons), nrow(parkinsons)), ] # Randomly get controls

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
rm(no_parkinsons, parkinsons, controls, train_controls, validate_controls)
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
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Parkinson's Disease", auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)

#Try a deeplearning model
dl_model <- h2o.deeplearning(
    model_id = "dl_model_tuned", 
    training_frame = train.hex, 
    validation_frame = validate.hex, 
    x = predictors, 
    y = response, 
    #variable_importances = T,
    hidden = c(200,200),          ## more hidden layers -> more complex interactions
    epochs = 100000,
    stopping_rounds=5,
    stopping_tolerance=0.01,
    stopping_metric = "misclassification"
) 

summary(dl_model)
# Graceful shutdown of cluster
h2o.shutdown(prompt = TRUE)