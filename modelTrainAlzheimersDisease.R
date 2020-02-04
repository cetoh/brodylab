library(h2o)
library(h2o4gpu)
library(ukbtools)
library(tidyverse)

# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
condensed <- condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed.txt", sep = " ")
# Dementia Alzheimer's Disease data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")

my_data <- select(my_ukb_data,eid,
                  datereported = date_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0,
                  sourcereported = source_of_report_of_f00_dementia_in_alzheimers_disease_f130837_0_0)
# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")

alzheimers <- all_data[!is.na(all_data[, "sourcereported"]),]
no_alzheimers <- all_data[is.na(all_data[, "sourcereported"]),]

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
validate <- validate[,!names(train) %in% c("ids","datereported")]

# Free up data 
rm(no_alzheimers, alzheimers, controls, train_controls, validate_controls)
rm(my_data, my_ukb_data, condensed)

# Load data into h2o

train.hex <- as.h2o(train, destination_frame = "train.hex")  
validate.hex <- as.h2o(validate, destination_frame = "validate.hex")

#
#  I usually stop here and goto http://localhost:54321/flow/index.html 
#  h2o runs a local webserver on port 54321, it offers a nice little interface.
#  you can run the AutoML from the web browser there and it has some nice features so that 
# you can monitor the progress of the training.

train$sourcereported <- as.factor(train$sourcereported)

response <- "sourcereported"
predictors <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
       "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
       "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")
model <- h2o.automl(x = predictors,
                    y = response,
                    training_frame = train.hex,
                    validation_frame = validate.hex,
                    nfolds=5,
                    max_runtime_secs = 36000)

#record the Leading model AUC in the dataset
leader <- model@leader
auc=h2o.auc(leader, train=FALSE, xval=TRUE)

# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main="Alzheimer's Disease")

# Graceful shutdown of cluster
h2o.shutdown(prompt = TRUE)