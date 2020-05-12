library(h2o)
library(h2o4gpu)
library(ukbtools)
library(tidyverse)

# Load h2o
h2o.init(nthreads=15)

## Create training and validation frames
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_allchr_condensed_4splits.txt", sep = " ")

# Get cancer and age related information
my_ukb_data_cancer <- ukb_df("ukb29274", path = "/data/ukbiobank/cancer")
my_data <- select(my_ukb_data_cancer, 
                      eid, 
                      yearBorn = year_of_birth_f34_0_0, 
                      sex_f31_0_0,
                      cancer_code_selfreported_f20001_0_0,
                      type_of_cancer_icd9_f40013_0_0,
                      type_of_cancer_icd9_f40013_1_0,
                      behaviour_of_cancer_tumour_f40012_0_0)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")

# Get control data
controldata <- all_data %>%  
    filter(sex_f31_0_0==0) %>% 
    filter(is.na(type_of_cancer_icd9_f40013_0_0) &
               is.na(type_of_cancer_icd9_f40013_1_0) &
               is.na(cancer_code_selfreported_f20001_0_0) &
               is.na(behaviour_of_cancer_tumour_f40012_0_0)
    ) %>%  #include only those cases that don't say anything that looks like cancer 
    select(ids,
           sex=sex_f31_0_0,
           yearBorn,
           cancer_selfreported=cancer_code_selfreported_f20001_0_0, 
           behavior=behaviour_of_cancer_tumour_f40012_0_0,
           icd9_0=type_of_cancer_icd9_f40013_0_0,
           icd9_1= type_of_cancer_icd9_f40013_1_0,
           contains("Chr")
    )

#
# This is the dataset of women with breast cancer.
# Note sex code = 0 is female, = 1 is male
#
cancerdata <- all_data %>%  
    filter(sex_f31_0_0==0 &
            (cancer_code_selfreported_f20001_0_0==1002 &
            (!is.na(type_of_cancer_icd9_f40013_0_0) |
            !is.na(type_of_cancer_icd9_f40013_1_0)))) %>% 
    select(ids,
           sex=sex_f31_0_0,
           yearBorn,
           cancer_selfreported=cancer_code_selfreported_f20001_0_0, 
           behavior=behaviour_of_cancer_tumour_f40012_0_0,
           icd9_0=type_of_cancer_icd9_f40013_0_0,
           icd9_1= type_of_cancer_icd9_f40013_1_0,
           contains("Chr")
    )

# Indicate cancers in data set
cancerdata$cancer_selfreported <- TRUE
controldata$cancer_selfreported <- FALSE

# Get breakdown of Schizophrenia patients by age
cancer_age <- table(cancerdata$yearBorn)

#Get non Bipolar Patients
no_cancer_initial <- controldata

# Randomly get non Bipolar patients for controls so that there is an equal amount based on age
# This will ensure that the controls are age-matched to the Bipolar sample
no_cancer <- data.frame(matrix(ncol = ncol(no_cancer_initial), nrow = 0))
colnames(no_cancer) <- colnames(no_cancer_initial)
for (i in 1:length(cancer_age)) {
    temp <- cancer_age[i]
    age_check <- as.numeric(names(temp))
    number_cases <- as.numeric(unname(temp))
    possible_controls <- no_cancer_initial[no_cancer_initial$yearBorn == age_check,]
    no_cancer <- rbind(no_cancer, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}

#Split for training and valdiation
ind <- sample(c(TRUE, FALSE), nrow(cancerdata), replace=TRUE, prob=c(0.7, 0.3)) # Random split
train <- cancerdata[ind, ]
validate <- cancerdata[!ind,]

controls <- no_cancer[sample(nrow(no_cancer), nrow(cancerdata), replace = TRUE), ] # Randomly get controls

train_controls <- controls[ind, ]
validate_controls <- controls[!ind, ]

# Combine controls with samples

train <- rbind(train, train_controls)
validate <- rbind(validate, validate_controls)

# Set response column to factor
train$cancer_selfreported <- as.factor(train$cancer_selfreported)
validate$cancer_selfreported <- as.factor(validate$cancer_selfreported)

#Remove unnecessary rows
train <- train[,!names(train) %in% c("ids")]
validate <- validate[,!names(validate) %in% c("ids")]

# Free up data 
rm(controldata, cancerdata, no_cancer_initial, controls, train_controls, validate_controls)
rm(my_data, my_ukb_data_cancer, all_data)

# Load data into h2o

train.hex <- as.h2o(train, destination_frame = "train.hex")  
validate.hex <- as.h2o(validate, destination_frame = "validate.hex")

#
#  I usually stop here and goto http://localhost:54321/flow/index.html 
#  h2o runs a local webserver on port 54321, it offers a nice little interface.
#  you can run the AutoML from the web browser there and it has some nice features so that 
# you can monitor the progress of the training.

#Response column
response <- "cancer_selfreported"
#Get Predictors
predictors <- colnames(train)
predictors <- predictors[! predictors %in% response] #Response cannot be a predictor
predictors <- predictors[! predictors %in% "yearBorn"] #Response cannot be a predictor
predictors <- predictors[! predictors %in% "sex"] #Response cannot be a predictor
predictors <- predictors[! predictors %in% "behavior"] #Response cannot be a predictor
predictors <- predictors[! predictors %in% "icd9_0"] #Response cannot be a predictor
predictors <- predictors[! predictors %in% "icd9_1"] #Response cannot be a predictor
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
plot(h2o.performance(leader,train=FALSE, xval=TRUE),type='roc',main=paste("Breast Cancer", auc))

# Print performance info of leader
leader@algorithm
h2o.performance(leader,train=FALSE, xval=TRUE)


# Graceful shutdown of cluster
h2o.shutdown(prompt = TRUE)
