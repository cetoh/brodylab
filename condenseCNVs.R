library(parallel)
library(data.table)
library(ukbtools)
library(tidyverse)

#Read in all the IDs
ids <- read.csv("/data/ukbiobank/patientIDSall.fam", header = FALSE, sep = " ")
ids <- ids[[1]] #This just gets a vector of the ids in a column

#Functions
processFile <- function(filepath) {
    #Create a file connection. We do this as reading in the whole file would take up too much memory.
    con <- file(filepath, "r")
  
    # Read first line
    line <- scan(file = con, nlines = 1, quiet = TRUE)

    # Create processing data frame
    tmp <- data.frame(line)
  
    #For each line, average it with the previous line to keep the data frame to one column
     #At the end this column will represent the averaged CNV value for each patient
    while ( TRUE ) {
        line <- scan(file = con, nlines = 1, quiet = TRUE, na.strings = '-0.0m\xfdW\006\xab\xa4\xc4C\x84\xe5', fileEncoding = "UTF-8") # Get the line
    
        # If the line is empty you've reached the end of the file
        if ( length(line) == 0 ) { 
            break
        }
    
        # Compare and average values
        tmp[is.na(tmp)] <- -1 # NA values replaced by -1 by spec
        
        tmp <- cbind(tmp, line)
        tmp <- rowMeans(tmp, na.rm = TRUE)
    }
  
    close(con)
  
    return(tmp)
}

processFileWithSplits <- function(filepath, numSplits) {
    #Get number of lines
    system(paste("echo 'Get number of rows...'"))
    total_records <- as.integer(system2("wc",
                                    args = c("-l",
                                        filepath,
                                        " | awk '{print $1}'"),
                                    stdout = TRUE))
    system(paste("echo 'Rows:",total_records,"'"))
  
    # Calculate split size and split locations
    split_size <- floor(total_records/numSplits)
    split_locs <- numeric(numSplits)
    for (i in 0:numSplits) {
        split_locs[i] <- (i * split_size)
    }
    #Ensure last value is the total_record value to prevent data loss
    split_locs <- replace(split_locs, length(split_locs), total_records) 
  
    #Create a file connection. We do this as reading in the whole file would take up too much memory.
    con <- file(filepath, "r")
  
    # Read first line
    line <- scan(file = con, nlines = 1, quiet = TRUE)
  
    # Create processing data frame
    tmp <- data.frame(line)
    results <- data.frame(matrix(NA, nrow = length(line), ncol = numSplits))
  
    #For each line, average it with the previous line to keep the data frame to one column
    #At the end this column will represent the averaged CNV value for each patient
    count <- 1
    current_split <- 1
    reset <- FALSE
    while ( TRUE ) {
        if(count %% 1000 == 0) system(paste("echo 'now processing:",count,"'"))
        line <- scan(file = con, nlines = 1, quiet = TRUE, na.strings = '-0.0m\xfdW\006\xab\xa4\xc4C\x84\xe5', fileEncoding = "UTF-8") # Get the line
    
        # If the line is empty you've reached the end of the file
        if ( length(line) == 0 ) { 
            break
        }
    
        # Compare and average values
        if (reset) {
            tmp <- data.frame(line)
        } else {
            tmp <- cbind(tmp, line)
            reset <- FALSE
        }
    
        tmp[is.na(tmp)] <- -1 # NA values replaced by -1 by spec
        tmp <- rowMeans(tmp, na.rm = TRUE)
    
        # If a split is reached put the current averages into the results dataframe
        # Then set tmp frame to reset
        if (count == split_locs[current_split]) {
            results[current_split] <- tmp[1]
            reset <- TRUE
            current_split <- current_split + 1
        }
    
        count <- count + 1
    }
  
    close(con)
  
    return(results)
}

condenseChromosome <- function(chrNum) {
    res <- processFile(paste("/data/ukbiobank/ukb_l2r_chr", chrNum, "_v2.txt", sep=""))

    return(res)
}

condenseChromosomeWithSplits <- function(chrNum, splits) {
    res <- processFileWithSplits(paste("/data/ukbiobank/ukb_l2r_chr", chrNum, "_v2.txt", sep=""), numSplits = splits)
  
    return(res)
}

writeResults <- function(chrNum) {
    write.table(save1[[chrNum]], paste("/data/ukbiobank/ukb_l2r_chr", chrNum, "_averaged.txt", sep=""), 
                sep = " ", row.names = FALSE, col.names = FALSE, dec = ".")
}

##### Processing Actions #####
# Get Number of Cores available
numCores <- detectCores()

#For each chromosome average the CNVs for each patient
# We will now run functions in a parallized manner to utilize all available cores 
# except 1. This will do all 22 chromosomes in parallel. This still took 8 hours on our system.
system.time(save1 <- mclapply(1:22, condenseChromosome, mc.cores = numCores - 1))
system.time(chrX <- mclapply("X", condenseChromosome, mc.cores = numCores - 1))

# Write results to file
system.time(mclapply(1:22, writeResults, mc.cores = numCores - 1))
system.time(mclapply("X", writeResults, mc.cores = numCores - 1))

# Unify results into final data frame
condensed <- data.table(ids, chr1 = save1[[1]], chr2 = save1[[2]], 
                        chr3 = save1[[3]], chr4 = save1[[4]], 
                        chr5 = save1[[5]], chr6 = save1[[6]], 
                        chr7 = save1[[7]], chr8 = save1[[8]], 
                        chr9 = save1[[9]], chr10 = save1[[10]], 
                        chr11 = save1[[11]], chr12 = save1[[12]], 
                        chr13 = save1[[13]], chr14 = save1[[14]],
                        chr15 = save1[[15]], chr16 = save1[[16]], 
                        chr17 = save1[[17]], chr18 = save1[[18]], 
                        chr19 = save1[[19]], chr20 = save1[[20]], 
                        chr21 = save1[[21]], chr22 = save1[[22]], 
                        chrX = chrX[[1]])

# Free up memory
rm(chrX, save1)

# Write final table to file
write.table(condensed, "/data/ukbiobank/ukb_l2r_ids_allchr_condensed.txt",
            sep = " ", row.names = FALSE, col.names = TRUE, dec = ".")


# This section of code pulls out the relevant data for the actual health outcomes
# using the ukbtools package

# Breast Cancer data
my_ukb_data_cancer <- ukb_df("ukb29274", path="/data/ukbiobank/cancer")

my_data_cancer <- select(my_ukb_data_cancer,eid,sex=sex_f31_0_0,
                         selfreported = cancer_code_selfreported_f20001_0_0,
                         ICD9_0=type_of_cancer_icd9_f40013_0_0,
                         ICD9_1=type_of_cancer_icd9_f40013_1_0,
                         behvaior=behaviour_of_cancer_tumour_f40012_0_0)

# Free up memory
rm(my_ukb_data_cancer)

# Dementia Alzheimer's Disease data
my_ukb_data <- ukb_df("ukb39651", path="/data/ukbiobank")

my_data <- select(my_ukb_data,eid,
                  datereported = date_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0,
                  sourcereported = source_of_report_of_f00_dementia_in_alzheimers_disease_f130837_0_0)

# Free up memory
rm(my_ukb_data)

# Merge with CNV data
all_data <- merge(condensed, my_data, by.x = "ids", by.y = "eid")

# Free up memory
rm(my_data, my_data_cancer)