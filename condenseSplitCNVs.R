library(parallel)
library(data.table)
library(ukbtools)
library(tidyverse)

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
        if(count %% 1000 == 0) system(paste("echo 'Now processing row:",count,"for ", filepath, "'"))
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
        if (count == split_locs[current_split] - 1) {
            results[current_split] <- tmp
            reset <- TRUE
            current_split <- current_split + 1
        }
    
        count <- count + 1
    }
  
    close(con)
  
    return(results)
}

condenseChromosomeWithSplits <- function(chrNum, splits) {
    res <- processFileWithSplits(paste("/data/ukbiobank/ukb_l2r_chr", chrNum, "_v2.txt", sep=""), numSplits = splits)
  
    return(res)
}

writeResultsWithSplits <- function(chrNum, splits) {
    for (i in 1:splits) {
        write.table(save1[[chrNum]][i], paste("/data/ukbiobank/ukb_l2r_chr", chrNum, "_averaged_", splits, "splits_split",i,".txt", sep=""), 
            sep = " ", row.names = FALSE, col.names = FALSE, dec = ".")
    }
  
}



# Run
numCores <- detectCores()
system.time(save1 <- mclapply(1:22, condenseChromosomeWithSplits, 
                              mc.cores = numCores - 1, splits = 4))
system.time(saveX <- mclapply("X", condenseChromosomeWithSplits, 
                              mc.cores = numCores - 1, splits = 4))

# Write to file
# Write results to file
system.time(mclapply(1:22, writeResultsWithSplits, 
                     mc.cores = numCores - 1, splits = 4))
for (i in 1:4) {
    system.time(write.table(saveX[[1]][i], paste("/data/ukbiobank/ukb_l2r_chrX_averaged_4splits_split",i,".txt", sep=""), 
                sep = " ", row.names = FALSE, col.names = FALSE, dec = "."))
}
