
#The following two commands remove any previously installed H2O packages for R.
if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
                 
#Next, we download packages that H2O depends on.

if (! ("methods" %in% rownames(installed.packages()))) { install.packages("methods") }
if (! ("statmod" %in% rownames(installed.packages()))) { install.packages("statmod") }
if (! ("stats" %in% rownames(installed.packages()))) { install.packages("stats") }
if (! ("graphics" %in% rownames(installed.packages()))) { install.packages("graphics") }
if (! ("RCurl" %in% rownames(installed.packages()))) { install.packages("RCurl") }
if (! ("rjson" %in% rownames(installed.packages()))) { install.packages("rjson") }
if (! ("tools" %in% rownames(installed.packages()))) { install.packages("tools") }
if (! ("utils" %in% rownames(installed.packages()))) { install.packages("utils") }

#Now we download, install and initialize the H2O package for R.
# Must use h2o v3.32.0.2 or higher for the explainability plots
# Run this section to install latest h2o package
if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
pkgs <- c("RCurl","jsonlite","ukbtools")
for (pkg in pkgs) {
    if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
}
install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R")))

library(h2o)
library(ukbtools)
library('stringr')

#Start a single-node instance of H2O using all available processor cores and reserve 8GB of memory (less should work too, check your logs for `WARN: Pausing to swap to disk; more memory may help` using the Flow GUI at [localhost:54321](http://localhost:54321), `Admin` -> `View Log` -> `SELECT LOG FILE TYPE: warn`).

h2oServer <- h2o.init(ip="localhost", port=54321, max_mem_size="8g", nthreads=-1)
#h2oServer <- h2o.init(ip="h2o-cluster", port=54321) # optional: connect to running H2O cluster

score_test_set=T  #disable if only interested in training throughput
run <- function(extra_params) {
    str(extra_params)
    print("Training.")
    model <- do.call(h2o.deeplearning, modifyList(list(x=predictors, y=response,
                                                       training_frame=train.hex, model_id="dlmodel"), extra_params))
    
    sampleshist <- model@model$scoring_history$samples
    samples <- sampleshist[length(sampleshist)]
    time <- model@model$run_time/1000
    auc <- h2o.auc(model)
    print(paste0("training samples: ", samples))
    print(paste0("training time   : ", time, " seconds"))
    print(paste0("training speed  : ", samples/time, " samples/second"))
    print(paste0("auc :",auc))
    
    if (score_test_set) {
        print("Scoring on test set.")
        ## Note: This scores full test set (10,000 rows) - can take time!
        p <- h2o.performance(model, validate.hex)
        cm <- h2o.confusionMatrix(p)
        test_error <- cm$Error[length(cm$Error)]
        print(paste0("test set error  : ", test_error))
    } else {
        test_error <- 1.0
    }
    h2o.rm("dlmodel")
    c(paste(names(extra_params), extra_params, sep = "=", collapse=" "), 
      samples, sprintf("%.3f", time), 
      sprintf("%.3f", samples/time), sprintf("%.3f", test_error), auc)
}
writecsv <- function(results, file) {
    table <- matrix(unlist(results), ncol = 6, byrow = TRUE)
    colnames(table) <- c("parameters", "training samples",
                         "training time", "training speed", "test set error", "auc")
    write.csv(table, file.path(workdir,file), 
              col.names = T, row.names=F, quote=T, sep=",")
}
workdir="/data/ukbiobank"
setwd(workdir)
###################################################
# Create a Model using only 64 Split X Chromosome #
###################################################

## Create training and validation frames
condensed <- read.csv("/data/ukbiobank/ukb_l2r_ids_chrX_condensed_64splits.txt", sep = " ")

# Schizophrenia data
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
# For example there are 5 patients born 1937 who have the disease so we will randomly grab 5 other 
# patients born 1937 who do not have the disease
no_schiz <- data.frame(matrix(ncol = ncol(no_schiz_initial), nrow = 0))
colnames(no_schiz) <- colnames(no_schiz_initial)
for (i in 1:length(schiz_age)) {
    temp <- schiz_age[i]
    age_check <- as.numeric(names(temp))
    number_cases <- as.numeric(unname(temp))
    possible_controls <- no_schiz_initial[no_schiz_initial$yearBorn == age_check,]
    no_schiz <- rbind(no_schiz, possible_controls[sample(nrow(possible_controls), number_cases, replace = TRUE), ])
}

schiz$datereported <- "Schizophrenia"
no_schiz$datereported <- "Normal"

ind <- sample(c(TRUE, FALSE), nrow(schiz), replace=TRUE, prob=c(0.8, 0.2)) # Random split

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

### First Study - Various Neural Network Topologies
#As a first study, we vary the network topology, and we also fix a number of training samples, specified by the number of epochs (we process approximately 10% of the dataset here, chosen at random). All other parameters are left at default values (per-weight adaptive learning rate, no L1/L2 regularization, no Dropout, Rectifier activation function). Note that we run shallow neural nets (1 hidden layer), and deep neural nets with 2 or 3 hidden layers. You can modify this and run other parameter combinations as well (and the same holds for all experiments here).

#The number of columns of the dataset (all numerical) translates directly into the size of the first layer of neurons (input layer), and hence significantly affects the size of the model. The number of output neurons is 10, one for each class (digit) probability. For the MNIST dataset, all Deep Learning models below will have 717 input neurons, as the other 67 pixel values are constant (white background) and thus ignored. The size of the model is given by the number of connections between the fully connected input+hidden+output neuron layers (and their bias values), and grows as the sum of the product of the number of connected neurons (~quadratically). Since training involves reading and updating all the model coefficients (weights and biases), the model complexity is linear in the size of the model, up to memory hierarchy effects in the x86 hardware. The speed of the memory subsystem (both in terms of latency and bandwidth) has a direct impact on training speed.

EPOCHS=.1 #increase if you have the patience (or run multi-node), shouldn't matter much for single-node
args <- list(
    list(hidden=c(64),             epochs=EPOCHS),
    list(hidden=c(128),            epochs=EPOCHS),
    list(hidden=c(256),            epochs=EPOCHS),
    list(hidden=c(512),            epochs=EPOCHS),
    list(hidden=c(1024),           epochs=EPOCHS),
    list(hidden=c(64,64),          epochs=EPOCHS),
    list(hidden=c(128,128),        epochs=EPOCHS),
    list(hidden=c(256,256),        epochs=EPOCHS),
    list(hidden=c(512,512),        epochs=EPOCHS),
    list(hidden=c(1024,1024),      epochs=EPOCHS),
    list(hidden=c(64,64,64),       epochs=EPOCHS),
    list(hidden=c(128,128,128),    epochs=EPOCHS),
    list(hidden=c(256,256,256),    epochs=EPOCHS),
    list(hidden=c(512,512,512),    epochs=EPOCHS),
    list(hidden=c(1024,1024,1024), epochs=EPOCHS)
)
writecsv(lapply(args, run), "network_topology.csv")

network_topology <- read.csv("/data/ukbiobank/network_topology.csv")
#We can plot the training speed (x-axis, more is faster) for the runs above, in the same order as the listing above. On the y-axis, we denote the overall training time in seconds.
plot(network_topology$training.speed, network_topology$training.time)

### Scoring Overhead
#During training, the model is getting scored on both the training and optional validation datasets at user-given intervals, with a user-given duty cycle (fraction of training time) and on user-specified subsamples of the datasets. Sampling can help reduce this scoring overhead significantly, especially for large network topologies. Often, a statistical scoring uncertainty of 1% is good enough to get a sense of model performance, and 10,000 samples are typically fine for that purpose. For multi-class problems, `score_validation_sampling` can be set to stratified sampling to maintain a reasonable class representation.

#For small networks and small datasets, the scoring time is usually negligible. The training data set is sub-sampled to `score_training_samples` (default: 10,000) rows (for use in scoring only), while a user-given validation dataset is used in full, without sampling (`score_validation_samples=0`) unless specified otherwise. This is important to understand, as it can lead to a large scoring overhead, especially for short overall training durations (no matter how few samples are trained), as scoring happens at least once (every `score_interval` (default: 5) seconds after the first MapReduce pass is complete, and at least once at the end of training), but no more than a `duty_cycle` (default: 10%) fraction of the total runtime.

#To illustrate this, we train the same model 7 times with different scoring selections.
args <- list(
    list(hidden=c(1024, 1024), epochs=EPOCHS, validation_frame=test_hex),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=60000, 
         score_duty_cycle=1, score_interval=1),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=60000),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=10000),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=1000),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100, 
         score_duty_cycle=0, score_interval=10000)
)
writecsv(lapply(args, run), "scoring_overhead.csv")

scoring_overhead <- read.csv("/data/ukbiobank/scoring_overhead.csv")
#As you can see, the overhead for scoring the entire training dataset at every opportunity (after each MapReduce pass, more on this below) can be significant (it's the same as forward propagation!). The default option (first run) with a specified (reasonably sized) validation dataset is a good compromise, and it's further possible to reduce scoring to a minimum. Test set classification errors are computed after training is finished, and are expected to be similar here (same model parameters) up to noise due to different initial conditions and (intentional) multi-threading race conditions (see Hogwild! reference above).
plot(scoring_overhead$training.samples, scoring_overhead$training.time)

### Adaptive vs Manual Learning Rate and Momentum
#Next, we compare adaptive (per-coefficient) learning rate versus manually specifying learning rate and optional momentum.
args <- list(
    list(hidden=c(1024, 1024), epochs=EPOCHS, 
         score_training_samples=100, adaptive_rate=T),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS,
         score_training_samples=100, adaptive_rate=T, 
         rho=0.95, epsilon=1e-6),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS,
         score_training_samples=100, adaptive_rate=F),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS,
         score_training_samples=100, adaptive_rate=F, 
         rate=1e-3, momentum_start=0.5, momentum_ramp=1e5, momentum_stable=0.99)
)
writecsv(lapply(args, run), "adaptive_rate.csv")

#It is clear that the fastest option is to run with a manual learning rate and momentum ramp

#Note that the lowest test set error was achieved these settings as well
adaptive_rate <- read.csv("/data/ukbiobank/adaptive_rate.csv")
plot(adaptive_rate$training.samples, adaptive_rate$training.time)

### Training Samples Per Iteration

#As explained earlier, model scoring (and multi-node model averaging) happens after every MapReduce step, the length of which is given in terms of the number of training samples via the parameters `train_samples_per_iteration`. It can be set to any value, -2, -1, 0, 1, 2, ... and as high as desired. Special values are -2 (auto-tuning), -1 (all available data on each node) and 0 (one epoch). Values of 1 and above specify the (approximate) number of samples to process per MapReduce iteration. If the value is less than the number of training data points per node, stochastic sampling without replacement is done on the node-local data (which can be the entire dataset if `replicate_training_data=T`). If the value is larger than the number of training data points per node, stochastic sampling with replacement is done on the node-local data. In either case, the number of actually processed samples can vary between runs due to the stochastic nature of sampling. In single-node operation, there is potentially some model scoring overhead in-between these MapReduce iterations (see above for dependency on parameters), but no model-averaging overhead as there is only one model per node (racily updated by multiple threads).
args <- list(
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100, 
         score_interval=1, train_samples_per_iteration=-2),
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100,
         score_interval=1, train_samples_per_iteration=100),  
    
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100,
         score_interval=1, train_samples_per_iteration=1000), 
    
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100,
         score_interval=1, train_samples_per_iteration=6000)  
)
writecsv(lapply(args, run), "train_samples_per_iteration.csv")

train_samples_per_iteration <- read.csv("/data/ukbiobank/train_samples_per_iteration.csv")
plot(train_samples_per_iteration$training.samples, train_samples_per_iteration$training.time)

### Different activation functions

#The activation function can make a big difference in runtime, as evidenced below. The reason for this is that H2O Deep Learning is optimized to take advantage of sparse activation of the neurons, and the Rectifier activation function (`max(0,x)`) is sparse and leads to less computation (at the expense of some accuracy, for which the speed often more than compensates for). Dropout (default: 0.5, randomly disable half of the hidden neuron activations for each training row) further increases the sparsity and makes training faster, and for noisy high-dimensional data, often also leads to better generalization (lower test set error) as it is a form of ensembling of many models with different network architectures. Since Dropout is a choice the user needs to make (at least for now), our default `Rectifier` activation function seems like a reasonable choice, especially since it also leads to a low test set error. We note that we omitted the `Maxout` activation here as we find it least useful (might require a revisit).

args <- list(
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100,
         train_samples_per_iteration=1000, activation="Rectifier"),
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100,
         train_samples_per_iteration=1000, activation="RectifierWithDropout"),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100, 
         train_samples_per_iteration=1000, activation="Tanh"),
    
    list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100, 
         train_samples_per_iteration=1000, activation="TanhWithDropout")
    
    #  list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100, 
    #       train_samples_per_iteration=1000, activation="Maxout"),
    
    #  list(hidden=c(1024, 1024), epochs=EPOCHS, score_training_samples=100, 
    #       train_samples_per_iteration=1000, activation="MaxoutWithDropout")
)
writecsv(lapply(args, run), "activation_function.csv")
activation_function <- read.csv("/data/ukbiobank/activation_function.csv")
plot(activation_function$training.speed, activation_function$test.set.error)

### External Large Deep Network Benchmark
#Recently, we learned of a [Deep Learning project that benchmarked against H2O](http://www.comp.nus.edu.sg/~dbsystem/singa/development/2015/01/29/compare-h2o/), and we tried to reproduce their results, especially since they claimed that H2O is slower.

#Since we know that [all the compute cores are 100% max'ed out with H2O Deep Learning](https://twitter.com/arnocandel/status/499715893505454080), or, as my colleagues call it, *They got Arno'd*... I had to see what's going on. And trust me, I have profiled the living daylight out of H2O Deep Learning, and added some [low-level optimizations](https://github.com/h2oai/h2o-3/blob/master/h2o-algos/src/main/java/hex/deeplearning/Neurons.java#L1070-L1109) for maximum throughput.

#For a direct comparison, we aim to use the same [H2O parameters](http://www.comp.nus.edu.sg/~dinhtta/MLP.java) as the external H2O Benchmark study (which was using Java and Scala code instead of R, but that doesn't make a difference in runtime, as the same H2O Java code is executed):
#Notes:
#    *  Unfortunately, there was no attempt made at comparing test set errors in the benchmark study. Throughput alone is not a good measure of model performance, see part 2 below for test set error comparisons and our model with world-record test set accuracy on the MNIST dataset.
#    * From the parameters used by the benchmark study, we can see that the Java field `_convert_to_enum` was left at its default value of `false`, leading to a regression model (1 output neuron with `MeanSquare` loss function), which is arguably not the best way to do digit classification. We will set the corresponding parameter in R with `do_classification=T`.
#    * We also note that the last *hidden* layer was specified to have 10 neurons, probably with the intention to emulate the 10-class probability distribution in the last hidden layer of 10 deep features. However, in H2O, the output layer is automatically added after the last hidden layer, and the user only specifies the `hidden` layers. For correctness, we will keep the first 5 hidden layers, and let H2O Deep Learning add the output layer with 10 neurons (and the `CrossEntropy` loss function) automatically (since we're also setting `do_classification=T`). We note that the difference in model size (and loss function) is a mere 11 numbers, and has no impact on training speed (final layers 500->10->1 vs 500->10).
#    * For unknown reasons, the benchmarkers *instrumented* the H2O .jar file to measure speed inside the mappers/reducers. This is neither necessary nor very useful, as the logs already print out the training speed based on pure wall clock time, which is easily verified.
#    * H2O Deep Learning currently only does online learning (mini-batch size of 1), while it appears that the other method used a [mini-batch size of 1000](http://www.comp.nus.edu.sg/~dbsystem/singa/examples/2015/01/20/deep-simple-mlp/). This algorithmic difference both affects computational performance and model accuracy, and makes comparisons between different mini-batch sizes less useful (but still interesting).
#    * Given our observations above, parameters such as the scoring overhead or the choice of activation function has to be considered carefully, especially since the model has over 12 million weights and scoring (forward propagation) is going to take a significant amount of time even on the validation frame with "only" 10,000 rows.
#    * A large part of the benchmark study was concerned with the distributed scalability (strong scaling) of H2O (at least in terms of throughput). We will address this topic in part 3, but as noted earlier, almost linear speedups can be achieved when adding compute nodes as the communication overhead is fully user-controllable.
#    * For our single-node comparison (`NODES=1`), the parameter `replicate_training_data` has no effect, so we can leave it at its default value (`TRUE`).
  
#We first run with the same parameters as the benchmark study, and then modify the parameters to improve the training speed successively.    

#score_test_set=F ## uncomment if test set error is not needed
EPOCHS=1
args <- list(
    list(hidden=c(2500,2000,1500,1000,500), epochs=EPOCHS, activation="Tanh", 
         train_samples_per_iteration=1500, validation_frame=test_hex),
    list(hidden=c(2500,2000,1500,1000,500), epochs=EPOCHS, activation="Tanh", 
         train_samples_per_iteration=1500, 
         score_training_samples=100, score_duty_cycle=0), 
    
    list(hidden=c(2500,2000,1500,1000,500), epochs=EPOCHS, activation="Tanh", 
         train_samples_per_iteration=1500, adaptive_rate=F, 
         score_training_samples=100, score_duty_cycle=0),
    
    list(hidden=c(2500,2000,1500,1000,500), epochs=EPOCHS, activation="Rectifier",
         train_samples_per_iteration=1500, max_w2=10,
         score_training_samples=100, score_duty_cycle=0),
    
    list(hidden=c(2500,2000,1500,1000,500), epochs=EPOCHS, activation="RectifierWithDropout",
         train_samples_per_iteration=1500, adaptive_rate=F, max_w2=10,
         score_training_samples=100, score_duty_cycle=0)
)
writecsv(lapply(args, run), "large_deep_net.csv")

large_deep_net <- read.csv("/data/ukbiobank/large_deep_net.csv")
plot(large_deep_net$training.speed, large_deep_net$test.set.error)
#Observations:
#    * **We cannot reproduce the slow performance (~10 samples per second) cited by the other benchmark study.** Using the same parameters, we measure a training speed of **91 samples per second** on a single i7 5820k consumer PC.
#  * When removing the scoring overhead by not specifying `validation_frame` and reducing the `score_training_samples` to 100 (while keeping the `Tanh` activation and `adaptive_rate=T`), we obtain training speeds of 100 samples per second (notably with H2O's mini-batch size of 1).
#  * When further switching to `adaptive_rate=F` (and no momentum), we obtain training speeds of 156 samples per second (albeit at the cost of accuracy).
#  * When further switching to the `Rectifier` activation function, speeds increase to **294 samples per second**. This is over *3 times faster than the original parameters and results in half the test set error!* Note that we used the `max_w2` parameter, which helps to improve the numerical stability for deep `Rectifier` networks.
#  * We observe training speeds of **520 samples per second while obtaining a lower test set error** than the original parameters, when using the `RectifierWithDropout` activation function (a good choice for MNIST, see world-record performance in Part 3 below). Hence, **H2O Deep Learning is faster on a single PC than the other method on 16 Xeon nodes** (400 samples per second with `Tanh` and AdaGrad), so further benchmark comparisons really have to take test set accuracy into account to be meaningful. Also, H2O is using a mini-batch size of 1, which means that there are many more random-access memory writes to the weights and biases.
#  * Configuring and training H2O Deep Learning models is done with a *1-liner* (function call), while other projects typically require [large configuration files](http://www.comp.nus.edu.sg/~dbsystem/singa/examples/2015/01/20/deep-simple-mlp/).
#  * H2O Deep Learning models (single- and multi-node) can be trained from R, Python, Java, Scala, JavaScript and via the Flow Web UI, as well as via JSON through the REST API. 

## Part 2 - What Really Matters: Test Set Error vs Training Time

#Now that we know how to effectively run H2O Deep Learning, we are ready to train a few models in an attempt to get good generalization performance (low test set error), which is what ultimately matters for useful Deep Learning models. Our goal is to find good models that can be trained in less than one minute, so we limit the parameter space to models that we expect to have high throughput. This is the point where some hyper-parameter search would be highly beneficial, but for simplicity, we refer to [our tutorial on hyper-parameter search](http://learn.h2o.ai/content/hands-on_training/deep_learning.html). We also add some new features such as input dropout and L1 regularization, but those are most likely to actually hurt test set performance for such short training periods, as the models barely start fitting the training data, and won't be in danger of overfitting yet. This sounds like a good topic for a future blog post...
EPOCHS=2
args <- list(
    list(epochs=EPOCHS),
    
    list(epochs=EPOCHS, activation="Tanh"),
    list(epochs=EPOCHS, hidden=c(512,512)),
    list(epochs=5*EPOCHS, hidden=c(64,128,128)),
    list(epochs=5*EPOCHS, hidden=c(512,512), 
         activation="RectifierWithDropout", input_dropout_ratio=0.2, l1=1e-5),
    
    list(epochs=5*EPOCHS, hidden=c(256,256,256), 
         activation="RectifierWithDropout", input_dropout_ratio=0.2, l1=1e-5),
    
    list(epochs=5*EPOCHS, hidden=c(200,200,200,200), 
         activation="RectifierWithDropout"),
    
    list(epochs=5*EPOCHS, hidden=c(200,200), 
         activation="RectifierWithDropout", input_dropout_ratio=0.2, l1=1e-5),
    
    list(epochs=5*EPOCHS, hidden=c(100,100,100), 
         activation="RectifierWithDropout", input_dropout_ratio=0.2, l1=1e-5),
    
    list(epochs=5*EPOCHS, hidden=c(100,100,100,100), 
         activation="RectifierWithDropout", input_dropout_ratio=0.2, l1=1e-5)
)
writecsv(lapply(args, run), "what_really_matters.csv")

what_really_matters <- read.csv("/data/ukbiobank/what_really_matters.csv")
plot(what_really_matters$test.set.error, what_really_matters$training.time)

### Summary from all models built so far
args <- list(
    list(hidden=c(256, 256), epochs=EPOCHS*10, activation="Rectifier",
         train_samples_per_iteration=1500, adaptive_rate=F, max_w2=10,
         score_training_samples=10, score_duty_cycle=0, rate=1e-3, momentum_start=0.5, momentum_ramp=1e5, momentum_stable=0.99),
    
    list(hidden=c(256, 256), epochs=EPOCHS*10, activation="RectifierWithDropout",
         train_samples_per_iteration=1500, adaptive_rate=F, max_w2=10,
         score_training_samples=100, score_duty_cycle=0, rate=1e-3, momentum_start=0.5, momentum_ramp=1e5, momentum_stable=0.99),
    
    list(hidden=c(256, 256), epochs=EPOCHS*10, activation="Tanh",
         train_samples_per_iteration=1500, adaptive_rate=F, max_w2=10,
         score_training_samples=100, score_duty_cycle=0, rate=1e-3, momentum_start=0.5, momentum_ramp=1e5, momentum_stable=0.99),
    
    list(hidden=c(512, 512), epochs=EPOCHS*10, activation="Rectifier",
         train_samples_per_iteration=1500, adaptive_rate=F, max_w2=10,
         score_training_samples=100, score_duty_cycle=0, rate=1e-3, momentum_start=0.5, momentum_ramp=1e5, momentum_stable=0.99),
    
    list(hidden=c(512, 512), epochs=EPOCHS*10, activation="RectifierWithDropout",
         train_samples_per_iteration=1500, adaptive_rate=F, max_w2=10,
         score_training_samples=100, score_duty_cycle=0, rate=1e-3, momentum_start=0.5, momentum_ramp=1e5, momentum_stable=0.99),
    
    list(hidden=c(512, 512), epochs=EPOCHS*10, activation="Tanh",
         train_samples_per_iteration=1500, adaptive_rate=F, max_w2=10,
         score_training_samples=100, score_duty_cycle=0, rate=1e-3, momentum_start=0.5, momentum_ramp=1e5, momentum_stable=0.99)
    
)

writecsv(lapply(args, run), "final_test.csv")
final_test <- read.csv("/data/ukbiobank/final_test.csv")
plot(final_test$test.set.error, final_test$training.time)