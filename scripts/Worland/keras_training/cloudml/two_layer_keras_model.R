library(keras)
library(tfruns)
library(tidyverse)
library(feather)
library(cloudml)

# set working directory
# setwd("scripts/Worland/keras_training/cloudml")
# cloudml::cloudml_train("two_layer_keras_model.R", config = "tuning.yml")
# trials <- job_collect('cloudml_2018_06_22_164321294', trials = 'all')
# trials <- job_trials("cloudml_2018_06_22_141531355") # use this one

# load data file
d <- read_feather('all_gage_data.feather')

f15 <- c("f0.02","f0.5","f05","f10","f20", "f30","f40","f50",
         "f60","f70","f80","f90","f95","f99.5","f99.98")

# grab names for later
ynames <- colnames(select(d,f15))

# create multioutput array
Y <- select(d,f15) %>%
  mutate_all(funs(.+2)) %>%
  mutate_all(funs(log10)) %>%
  as.matrix() %>%
  split(.,col(.)) %>%
  unname() %>%
  keras_array()

# create predictor matrix
X <- select(d,major:flood_storage) %>%
  mutate_all(funs(as.numeric(as.factor(.)))) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.)))  %>%
  as.matrix()

# establish flags
FLAGS <- flags(
  flag_integer("dense_units1", 40, "Dense units in first layer"),
  flag_numeric("dropout1", 0.1, "Dropout in first layer"),
  flag_integer("dense_units2", 30, "Dense units in second layer"),
  flag_numeric("dropout2", 0.1, "Dropout in second layer"),
  flag_numeric("learning_rate", 0.001, "learning rate for rmsprop"),
  flag_numeric("epochs", 75, "number of epochs"),
  flag_numeric("batch_size", 25, "batch size for training")
)

# Build model function for multiple outputs
build_model <- function(){
  input <- layer_input(shape=dim(X)[2],name="basinchars")

  base_model <- input  %>%
    layer_dense(units = FLAGS$dense_units1,
                activation="relu") %>%
    layer_dropout(rate=FLAGS$dropout1) %>%
    layer_dense(units = FLAGS$dense_units2,
                activation="relu") %>%
    layer_dropout(rate=FLAGS$dropout2)

  for(i in 1:length(ynames)){
    y <- ynames[i]
    outstring <- paste0(
      sprintf("%s <- base_model %%>%%", y),
      sprintf(" layer_dense(units = 1, activation='relu', name='%s')",y)
    )
    eval(parse(text=outstring))
  }

  Ylist <- paste0("list(",paste(ynames,sep="",collapse=","),")")
  model <- keras_model(input,eval(parse(text=Ylist))) %>%
    compile(optimizer = optimizer_rmsprop(lr = FLAGS$learning_rate),
            loss="mse",
            metrics="mae")

  return(model)
}

# build model
model <- build_model()

# fit model
model_fit <- model %>%
  fit(x=X,
      y=Y,
      epochs=FLAGS$epochs,
      batch_size = FLAGS$batch_size,
      validation_split = 0.25,
      verbose=0)

