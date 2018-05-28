
library(keras)

# load data
d <- read_feather("data/gage/all_gage_data.feather")

Y <- select(d,f0.02:f99.98) %>%
  mutate_all(funs(.+2)) %>%
  mutate_all(funs(log10)) 

X <- select(d,ppt_mean:tot_rdx) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

# create train and test set
set.seed(3)
ind <- sample(2, nrow(Y), replace=T, prob=c(0.9, 0.1))

# Ytrain
Ytrain <- Y[ind==1,] %>% 
  as.matrix() %>%
  split(.,col(.)) %>%
  unname() %>%
  keras_array()

# Xtrain
Xtrain <- X[ind==1,]

# Ytest
Ytest <- Y[ind==2,] %>% 
  as.matrix() %>%
  split(.,col(.)) %>%
  unname() %>%
  keras_array()

# X test
Xtest <- X[ind==2,]

# build model ----
# caclulate weights
my_fun <- function(x){mean(Y$f0.02/x)}
weights <- as.numeric(round(apply(Y,2,my_fun),3))^2

# add covariates
input <- layer_input(shape=dim(X)[2],name="basinchars")

base_model <- input  %>%
  layer_dense(units = 34,activation="relu") %>%
  layer_dropout(rate=0.1) 

for(i in 1:dim(Y)[2]){
  y <- colnames(Y)[i]
  outstring <- paste0(
    sprintf("%s <- base_model %%>%%", y), 
    sprintf(" layer_dense(units = 1, activation='relu', name='%s')",y)
  )
  eval(parse(text=outstring))
}

Ylist <- paste0("list(",paste(colnames(Y),sep="",collapse=","),")")
model <- keras_model(input,eval(parse(text=Ylist))) %>%
  compile(optimizer = "rmsprop",
          loss="mse",
          metrics="mae",
          loss_weights=weights)

model_fit <- model %>% 
  fit(x=Xtrain,
      y=Ytrain, 
      epochs=100, 
      batch_size=50,
      #validation_data = list(Xtest, Ytest), 
      validation_split = 0.3,
      verbose=0)

# check model fit ----
data.frame(loss=sqrt(model_fit$metrics$val_loss),
           epoch=1:100) %>%
  ggplot() + 
  geom_line(aes(epoch,loss)) +
  theme_bw()

# predict values for test set
Yhat <- predict(model, Xtest) %>%
  data.frame() %>%
  setNames(colnames(Y))

# plot est vs obs
Ytest <- Y[ind==2,] %>%
  rownames_to_column() %>%
  gather(variable,obs,-rowname)

est_obs <- Yhat %>%
  rownames_to_column() %>%
  gather(variable,yhat,-rowname) %>%
  left_join(Ytest,by=c("rowname","variable")) %>%
  mutate(obs = (10^obs)-2,
         yhat = round((10^yhat)-2,2),
         variable = as.numeric(substring(variable, 2))) %>%
  gather(model,value,-rowname,-variable) %>%
  filter(rowname %in% as.character(sample(1:306,12)))

ggplot(est_obs,aes(variable,value,linetype=model)) +
  geom_line() +
  #geom_point() +
  facet_wrap(~rowname,scales="free_y") +
  scale_y_log10() +
  labs(x="Exceedence probability",y="log10(Q)") +
  theme_bw()

hold <- Yhat %>%
  rownames_to_column() %>%
  gather(variable,yhat,-rowname) %>%
  mutate(variable = as.numeric(substring(variable, 2)))

ggplot(hold) + geom_line(aes(variable,yhat,group=rowname)) 

# counts monotonic violations
not_increasing <- function(x){sum(cummax(x)!=x)}

# sum violations
sum(apply(Yhat,1,not_increasing))
sum(apply(Ytest,1,not_increasing))

