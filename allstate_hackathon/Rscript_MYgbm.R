library(data.table) # tutorial assumes  v1.9.7
library(gbm) # CRAN v2.1.1
setwd("/Users/dgrayson1/Desktop/allstate_hackathon")
rm(list=ls())

data_all <- fread('train.csv')
ids_test <- fread('test.csv')

data_all[, dummy_counter:=1]

# find the first/last event for each policy in training data
data_all[, max_timestamp_by_policy:=max(timestamp), by=id]
data_all[, min_timestamp_by_policy:=min(timestamp), by=id]

# featurize data, without the last element, by looking at counts of events
data_all_allButLast <- subset(data_all, timestamp != max_timestamp_by_policy)
data_all_ohe <- dcast(data_all_allButLast, 
                      id ~ event, 
                      fun.aggregate=sum, 
                      value.var="dummy_counter")

# create a feature that counts how many events occured in total before the last event
data_all_ohe[, total_events:=apply(data_all_ohe[, -c("id"), with=FALSE], 1, sum)]

# create feature: first event for each policy
data_all_first_element = subset(data_all, timestamp==min_timestamp_by_policy)
data_all_first_element$timestamp = NULL
data_all_first_element$dummy_counter = NULL
data_all_first_element$max_timestamp_by_policy = NULL
data_all_first_element$min_timestamp_by_policy = NULL
setnames(data_all_first_element, "event", "first_event")
# add first event into training data
data_all_ohe <- merge(data_all_ohe, data_all_first_element, by=c("id"))
data_all_ohe$first_event <- as.factor(data_all_ohe$first_event)

# create feature: 2nd_to_last event (tminus1event) for each policy
data_all_tminus1event = subset(data_all, timestamp == (max_timestamp_by_policy-1))
data_all_tminus1event$timestamp = NULL
data_all_tminus1event$dummy_counter = NULL
data_all_tminus1event$max_timestamp_by_policy = NULL
data_all_tminus1event$min_timestamp_by_policy = NULL
setnames(data_all_tminus1event, "event", "tminus1event")
# add tminus1 event into training data
data_all_ohe <- merge(data_all_ohe, data_all_tminus1event, by=c("id"))
data_all_ohe$tminus1event <- as.factor(data_all_ohe$tminus1event)

# create feature: 3rd_to_last event (tminus2event) for each policy
#data_all_tminus2event = subset(data_all, timestamp == (max_timestamp_by_policy-2))
#data_all_tminus2event$timestamp = NULL
#data_all_tminus2event$dummy_counter = NULL
#data_all_tminus2event$max_timestamp_by_policy = NULL
#data_all_tminus2event$min_timestamp_by_policy = NULL
#setnames(data_all_tminus2event, "event", "tminus2event")
# add tminus2 event into training data
#data_all_ohe <- merge(data_all_ohe, data_all_tminus2event, by=c("id"))
#data_all_ohe$tminus2event <- as.factor(data_all_ohe$tminus2event)

# create target variable: last event for each policy
data_all_last_element = subset(data_all, timestamp == max_timestamp_by_policy)
data_all_last_element$timestamp = NULL
data_all_last_element$dummy_counter = NULL
data_all_last_element$max_timestamp_by_policy = NULL
data_all_last_element$min_timestamp_by_policy = NULL
setnames(data_all_last_element, "event", "last_event")

# add dependent variable (last event) into dataframe
data_all_merged <- merge(data_all_ohe, data_all_last_element, by=c("id"))

#NOW, get subset of training data not part of test data subjects
data_train_merged <- subset(data_all_merged, !(id %in% ids_test$id))

# train a gbm model. For faster run times, lower n.trees
#you probably want to cut off n.trees where Improve starts going negative
model <- gbm(last_event ~ . - id, 
             distribution="multinomial", 
             data=data_train_merged, 
             interaction.depth=4,
             n.trees=100,
             shrinkage=0.1,
             train.fraction=0.5,
             verbose=TRUE) #,cv.folds = 5

#look at in-sample model error as fcn of iterations 
gbm.perf(model,oobag.curve = TRUE, overlay = TRUE)

#look at relative contributions of each predictor
summary(model, las=2)

##############################
##############################
# create the test dataset
data_test_ohe <- subset(data_all_ohe, id %in% ids_test$id)

# score the model
predictions_raw <- predict(model, data_test_ohe, type="response")

#package it into a format that can be scored
predictions <- as.data.table(predictions_raw[,,1])
setnames(predictions, names(predictions), colnames(predictions_raw))
predictions[, id:=data_test_ohe$id]

pred_columns <- colnames(predictions)
pred_columns <- c('id', sort(pred_columns[-length(pred_columns)]))
setcolorder(predictions, pred_columns)
setnames(predictions, c(pred_columns[1], 
                        paste('event_', pred_columns[-1], sep='')))

write.csv(predictions, file="Rpredictions.csv", quote=TRUE, row.names=FALSE)