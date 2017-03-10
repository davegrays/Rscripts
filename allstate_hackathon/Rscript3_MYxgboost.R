library(data.table) # tutorial assumes  v1.9.7
library(gbm) # CRAN v2.1.1
library(xgboost)
setwd("/Users/dgrayson1/Desktop/allstate_hackathon")
rm(list=ls())

data_all <- fread('train.csv')
ids_test <- fread('test.csv')

data_all[, dummy_counter:=1]

data_train <- subset(data_all, !(id %in% ids_test$id))

# find the first/last event for each policy in training data
data_train[, max_timestamp_by_policy:=max(timestamp), by=id]
data_train[, min_timestamp_by_policy:=min(timestamp), by=id]

data_train_last_element = subset(data_train, timestamp == max_timestamp_by_policy)
data_train_last_element$timestamp = NULL
data_train_last_element$dummy_counter = NULL
data_train_last_element$max_timestamp_by_policy = NULL
data_train_last_element$min_timestamp_by_policy = NULL
setnames(data_train_last_element, "event", "last_event")

data_train_first_element = subset(data_train, timestamp==min_timestamp_by_policy)
data_train_first_element$timestamp = NULL
data_train_first_element$dummy_counter = NULL
data_train_first_element$max_timestamp_by_policy = NULL
data_train_first_element$min_timestamp_by_policy = NULL
setnames(data_train_first_element, "event", "first_event")

# featurize data, without the last element, by looking at counts of events
data_train_allButLast <- subset(data_train, timestamp != max_timestamp_by_policy)
data_train_ohe <- dcast(data_train_allButLast, 
                        id ~ event, 
                        fun.aggregate=sum, 
                        value.var="dummy_counter")

# create a feature that counts how many events occured in total before the last event
data_train_ohe[, total_events:=apply(data_train_ohe[, -c("id"), with=FALSE], 1, sum)]

data_train_ohe <- merge(data_train_ohe, data_train_first_element, by=c("id"))
data_train_ohe$first_event <- as.factor(data_train_ohe$first_event)

# create feature: 2nd_to_last event (tminus1event) for each policy
data_train_tminus1event = subset(data_train, timestamp == (max_timestamp_by_policy-1))
data_train_tminus1event$timestamp = NULL
data_train_tminus1event$dummy_counter = NULL
data_train_tminus1event$max_timestamp_by_policy = NULL
data_train_tminus1event$min_timestamp_by_policy = NULL
setnames(data_train_tminus1event, "event", "tminus1event")
# add tminus1 event into training data
data_train_ohe <- merge(data_train_ohe, data_train_tminus1event, by=c("id"))
data_train_ohe$tminus1event <- as.factor(data_train_ohe$tminus1event)

# create feature: 3rd_to_last event (tminus2event) for each policy
data_train_tminus2event = subset(data_train, timestamp == (max_timestamp_by_policy-2))
data_train_tminus2event$timestamp = NULL
data_train_tminus2event$dummy_counter = NULL
data_train_tminus2event$max_timestamp_by_policy = NULL
data_train_tminus2event$min_timestamp_by_policy = NULL
setnames(data_train_tminus2event, "event", "tminus2event")
# add tminus2 event into training data
data_train_ohe <- merge(data_train_ohe, data_train_tminus2event, by=c("id"))
data_train_ohe$tminus2event <- as.factor(data_train_ohe$tminus2event)

# create feature: 4th_to_last event (tminus3event) for each policy
data_train_tminus3event = subset(data_train, timestamp == (max_timestamp_by_policy-3))
data_train_tminus3event$timestamp = NULL
data_train_tminus3event$dummy_counter = NULL
data_train_tminus3event$max_timestamp_by_policy = NULL
data_train_tminus3event$min_timestamp_by_policy = NULL
setnames(data_train_tminus3event, "event", "tminus3event")
# add tminus3 event into training data
data_train_ohe <- merge(data_train_ohe, data_train_tminus3event, by=c("id"), all = TRUE)
data_train_ohe$tminus3event <- as.factor(data_train_ohe$tminus3event)

#these need to be one-hot encoded!
#i.e.
#
# totalnumeventsGRTRthan3, num30018events_iftotalGRTRthan4, num30021"", num30024"", ..., firsteventiftotalGRTRthan2 is 30018, "" is 30021, "" is 30024, ..., tminus1event is 30018, "" is 30021, "" 30024, ..., tminus2event is 30018, "" is 30021, "" 30024, ...,
# ...
#
# translates to
#
# 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
# ...

data_train_ohe$id <- NULL
data_train_last_element$id <- NULL

data_train_last_element[data_train_last_element == 30018] <- 0
data_train_last_element[data_train_last_element == 30021] <- 1
data_train_last_element[data_train_last_element == 30024] <- 2
data_train_last_element[data_train_last_element == 30027] <- 3
data_train_last_element[data_train_last_element == 30039] <- 4
data_train_last_element[data_train_last_element == 30042] <- 5
data_train_last_element[data_train_last_element == 30045] <- 6
data_train_last_element[data_train_last_element == 30048] <- 7
data_train_last_element[data_train_last_element == 36003] <- 8
data_train_last_element[data_train_last_element == 45003] <- 9

model <- xgboost(data = data_train_ohe,
                label = data_train_last_element,
                max.depth = 6,
                objective = "multi:softmax",
                num_class = 10,
                verbose = 2) #eta = 1, nthread = 2, nround = 2,
                 



##try n.trees=2000 and shrinkage=0.005
##try n.trees=200 and shrinkage=0.05
##try n.trees=100 and shrinkage=0.1

## if you are not using train.fraction, you need to update the n.trees further below in the model prediciton
## if you are using it, then you need to remove n.trees further below

#also must edit the output filename to include number of trees and such




#look at model error as fcn of iterations 
gbm.perf(model,oobag.curve = TRUE, overlay = TRUE)

#look at relative contributions of each predictor
#summary(model, las=2)

# create the test dataset
data_test <- subset(data_all, id %in% ids_test$id)
data_test[, min_timestamp_by_policy:=min(timestamp), by=id]
data_test[, max_timestamp_by_policy:=max(timestamp), by=id]

# create the first event feature
data_test_first_event <- subset(data_test, timestamp==min_timestamp_by_policy)
data_test_first_event$timestamp = NULL
data_test_first_event$dummy_counter = NULL
data_test_first_event$min_timestamp_by_policy = NULL
setnames(data_test_first_event, "event", "first_event")

data_test_ohe <- dcast(data_test, id~event, fun.aggregate=sum, value.var="dummy_counter")
# add the total events feature
data_test_ohe[, total_events:=apply(data_test_ohe[, -c("id"), with=FALSE], 1, sum)]

# add the first event feature
data_test_ohe <- merge(data_test_ohe, data_test_first_event, by=c("id"))

# create feature: 2nd_to_last event (tminus1event) for each policy
data_test_tminus1event = subset(data_test, timestamp == max_timestamp_by_policy)
data_test_tminus1event$timestamp = NULL
data_test_tminus1event$dummy_counter = NULL
data_test_tminus1event$max_timestamp_by_policy = NULL
data_test_tminus1event$min_timestamp_by_policy = NULL
setnames(data_test_tminus1event, "event", "tminus1event")
# add tminus1 event into training data
data_test_ohe <- merge(data_test_ohe, data_test_tminus1event, by=c("id"))
data_test_ohe$tminus1event <- as.factor(data_test_ohe$tminus1event)

# create feature: 3rd_to_last event (tminus2event) for each policy
data_test_tminus2event = subset(data_test, timestamp == (max_timestamp_by_policy-1))
data_test_tminus2event$timestamp = NULL
data_test_tminus2event$dummy_counter = NULL
data_test_tminus2event$max_timestamp_by_policy = NULL
data_test_tminus2event$min_timestamp_by_policy = NULL
setnames(data_test_tminus2event, "event", "tminus2event")
# add tminus2 event into training data
data_test_ohe <- merge(data_test_ohe, data_test_tminus2event, by=c("id"))
data_test_ohe$tminus2event <- as.factor(data_test_ohe$tminus2event)

# create feature: 4th_to_last event (tminus3event) for each policy
data_test_tminus3event = subset(data_test, timestamp == (max_timestamp_by_policy-3))
data_test_tminus3event$timestamp = NULL
data_test_tminus3event$dummy_counter = NULL
data_test_tminus3event$max_timestamp_by_policy = NULL
data_test_tminus3event$min_timestamp_by_policy = NULL
setnames(data_test_tminus3event, "event", "tminus3event")
# add tminus3 event into testing data
data_test_ohe <- merge(data_test_ohe, data_test_tminus3event, by=c("id"), all = TRUE)
data_test_ohe$tminus3event <- as.factor(data_test_ohe$tminus3event)

# score the model
predictions_raw <- predict(model, data_test_ohe, type="response", n.trees=800)
#package it into a format that can be scored
predictions <- as.data.table(predictions_raw[,,1])
setnames(predictions, names(predictions), colnames(predictions_raw))
predictions[, id:=data_test_ohe$id]

pred_columns <- colnames(predictions)
pred_columns <- c('id', sort(pred_columns[-length(pred_columns)]))
setcolorder(predictions, pred_columns)
setnames(predictions, c(pred_columns[1], 
                        paste('event_', pred_columns[-1], sep='')))

write.csv(predictions, file="results_Rscript3_800trees_ALLTRAIN.csv", quote=TRUE, row.names=FALSE)
