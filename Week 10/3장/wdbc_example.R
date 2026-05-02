
#install.packages("FastKNN")
library(FastKNN)
library(gmodels)

setwd("G:/내 드라이브/4 Teaching_GSPH/2026_MLstudy/Week 10_Excercise with clinical data/ai")
wdbc <- read.csv("G:/내 드라이브/4 Teaching_GSPH/2026_MLstudy/Week 10_Excercise with clinical data/ai/3장/wdbc.csv")
table(wdbc$Class) 

wdbc <- wdbc[sample(1:nrow(wdbc)),]
wdbc_p <- wdbc[ ,3:32]
normalize <-function(x) { return ((x-min(x))/(max(x)-min(x))) }
wdbc_p <- as.data.frame(lapply(wdbc_p[1:30], normalize))

wdbc_train <- wdbc_p[1:467, ]
wdbc_test <- wdbc_p[468:567, ]

wdbc_train_labels <- wdbc[1:467, 2]
wdbc_test_labels <- wdbc[468:567, 2]

wdbc_distance <- Distance_for_KNN_test(wdbc_test, wdbc_train)
wdbc_prediction <- knn_test_function(wdbc_train, wdbc_test, wdbc_distance, wdbc_train_labels, k=11)

result_table <- data.frame(wdbc[468:567, 1], wdbc_test_labels, wdbc_prediction)
write.csv(result_table, file= "wdbc_result.csv")

CrossTable(wdbc_test_labels, wdbc_prediction)
