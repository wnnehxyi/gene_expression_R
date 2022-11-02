data <- read.delim("MLL_train.txt")
name1 <- data[,1]
del_index <- grep("AFFX", name1)  #抓取column中含AFFX字眼
data1_nna <- data[-del_index,]        #去除AFFX

nor <- function(x){            #Normolization
  min.x <- min(x)
  max.x <- max(x)
  a <- (x-min.x)/(max.x-min.x)
  return(a)
}

data_nor <- t(apply(data1_nna[,-c(1,2)],MARGIN = 1,nor))

group1 <- grep("ALL",colnames(data_nor))  #分類
group2 <- grep("MLL",colnames(data_nor))
group3 <- grep("AML",colnames(data_nor))

test_func <- function(data, mode){
  if(mode ==1){
    a <- t.test(data[group1], data[c(group2,group3)])
    return(a$statistic)
  }
  if(mode ==2){
    a <- t.test(data[group2], data[c(group1,group3)])
    return(a$statistic)
  }
  if(mode ==3){
    a <- t.test(data[group3], data[c(group1,group2)])
    return(a$statistic)
  }
}  #計算t-test

g1 <- apply(data_nor,MARGIN = 1, test_func, mode=1)
g2 <- apply(data_nor,MARGIN = 1, test_func, mode=2)
g3 <- apply(data_nor,MARGIN = 1, test_func, mode=3)

#合併新的資料表，並找出三類各別t-test計算出的最大值
select_data <- data.frame(cbind(data1_nna[,1], g1,g2,g3))
g1_select <- select_data[order(select_data$g1,decreasing = T),1][1:10]
g2_select <- select_data[order(select_data$g2,decreasing = T),1][1:10]
g3_select <- select_data[order(select_data$g3,decreasing = T),1][1:10]

BiocManager::install("e1071")
library(e1071)

####A 3-class classification problem using all genes
data_used <- data.frame(t(data1_nna[,-c(1,2)]))   #使用去除NA和AFFX的data
label <- factor(c(rep("ALL",20), rep("MLL",17),rep("AML",20))) #標記
data_used <- cbind(data_used, label)  #資料合併
model1 <- svm(label~., data=data_used)  #進行binary classifier
prediction <- predict(model1, data_used)  #預測模型
table(real = data_used$label, predict = prediction) #建立confusion matrix

####Three 2-class classification problems, each case using all genes
###ALL vs others
data_used <- data.frame(t(data1_nna[,-c(1,2)]))
label <- factor(c(rep("ALL",20), rep("other",37)))
data_used <- cbind(data_used, label)
model2 <- svm(label~., data=data_used)
prediction <- predict(model2, data_used)
table(real = data_used$label, predict = prediction)

###MLL vs others
data_used <- data.frame(t(data1_nna[,-c(1,2)]))
label <- factor(c(rep("other",20), rep("MLL",17),rep("other",20)))
data_used <- cbind(data_used, label)
model3 <- svm(label~., data=data_used)
prediction <- predict(model3, data_used)
table(real = data_used$label, predict = prediction)

###AML vs others
data_used <- data.frame(t(data1_nna[,-c(1,2)]))
label <- factor(c(rep("other",20), rep("other",17),rep("AML",20)))
data_used <- cbind(data_used, label)
model4 <- svm(label~., data=data_used)
prediction <- predict(model4, data_used)
table(real = data_used$label, predict = prediction)

####Three 2 class classification problems, each case using 10 selected genes
###ALL vs others 
select_gene <- order(select_data$g1,decreasing = T)[1:10] #降冪排序，選前10個ALL
data_used <- data.frame(t(data1_nna[select_gene,-c(1,2)]))
label <- factor(c(rep("ALL",20), rep("other",37)))
data_used <- cbind(data_used, label)
model5 <- svm(label~., data=data_used)
prediction <- predict(model5, data_used)
table(real = data_used$label, predict = prediction)

###MLL vs others
select_gene <- order(select_data$g2,decreasing = T)[1:10]
data_used <- data.frame(t(data1_nna[select_gene,-c(1,2)]))
label <- factor(c(rep("other",20), rep("MLL",17),rep("other",20)))
data_used <- cbind(data_used, label)
model6 <- svm(label~., data=data_used)
prediction <- predict(model6, data_used)
table(real = data_used$label, predict = prediction)

###AML vs others
select_gene <- order(select_data$g3,decreasing = T)[1:10]
data_used <- data.frame(t(data1_nna[select_gene,-c(1,2)]))
label <- factor(c(rep("other",20), rep("other",17),rep("AML",20)))
data_used <- cbind(data_used, label)
model7 <- svm(label~., data=data_used)
prediction <- predict(model7, data_used)
table(real = data_used$label, predict = prediction)

####A 3-class classification problem using the combination of 30 selected genes
select_gene <- order(select_data$g1,g2,g3,decreasing = T)[1:30]
data_used <- data.frame(t(data1_nna[select_gene,-c(1,2)]))
label <- factor(c(rep("ALL",20), rep("MLL",17),rep("AML",20)))
data_used <- cbind(data_used, label)
model8 <- svm(label~., data=data_used)
prediction <- predict(model8, data_used)
table(real = data_used$label, predict = prediction)

