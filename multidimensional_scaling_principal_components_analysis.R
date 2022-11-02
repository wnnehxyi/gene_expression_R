library(ggplot2)
data_HW3 <- read.delim("SRBCT_train.txt") 
data1_HW3 <- as.matrix(data_HW3[,-c(1:2)])  
name1_HW3 <- data_HW3[,1]    
nor <- function(x){   
  min.x <- min(x)
  max.x <- max(x)
  (x-min.x)/(max.x-min.x)
}

data1_nor_HW3 <- t(apply(data1_HW3,1,nor))  

#Scatter plot for Image ID "770394" in HW3
A <- data1_nor_HW3[which(name1_HW3 == "770394"),] 

color <- rep(0,63)
color[grep("EWS",names(A))] <- 1
color[grep("BL",names(A))] <- 2
color[grep("NB",names(A))] <- 3
color[grep("RMS",names(A))] <-4

use_data <- data.frame(A)
use_data[,2]<-c(1:63)   #一共有63個樣本(編號)
names(use_data) <-c("genes_expression_values","Samples")

ggplot(use_data,aes(x=Samples,y=genes_expression_values))+
  geom_point(shape=16,size=4,col=color)+
  geom_smooth(method = lm)


#Scatter plot for Image ID "770394" vs"236382"
use_data1 <- data.frame(A)
use_data1[,2] <- data1_nor_HW3[which(name1_HW3=="236282"),]
names(use_data1) <-c("ID_236282","ID_770394")

ggplot(use_data1,aes(x=ID_770394,y=ID_236282))+
  geom_point(shape=16,size=4,col=color)


##MDS and PCA
data <- read.delim("MLL_train.txt")
name1 <- data[,1]
del_index <- grep("AFFX", name1)
data_noAF <- data[-del_index,]

nor <- function(x){
  min.x <- min(x)
  max.x <- max(x)
  a <- (x-min.x)/(max.x-min.x)
  return(a)
}

data_nor <- t(apply(data_noAF[,-c(1,2)],MARGIN = 1,nor))

group1 <- grep("ALL",colnames(data_nor))
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
}  

g1 <- apply(data_nor,MARGIN = 1, test_func, mode=1)
g2 <- apply(data_nor,MARGIN = 1, test_func, mode=2)
g3 <- apply(data_nor,MARGIN = 1, test_func, mode=3)

select_data <- data.frame(cbind(data_noAF[,1], g1,g2,g3))
g1_select <- select_data[order(select_data$g1,decreasing = T),1][1:10]
g2_select <- select_data[order(select_data$g2,decreasing = T),1][1:10]
g3_select <- select_data[order(select_data$g3,decreasing = T),1][1:10]
select_gene <- order(select_data$g1,g2,g3,decreasing = T)[1:30]
data_used_plot <- data.frame(t(data_noAF[select_gene,-c(1,2)]))

#MDS
data_mds1<-cmdscale(dist(data_used_plot[,1:30]),2)
data_MDS <- data_mds1[,1]
data_MDS_frame <- data.frame(data_MDS)
data_MDS_frame[,2] <- data_mds1[,2]

color <- rep(1,3)
color[grep("ALL",rownames(data_MDS_frame))] <- 1
color[grep("MLL",rownames(data_MDS_frame))] <- 2
color[grep("AML",rownames(data_MDS_frame))] <- 3

ggplot(data_MDS_frame,aes(x = data_MDS,y = V2))+
  geom_point(size=3,color=color)+
  xlab("1st direction")+
  ylab("2nd direction")


#PCA
data.prcomp1<-prcomp(data_used_plot[,1:30], scale=T, retx=T)
data_PCA <- data_mds1[,1]
data_PCA_frame <- data.frame(data_MDS)
data_PCA_frame[,2] <- data_mds1[,2]

ggplot(data_PCA_frame,aes(x = data_MDS,y = V2))+
  geom_point(size=3,color=color)+     
  xlab("1st PC")+
  ylab("2nd PC")

