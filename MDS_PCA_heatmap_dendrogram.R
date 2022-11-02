data <- read.delim("MLL_train.txt")
name1 <- data[,1]
del_index <- grep("AFFX", name1)
data1 <- data[-del_index,]

nor <- function(x){
  min.x <- min(x)
  max.x <- max(x)
  a <- (x-min.x)/(max.x-min.x)
  return(a)
}

data_nor <- t(apply(data1[,-c(1,2)],MARGIN = 1,nor))

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

select_data <- data.frame(cbind(data1[,1], g1,g2,g3))
g1_select <- select_data[order(select_data$g1,decreasing = T),1][1:10]
g2_select <- select_data[order(select_data$g2,decreasing = T),1][1:10]
g3_select <- select_data[order(select_data$g3,decreasing = T),1][1:10]
select_gene <- order(select_data$g1,g2,g3,decreasing = T)[1:30]
data_used_plot <- data.frame(t(data1[select_gene,-c(1,2)]))

#MDS
data.mds1<-cmdscale(dist(data_used_plot[,1:30]),2)
plot(c(min(data.mds1[,1]), max(data.mds1[,1])), c(min(data.mds1[,2]), max(data.mds1[,2])), type="n", xlab="1st direction", ylab="2nd direction")
text(data.mds1[1:20,1], data.mds1[1:20,2], 1:20,col=1)
text(data.mds1[21:37,1], data.mds1[21:37,2], 21:37,col=2)
text(data.mds1[38:57,1], data.mds1[38:57,2], 38:57,col=3)

#PCA
data.prcomp1<-prcomp(data_used_plot[,1:30], scale=T, retx=T)
plot(c(min(data.prcomp1$x[,1]), max(data.prcomp1$x[,1])), c(min(data.prcomp1$x[,2]), max(data.prcomp1$x[,2])), type="n", xlab="1st PC", ylab="2nd PC")
text(data.prcomp1$x[1:20,1], data.prcomp1$x[1:20,2], 1:20,col=1)
text(data.prcomp1$x[21:37,1], data.prcomp1$x[21:37,2], 21:37,col=2)
text(data.prcomp1$x[38:57,1], data.prcomp1$x[38:57,2], 38:57,col=3)

#data.frame(heatmap and dendrogram)
data_used <- data.frame(t(data1[select_gene,-c(2)]))
colnames(data_used) <- data_used[1,]
A <- data_used[1,]
data_used <- data_used[-1,]
data_used1 <- t(apply(data_used,1,as.numeric))
colnames(data_used1) <- A

#heatmap and dendrogram
library(gplots)
data.kmeans<-kmeans(as.matrix(data_used1), 10)
clustered.data<-as.matrix(data_used1[data.kmeans$cluster==1,1:30]) 
for(i in 2:10)
{
  clustered.data<-rbind(clustered.data, as.matrix(data_used1[data.kmeans$cluster==i, 1:30]))
}
my_dist <- function(x){
  dist(x,method="euclidean")
}
my_hclust<- function(x){
  hclust(x,method="average")
}
heatmap.2(as.matrix(clustered.data), scale= "none",
          col=rgb(c(rep(0, 16), 0, 1:16)/16, c(16:1, 0, rep(0, 16))/16, rep(0, 33)),
          distfun=my_dist, hclustfun=my_hclust, dendrogram="row")
