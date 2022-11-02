data <- read.delim("SRBCT_train.txt") #讀檔
data1 <- as.matrix(data[,-c(1:2)])  #將檔案第一、二個column去除
name1 <- data[,1]     #name1為原始檔中第一個column的內容
nor <- function(x){   #建立一個function進行Normalization
  min.x <- min(x)
  max.x <- max(x)
  (x-min.x)/(max.x-min.x)
}

data1_nor <- t(apply(data1,1,nor))  #進行transpose

#選出所要的Image.Id，並將各群指定顏色
A <- data1_nor[which(name1 == "770394"),] 
B <- data1_nor[which(name1 == "814260"),]
C <- data1_nor[which(name1 == "491565"),]
color <- rep(0,63)
color[grep("EWS",names(A))] <- 1
color[grep("BL",names(A))] <- 2
color[grep("NB",names(A))] <- 3
color[grep("RMS",names(A))] <-4

#畫三張圖並加上圖例
plot(A,xlab = "Samples",ylab = "genes expression values",pch=color, col=color)   
legend("topright", legend=c("EWS","BL","NB","RMS"), pch = c(1:4), col = c(1:4),bty="n")
plot(B,xlab = "Samples",ylab = "genes expression values",pch=color, col=color)
legend("topright", legend=c("EWS","BL","NB","RMS"), pch = c(1:4), col = c(1:4),bty="n")
plot(C,xlab = "Samples",ylab = "genes expression values",pch=color, col=color)
legend("topright", legend=c("EWS","BL","NB","RMS"), pch = c(1:4), col = c(1:4),bty="n")

#找出另外三個Image.Id，並額外畫兩張versus圖
D <- data1_nor[which(name1 == "236282"),]
E <- data1_nor[which(name1 == "812105"),]
G <- data1_nor[which(name1 == "784224"),]
plot(A,D, pch=color, col = color, xlab = "770394", ylab = "236282")
legend("topright", legend=c("EWS","BL","NB","RMS"), pch = c(1:4), col = c(1:4),bty="n")
plot(E,G, pch=color, col = color, xlab = "812105", ylab = "784224")
legend("topright", legend=c("EWS","BL","NB","RMS"), pch = c(1:4), col = c(1:4),bty="n")


