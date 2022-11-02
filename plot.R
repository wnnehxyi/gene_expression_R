data <- read.delim("SRBCT_train.txt") #Ū��
data1 <- as.matrix(data[,-c(1:2)])  #�N�ɮײĤ@�B�G��column�h��
name1 <- data[,1]     #name1����l�ɤ��Ĥ@��column�����e
nor <- function(x){   #�إߤ@��function�i��Normalization
  min.x <- min(x)
  max.x <- max(x)
  (x-min.x)/(max.x-min.x)
}

data1_nor <- t(apply(data1,1,nor))  #�i��transpose

#��X�ҭn��Image.Id�A�ñN�U�s���w�C��
A <- data1_nor[which(name1 == "770394"),] 
B <- data1_nor[which(name1 == "814260"),]
C <- data1_nor[which(name1 == "491565"),]
color <- rep(0,63)
color[grep("EWS",names(A))] <- 1
color[grep("BL",names(A))] <- 2
color[grep("NB",names(A))] <- 3
color[grep("RMS",names(A))] <-4

#�e�T�i�Ϩå[�W�Ϩ�
plot(A,xlab = "Samples",ylab = "genes expression values",pch=color, col=color)   
legend("topright", legend=c("EWS","BL","NB","RMS"), pch = c(1:4), col = c(1:4),bty="n")
plot(B,xlab = "Samples",ylab = "genes expression values",pch=color, col=color)
legend("topright", legend=c("EWS","BL","NB","RMS"), pch = c(1:4), col = c(1:4),bty="n")
plot(C,xlab = "Samples",ylab = "genes expression values",pch=color, col=color)
legend("topright", legend=c("EWS","BL","NB","RMS"), pch = c(1:4), col = c(1:4),bty="n")

#��X�t�~�T��Image.Id�A���B�~�e��iversus��
D <- data1_nor[which(name1 == "236282"),]
E <- data1_nor[which(name1 == "812105"),]
G <- data1_nor[which(name1 == "784224"),]
plot(A,D, pch=color, col = color, xlab = "770394", ylab = "236282")
legend("topright", legend=c("EWS","BL","NB","RMS"), pch = c(1:4), col = c(1:4),bty="n")
plot(E,G, pch=color, col = color, xlab = "812105", ylab = "784224")
legend("topright", legend=c("EWS","BL","NB","RMS"), pch = c(1:4), col = c(1:4),bty="n")

