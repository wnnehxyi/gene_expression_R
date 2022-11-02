aa <- read.delim("SRBCT_train.txt")  #讀檔
bb <- as.matrix(aa[,-c(1,2)])  #將表格中第一、二個column去除
mean_gene <- apply(bb,1,mean)  #計算每一個row的平均值(1為row)
sd_gene <- apply(bb,1,sd)      #計算每一個row的標準差(1為row)
mean_sample <- apply(bb,2,mean)  #計算每一個column的平均值(2為column)
sd_sample <- apply(bb,2,sd)  #計算每一個column的標準差(2為column)
cc <- cbind(aa,sd_gene)   #將標準差的每個row做合併(疊加所有row)
dd <- cc[order(cc$sd_gene,decreasing = 1),]  #合併後的row以降冪方式進行排列，並以index為排序
top30 <- dd$Image.Id[1:round(2308*0.3,0)] #將前30%大的標準差，以Image.Id的方式呈現

