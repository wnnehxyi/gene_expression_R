aa <- read.delim("SRBCT_train.txt")  #Ū��
bb <- as.matrix(aa[,-c(1,2)])  #�N���椤�Ĥ@�B�G��column�h��
mean_gene <- apply(bb,1,mean)  #�p��C�@��row��������(1��row)
sd_gene <- apply(bb,1,sd)      #�p��C�@��row���зǮt(1��row)
mean_sample <- apply(bb,2,mean)  #�p��C�@��column��������(2��column)
sd_sample <- apply(bb,2,sd)  #�p��C�@��column���зǮt(2��column)
cc <- cbind(aa,sd_gene)   #�N�зǮt���C��row���X��(�|�[�Ҧ�row)
dd <- cc[order(cc$sd_gene,decreasing = 1),]  #�X�᪺֫row�H�����覡�i��ƦC�A�åHindex���Ƨ�
top30 <- dd$Image.Id[1:round(2308*0.3,0)] #�N�e30%�j���зǮt�A�HImage.Id���覡�e�{
