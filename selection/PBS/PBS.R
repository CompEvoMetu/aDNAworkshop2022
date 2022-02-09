#!/usr/bin/env Rscript

af_euas=read.delim("af_euas_norep.windowed.weir.fst")
af_euas=af_euas[,c(2,4,5)]
colnames(af_euas)=c("POS","n_var","af_euas")

af_eas=read.delim("af_eas_norep_windowed.weir.fst")
af_eas=af_eas[,c(2,4,5)]
colnames(af_eas)=c("POS","n_var","af_eas")

eas_euas=read.delim("eas_euas_norep_windowed.weir.fst")
eas_euas=eas_euas[,c(2,4,5)]
colnames(eas_euas)=c("POS","n_var","eas_euas")

merged=merge(af_euas,af_eas,by="POS")
merged=merge(merged,eas_euas,by="POS")

Fst=merged[,c(3,5,7)]
Fst[Fst<0]=0

t_Fst=-log10(1-Fst)

merged$PBS=(t_Fst$af_eas+t_Fst$eas_euas-t_Fst$af_euas)/2

jpeg("plot.jpeg",width = 1000,height = 500)
plot(merged$PBS~merged$POS,pch=20,col="gray45",xlab="POS",ylab="PBS",las=1)
abline(h=mean(merged$PBS)+3*sqrt(var(merged$PBS)))
dev.off()

