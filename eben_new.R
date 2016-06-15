setwd('~/Dropbox (Personal)/Coloncancer/');
library('EBEN');


mi <- read.table("bc_matrix.txt",header=T);
mi <- as.matrix(mi);
mi2 <- mi[,which(mi[nrow(mi),]!='NA')];
#write.table(mi2,"mi2",sep='\t',quote=F,col.names = T,row.names = F);
target <- as.matrix(mi2[nrow(mi2),(2:ncol(mi2))]);
target <- as.numeric(target);
target1 <- log(target,base=exp(1));
x <- mi2[(1:nrow(mi2)-1),];
x11 <- x[,2:ncol(x)];
x11 <- matrix(as.numeric(x11),nrow(x11));

x1 <- NULL;
for(i in 1:nrow(x11)){
  if(sum(as.numeric(x11[i,])!=0)){x1 <- rbind(x1,x[i,]);}
}

x2 <- NULL;
criteria <- trunc((ncol(x1)-1) *0.8);  
for(i in 1:nrow(x1)){
  if(sum(as.numeric(x1[i,(2:ncol(x1))])!=0) > criteria){
    x2 <- rbind(x2,x1[i,]);
  }
}
colnames(x2) <- colnames(mi2);
new_matrix_miRNA <- rbind(x2,c("stage",target));
#new_matrix_miRNA.txt is the final filtered data: 377*234 ( including header and row names);

############standarlization:
#temp <- x2[(2:nrow(x2)),];
#temp <- matrix(as.numeric(temp),nrow(temp));
#x3 <- scale(temp, center = TRUE, scale = apply(temp, 2, sd, na.rm = TRUE));


##########quantile normalization:
x3 <- x2[,2:ncol(x2)];
for( sl in 1:nrow(x3) ) {
  mat = matrix(as.numeric(x3[sl,]),1);
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(x3)+1));
  x3[sl,] = mat;
}
rm(sl, mat);
##new <- cbind(x2[,1],x3); 'send to Benika'
##write.table(new,"new_miRNA",quote=F,sep='\t',col.names = T,row.names = F);

###########EBEN:
set.seed(1);
x4 <- matrix(as.numeric(x3),nrow = nrow(x3));
#write.table(x4,'x4',col.names = F,row.names = F,sep='\t');
CV = EBelasticNet.GaussianCV(t(x4), target1, nFolds = 5,Epis = "no");
Blup1 = EBelasticNet.Gaussian(t(x4), target1,lambda = CV$Lambda_optimal,alpha = CV$Alpha_optimal, Epis = "no",verbose = 0)

####### substract the main effect:
main <- read.table("re_res_norm_main05",head=F);

x5 <- t(x4);
index_main <- main[,1];
effect_main <- main[,2];
target_new <- as.matrix(target1) - x5[,index_main] %*% (as.matrix(effect_main));
set.seed(1);
CV_epis = EBelasticNet.GaussianCV(t(x4), target_new, nFolds = 5,Epis = "yes");
Blup_epis = EBelasticNet.Gaussian(t(x4), target_new,lambda =  CV_epis$Lambda_optimal,alpha = CV_epis$Alpha_optimal, Epis = "yes",verbose = 0)

########final_run:
miRNA<-read.table('x4',header=F,stringsAsFactors = F);
mir <- as.matrix(t(x3));
mir <- matrix(as.numeric(mir),nrow = nrow(mir));
main_epi_miR_id <-  read.table("main_epi_miR",header=F);
main_epi_miR_id<- matrix(as.numeric(main_epi_miR_id),nrow = nrow(main_epi_miR_id));
new_x6 <- NULL;
for ( i in 1:nrow(main_epi_miR_id)) {
  if (main_epi_miR_id[i,1]==main_epi_miR_id[i,2]){
    new_x6 <- cbind(new_x6,mir[,main_epi_miR_id[i,1]]);
  }
  if (main_epi_miR_id[i,1]!=main_epi_miR_id[i,2]){
    col <- mir[,main_epi_miR_id[i,1]] * mir[,main_epi_miR_id[i,2]];
    new_x6 <- cbind(new_x6,col);
  }
}

new_x7 <- t(new_x6);
for( sl in 1:nrow(new_x7) ) {
  mat = matrix(as.numeric(new_x7[sl,]),1);
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(new_x7)+1));
  new_x7[sl,] = mat;
}
rm(sl, mat);

new_x8 <- t(new_x7);
set.seed(1);
CV_full = EBelasticNet.GaussianCV(new_x8, target1, nFolds = 5,Epis = "no");
Blup_full = EBelasticNet.Gaussian(new_x8, target1,lambda =  CV_full$Lambda_optimal,alpha = CV_full$Alpha_optimal, Epis = "no",verbose = 0)

final <- read.table("final_epi",header=F);
newindex <- cbind(c(1:nrow(miRNA)),rownames(miRNA));

idma <- matrix(NA,nrow = nrow(final),2);
for(i in 1:nrow(final)){
  idma[i,1] = newindex[final[i,1],2];
  idma[i,2] = newindex[final[i,2],2];
}
write.table(idma,"idma01",quote=F,sep='\t',row.names = F,col.names = F);