{
  library(parallel)
}
#=======Read files========
n1 <- readline(prompt="Enter file address: ")
no_cores<-readline(prompt="Enter number of CPUs: ")
genome<-readline(prompt="Enter the source of genome(hg19/hg18): ")
no_cores<-as.double(no_cores)
input_case<- read.csv(n1,header = T)
input_control<- read.csv('/Users/soudbeh/Desktop/Autism/COOPER_cnv_control.csv',header = T)
#=======find total number of Duplication and Deletion========
no_case_dup<- nrow(subset(input_case,input_case$cnv_type=='Dup'))
no_case_del<-nrow(subset(input_case,input_case$cnv_type=='Del'))
no_control_dup<- nrow(subset(input_control,input_control$cnv_type=='Dup'))
no_control_del<- nrow(subset(input_control,input_control$cnv_type=='Del'))
#=======hg18 and hg19==========
hg18<- data.frame(chr1=247249719,
                  chr2=242951149,
                  chr3=199501827,
                  chr4=191273063,
                  chr5=180857866,
                  chr6=170899992,
                  chr7=158821424,
                  chr8=146274826,
                  chr9=140273252,
                  chr10=135374737,
                  chr11=134452384,
                  chr12=132349534,
                  chr13=114142980,
                  chr14=106368585,
                  chr15=100338915,
                  chr16=88827254,
                  chr17=78774742,
                  chr18=76117153,
                  chr19=63811651,
                  chr20=62435964,
                  chr21=46944323,
                  chr22=49691432,
                  chrX=154913754,
                  chrY=57772954,
                  chrM=16571
)
hg19<- data.frame(chr1=249250621,
                  chr2=243199373,
                  chr3=198022430,
                  chr4=191154276,
                  chr5=180915260,
                  chr6=171115067,
                  chr7=159138663,
                  chr8=146364022,
                  chr9=141213431,
                  chr10=135534747,
                  chr11=135006516,
                  chr12=133851895,
                  chr13=115169878,
                  chr14=107349540,
                  chr15=102531392,
                  chr16=90354753,
                  chr17=81195210,
                  chr18=78077248,
                  chr19=59128983,
                  chr20=63025520,
                  chr21=48129895,
                  chr22=51304566,
                  chrX=155270560,
                  chrY=59373566,
                  chrM=16569
)
#==========seprate Dup and Del================
case_del<-subset(input_case,input_case$cnv_type=='Del')
case_dup<-subset(input_case,input_case$cnv_type=='Dup')
control_del<-subset(input_control,input_control$cnv_type=='Del')
control_dup<-subset(input_control,input_control$cnv_type=='Dup')
rm(input_case,input_control)
gc()

#==========Create Matrix=======================
if (genome == 'hg18'){
  chr1_del<-matrix(0, nrow = hg18$chr1, ncol = 6)
}else{
  chr1_del<-matrix(0, nrow = hg19$chr1, ncol = 6)
}
if (genome == 'hg18'){
  chr1_dup<-matrix(0, nrow = hg18$chr1, ncol = 6)
}else{
  chr1_dup<-matrix(0, nrow = hg19$chr1, ncol = 6)
}
chr1_case_no_del<-nrow(subset(case_del,case_del$chr=='chr1'))
chr1_control_no_del<-nrow(subset(control_del,control_del$chr=='chr1'))
chr1_case_no_dup<-nrow(subset(case_dup,case_dup$chr=='chr1'))
chr1_control_no_dup<-nrow(subset(control_dup,control_dup$chr=='chr1'))
#==========Define function count===============
{
  count_del1<-function(x){
    chr1_del[x,1]<<-x
  }
  count_dup1<-function(x){
    chr1_dup[x,1]<<-x
  }
  count_del2<-function(x){
    chr1_del[x,2]<<-nrow(subset(case_del,case_del$start<=x & case_del$end >= x))
  }
  count_dup2<-function(x){
    chr1_dup[x,2]<<-nrow(subset(case_dup,case_dup$start<=x & case_dup$end >= x))
  }
  count_del3<-function(x){
    chr1_del[x,3]<<-nrow(subset(control_del,control_del$start<=x & control_del$end >= x))
  }
  count_dup3<-function(x){
    chr1_dup[x,3]<<-nrow(subset(control_dup,control_dup$start<=x & control_dup$end >= x))
  }
  count_del4<-function(x){
    fisher_matrix <<- matrix(c(chr1_del[x,2],chr1_del[x,3],chr1_case_no_del- chr1_del[x,2],chr1_control_no_del -chr1_del[x,3]),nrow = 2,ncol = 2)
    test <<- fisher.test(fisher_matrix,alternative = "two.sided",conf.level =0.9)
    chr1_del[x,4]<<- test$p.value
  }
  count_dup4<-function(x){
    fisher_matrix <<- matrix(c(chr1_dup[x,2],chr1_dup[x,3],chr1_case_no_dup- chr1_dup[x,2],chr1_control_no_dup -chr1_dup[x,3]),nrow = 2,ncol = 2)
    test <<-fisher.test(fisher_matrix,alternative = "two.sided",conf.level =0.9)
    chr1_dup[x,4]<<- test$p.value
  }
  count_del5<-function(x){
    fisher_matrix <<- matrix(c(chr1_del[x,2],chr1_del[x,3],chr1_case_no_del- chr1_del[x,2],chr1_control_no_del -chr1_del[x,3]),nrow = 2,ncol = 2)
    test<<-fisher.test(fisher_matrix,alternative = "greater",conf.level =0.9)
    chr1_del[x,5]<<- test$p.value
  }
  count_dup5<-function(x){
    fisher_matrix <<- matrix(c(chr1_dup[x,2],chr1_dup[x,3],chr1_case_no_dup- chr1_dup[x,2],chr1_control_no_dup -chr1_dup[x,3]),nrow = 2,ncol = 2)
    test<<-fisher.test(fisher_matrix,alternative = "greater",conf.level =0.9)
    chr1_dup[x,5]<<-test$p.value
  }
  count_del6<-function(x){
    fisher_matrix <<- matrix(c(chr1_del[x,2],chr1_del[x,3],chr1_case_no_del- chr1_del[x,2],chr1_control_no_del -chr1_del[x,3]),nrow = 2,ncol = 2)
    test<-fisher.test(fisher_matrix,alternative = "less",conf.level =0.9)
    chr1_del[x,6]<<- test$p.value
  }
  count_dup6<-function(x){
    fisher_matrix <<- matrix(c(chr1_dup[x,2],chr1_dup[x,3],chr1_case_no_dup- chr1_dup[x,2],chr1_control_no_dup -chr1_dup[x,3]),nrow = 2,ncol = 2)
    test<<-fisher.test(fisher_matrix,alternative = "less",conf.level =0.9)
    chr1_dup[x,6]<<- test$p.value
  }
}

#============chromosome count==================

cl <- makeCluster(no_cores)

if (genome == 'hg18'){
  t<-hg18$chr1
  clusterExport(cl=cl, varlist=c('chr1_control_no_dup','chr1_case_no_dup','chr1_control_no_del','chr1_case_no_del','control_del','control_dup','case_del','case_dup',"count_del1", 'count_del2','count_del3','count_del4','count_del5','count_del6','count_dup1','count_dup2','count_dup3','count_dup4','count_dup5','count_dup6','chr1_dup','chr1_del'))
  save<-parLapply(cl, 1:t, count_del1)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,1]<-save
  save<-parLapply(cl, 1:t, count_del2)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,2]<-save
  save<-parLapply(cl, 1:t, count_del3)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,3]<-save
  save<-parLapply(cl, 1:t, count_del4)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,4]<-save
  save<-parLapply(cl, 1:t, count_del5)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,5]<-save
  save<-parLapply(cl, 1:t, count_del6)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,6]<-save
  save<-parLapply(cl, 1:t, count_dup1)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,1]<-save
  save<-parLapply(cl, 1:t, count_dup2)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,2]<-save
  save<-parLapply(cl, 1:t, count_dup3)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,3]<-save
  save<-parLapply(cl, 1:t, count_dup4)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,4]<-save
  save<-parLapply(cl, 1:t, count_dup5)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,5]<-save
  save<-parLapply(cl, 1:t, count_dup6)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,6]<-save
}else{
  t<-hg19$chr1
  clusterExport(cl=cl, varlist=c('chr1_control_no_dup','chr1_case_no_dup','chr1_control_no_del','chr1_case_no_del','control_del','control_dup','case_del','case_dup',"count_del1", 'count_del2','count_del3','count_del4','count_del5','count_del6','count_dup1','count_dup2','count_dup3','count_dup4','count_dup5','count_dup6','chr1_dup','chr1_del'))
  save<-parLapply(cl, 1:t, count_del1)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,1]<-save
  save<-parLapply(cl, 1:t, count_del2)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,2]<-save
  save<-parLapply(cl, 1:t, count_del3)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,3]<-save
  save<-parLapply(cl, 1:t, count_del4)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,4]<-save
  save<-parLapply(cl, 1:t, count_del5)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,5]<-save
  save<-parLapply(cl, 1:t, count_del6)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_del[,6]<-save
  save<-parLapply(cl, 1:t, count_dup1)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,1]<-save
  save<-parLapply(cl, 1:t, count_dup2)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,2]<-save
  save<-parLapply(cl, 1:t, count_dup3)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,3]<-save
  save<-parLapply(cl, 1:t, count_dup4)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,4]<-save
  save<-parLapply(cl, 1:t, count_dup5)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,5]<-save
  save<-parLapply(cl, 1:t, count_dup6)
  save<- as.data.frame(save)
  save<-t(save)
  chr1_dup[,6]<-save
}


