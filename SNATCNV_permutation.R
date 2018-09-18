#=============Permutation================
number_of_permutation<- 500000
#max_cnv<- max(chr2_del[,2],chr2_del[,3],chr2_dup[,2],chr2_dup[,3])
#no_case_dup<- nrow(subset(input_case,input_case$CNV.Type=='Del'))
#no_control_dup<- nrow(subset(input_control,input_control$CNV.Type=='Del'))su:774-nu:1081
max_cnv<-179
no_case_del<-785
no_control_del<-1094
ASD_random_del_chr2<- matrix(, nrow =no_case_del+no_control_del , ncol = 2)
ASD_random_del_chr2[1:no_case_del,1]<-'case'
tt<-no_case_del+1
ttt<-nrow(ASD_random_del_chr2)
print(ttt)
ASD_random_del_chr2[tt:ttt,1]<-'control'
CNVarray_report <- matrix(1,nrow=number_of_permutation, ncol = max_cnv)

for (numberCNV in 1:max_cnv) {
  ASD_random_del_chr2[,2]<- runif(ttt,0,1)
  ASD_random_del_chr2<-ASD_random_del_chr2[order(ASD_random_del_chr2[,2]),]
  print(paste0(numberCNV,'th permutation round \n'))
  for (i in 1:number_of_permutation){
    r_case_control<- sample(2:ttt, 1)
    ASD<-ASD_random_del_chr2[1:r_case_control,]
    CNV_case_positive<- nrow(subset(ASD,ASD[,1]=='case'))
    CNV_control_positive <-nrow(subset(ASD,ASD[,1]=='control'))
    CNV_case_negative <- no_case_del -CNV_case_positive
    CNV_control_negative <- no_control_del - CNV_control_positive
    contingency_table<-matrix(c(CNV_case_positive,CNV_control_positive,CNV_case_negative,CNV_control_negative),nrow = 2)
    test<-fisher.test(contingency_table,alternative = "greater",conf.level =0.9)
    pval<-test$p.value
    CNVarray_report[i,numberCNV]<-pval
  }
}
print('write')
path<-'/home/san/halinejad/Desktop/Dashti'
path<-paste0(path,'/CNVarray_report_del.cnv')
write.csv(CNVarray_report,path)
