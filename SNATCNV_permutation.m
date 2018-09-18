clear;
%%%%%%%%%%%%%%%%%%%%%%%%% SETTING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
number_CASE_DEL=12721;
number_CASE_DUP=4441;
number_CONTROL_DEL=11924;
number_CONTROL_DUP=5015;
number_case=number_CASE_DUP; %specify this parameter with deletion of duplication values
number_control=number_CONTROL_DUP; %specify this parameter with deletion of duplication values
resultname='../Permutation_ASD_dup_100k.txt';
max_cnv=372; %length(CNVnumber_org);
%%%%%%%%%%%%%%%%%%%%%%%%% SETTING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%

CNVarray_report(100001,max_cnv)=0.000;
for i =1:max_cnv
CNVarray_report(1,i)=i;
end

ASD_random_temp=cell(number_case+number_control,2);
ASD_random_temp(1:number_case,1)=cellstr('case');
ASD_random_temp(number_case+1:number_case+number_control,1)=cellstr('control');

for counter = 1: max_cnv
    numberCNV=CNVarray_report(1,counter);
    ASD_random_temp(1:number_case+number_control,2)=num2cell(rand(number_case+number_control,1));
    ASD_random_temp=sortrows(ASD_random_temp,+2);
    for i = 2: 100001
        fprintf(' %1.0f , %1.0f\n',numberCNV,i);
        r_case_control = randi([1 number_case+number_control],numberCNV,1);
        CNV_case_positive=nnz(find(strcmp(ASD_random_temp(r_case_control,1), cellstr('case'))));
        CNV_control_positive=nnz(find(strcmp(ASD_random_temp(r_case_control,1), cellstr('control'))));
        CNV_case_negative=number_case - CNV_case_positive;
        CNV_control_negative=number_control - CNV_control_positive;
        
        contingency_table = [CNV_case_positive, CNV_case_negative; CNV_control_positive, CNV_control_negative];
        [~,p_value_right]= fishertest(contingency_table, 'Tail','right', 'Alpha', 0.05); 
        CNVarray_report(i,counter)=p_value_right;
    end
end

%%% export permutation
fprintf('exporting permutation result...\n')
dlmwrite(resultname,CNVarray_report,'delimiter', '\t');
