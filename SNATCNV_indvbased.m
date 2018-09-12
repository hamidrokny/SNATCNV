clear;
%%% GENOME_BUILD initialization
fprintf('initializing...\n')
chr_initialization;
%%%%%%%%%%%%%%%%%%%%%%%%% SETTING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
CNV_type='Del';
GENOME_BUILD=GENOME_BUILD_19; %genome build
PATH_INPUT='../input';
PATH_OUTPUT='../output';
FILE_CASE_NAME='ASD_cnv_case.txt';
FILE_CONTROL_NAME='ASD_control_case.txt';
NUMBER_OF_CHR=24;
%%%%%%%%%%%%%%%%%%%%%%%%% SETTING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%

%%% import case CNV
fname=strcat(PATH_INPUT,'/',FILE_CASE_NAME);
formatString='%s %s %s %s %s %s %s %*[^\n]';
fid=fopen(fname, 'r');
CASE_cnv=textscan(fid, formatString, 'delimiter', '\t');
fclose(fid);
%%% import control CNV
fname=strcat(PATH_INPUT,'/',FILE_CASE_NAME);
formatString='%s %s %s %s %s %s %s %*[^\n]';
fid=fopen(fname, 'r');
CONTROL_cnv=textscan(fid, formatString, 'delimiter', '\t');
fclose(fid);

fname=strcat(PATH_OUTPUT,'/','significant_regions_based_on_indv_',CNV_type);

regions_report(1,8)=0;
row_counter=0;
for chr_id=22 : NUMBER_OF_CHR
    %%% count number of case and control
    indx=strcmp(CASE_cnv{1,6}, CNV_type);
    number_case=length(unique(CASE_cnv{1,5}(indx))); % indv based
    indx=strcmp(CONTROL_cnv{1,6}, CNV_type);
    number_control=length(unique(CONTROL_cnv{1,5}(indx))); % indv based

    CNVarray_case=0;
    CNVarray_control=0;
    maskarray=0;
    CNVarray_case(cell2mat(GENOME_BUILD(chr_id,3)),1)=0;
    CNVarray_control(cell2mat(GENOME_BUILD(chr_id,3)),1)=0;
    maskarray(cell2mat(GENOME_BUILD(chr_id,3)),1)=0;

    %%% filter on case chr
    indx=strcmp(CASE_cnv{1,1}, GENOME_BUILD(chr_id));  
    temp = cellfun(@(x) x(indx), CASE_cnv, 'UniformOutput', false);
    %%% filter on case type
    indx=strcmp(temp{1,6}, CNV_type);
    number_chr_type=nnz(indx);
    SCH_chr_type = cellfun(@(x) x(indx), temp, 'UniformOutput', false);
    %%% count number of overlapping CNVs 
    fprintf('chr %1.0f - counting overlapping CNVs in case...\n', chr_id)
    for i = 1 :  number_chr_type
        maskarray(:,1)=0;
        mask_start=str2double(cell2mat(SCH_chr_type{1,3}(i)));
        mask_end=str2double(cell2mat(SCH_chr_type{1,4}(i)));
        maskarray(mask_start:mask_end)=1;
        CNVarray_case(mask_start:mask_end)= CNVarray_case(mask_start:mask_end)+ maskarray(mask_start:mask_end);
    end

    %%% filter of control chr
    indx=strcmp(CONTROL_cnv{1,1}, GENOME_BUILD(chr_id));
    temp = cellfun(@(x) x(indx), CONTROL_cnv, 'UniformOutput', false);   
    %%% filter on control Dup
    indx=strcmp(temp{1,6}, CNV_type);
    number_chr_type=nnz(indx);
    SCH_chr_type = cellfun(@(x) x(indx), temp, 'UniformOutput', false);  
    %%% count number of overlapping CNVs
    fprintf('chr %1.0f - counting overlapping CNVs in control...\n', chr_id)
    for i = 1 :  number_chr_type
        maskarray(:,1)=0;
        mask_start=str2double(cell2mat(SCH_chr_type{1,3}(i)));
        mask_end=str2double(cell2mat(SCH_chr_type{1,4}(i)));
        maskarray(mask_start:mask_end)=1;
        CNVarray_control(mask_start:mask_end)= CNVarray_control(mask_start:mask_end)+ maskarray(mask_start:mask_end);
    end

    maskarray=0;
    diff_row=0;
    fprintf('calculating P value...\n')
    diff_case=find(diff(CNVarray_case) ~= 0);
    diff_control=find(diff(CNVarray_control) ~= 0);
    %%% merging two array
    [sA1, sA2] = size(diff_case);
    [sB1, sB2] = size(diff_control);
    diff_row(sA1+1:sA1+sB1,1:sB2)=diff_control;
    diff_row(1:sA1, 1:sA2)= diff_case;
    diff_row=unique(diff_row);
    number_of_position=length(diff_row);

    CNVarray_report=0;
    CNVarray_report(number_of_position+1,8)=0;

    %%% p-value calculation - start %%%
    region_start=1;
    for i = 1 : number_of_position
        fprintf('%1.0f , %1.0f of %1.0f \n',chr_id, i, number_of_position);
        region_end=diff_row(i,1);
    
        CNV_case_positive=CNVarray_case(region_start);
        CNV_case_negative=number_case - CNV_case_positive;
        CNV_control_positive=CNVarray_control(region_start);
        CNV_control_negative=number_control - CNV_control_positive;  
        contingency_table = [CNV_case_positive, CNV_case_negative; CNV_control_positive, CNV_control_negative];
    
        [~,p_value_right]= fishertest(contingency_table, 'Tail','right', 'Alpha', 0.05); 
        [~,p_value_both]= fishertest(contingency_table, 'Tail','both', 'Alpha', 0.05); 
        [~,p_value_left]= fishertest(contingency_table, 'Tail','left', 'Alpha', 0.05); 
    
        CNVarray_report(i,1)=chr_id;
        CNVarray_report(i,2)=region_start;
        CNVarray_report(i,3)=region_end;
        CNVarray_report(i,4)=CNV_case_positive;
        CNVarray_report(i,5)=CNV_control_positive;
        CNVarray_report(i,6)=p_value_right;
        CNVarray_report(i,7)=p_value_left;
        CNVarray_report(i,8)=p_value_both;
    
        region_start=diff_row(i,1)+1;
    end
    %%% for the last region
    region_start=diff_row(i,1)+1;
    region_end=cell2mat(GENOME_BUILD(chr_id,3)); %end position in the chr

    CNV_case_positive=CNVarray_case(region_start);
    CNV_case_negative=number_case - CNV_case_positive;
    CNV_control_positive=CNVarray_control(region_start);
    CNV_control_negative=number_control - CNV_control_positive;  
    contingency_table = [CNV_case_positive, CNV_case_negative; CNV_control_positive, CNV_control_negative];
    
    [~,p_value_right]= fishertest(contingency_table, 'Tail','right', 'Alpha', 0.05); 
    [~,p_value_both]= fishertest(contingency_table, 'Tail','both', 'Alpha', 0.05); 
    [~,p_value_left]= fishertest(contingency_table, 'Tail','left', 'Alpha', 0.05); 
    
    CNVarray_report(i+1,1)=chr_id;
    CNVarray_report(i+1,2)=region_start;
    CNVarray_report(i+1,3)=region_end;
    CNVarray_report(i+1,4)=CNV_case_positive;
    CNVarray_report(i+1,5)=CNV_control_positive;
    CNVarray_report(i+1,6)=p_value_right;
    CNVarray_report(i+1,7)=p_value_left;
    CNVarray_report(i+1,8)=p_value_both;
    %%% p-value calculation - end %%%

    regions_report(row_counter+1:row_counter+length(CNVarray_report),:)=CNVarray_report;
    row_counter=length(regions_report);
    
end

%%% export association
fprintf('exporting associations...\n')
fid=fopen(fname, 'a+');
for i= 1 : length(regions_report)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%2.3e\t%2.3e\t%2.3e\n', regions_report(i,:));
end
fclose(fid);

