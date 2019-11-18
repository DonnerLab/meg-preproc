% Calculate mean and std of rejected components

[~,drug_H_info,~] = xlsread('/mnt/homes/home024/jschipp/Surprise_Drug/drug_H_info.xlsx');
subjects = drug_H_info(:,1);

count = 1;
for s = 1:length(subjects)
    file_names = dir(['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/comp_ICA/' subjects{s} '/*rej*']);
    
    for f = 1:length(file_names)
        load(['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/comp_ICA/' subjects{s} '/' file_names(f).name])
        comp_num(count) = length(rejComps);
        count = count+1;
    end
    
end

mean_rej = mean(comp_num);
std_rej = std(comp_num);
  
  

