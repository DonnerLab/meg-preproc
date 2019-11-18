clear, close all

addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
addpath '/mnt/homes/home024/jschipp/Surprise_Drug/meg_preprocessing'
addpath '/mnt/homes/home024/jschipp/Surprise_Drug/meg_data/'
ft_defaults

% Go to data folder
cd /mnt/homes/home024/jschipp/Surprise_Drug/meg_data/

% Load important information about files
% files contains the complete names of all files that must be processed
% info_EL_blocks: matrix with 5 columns
% col 1 - first block to process with regard to behavioral data (normally 1 or 5)
% col 2 - last block to process with regard to behavioral data (normally 4 or 8)
% col 3 - number of blocks in file (normally 4)
% col 4 - run (1 or 2)
% col 5 - eyelink data (1: usable, 0: use veog instead)
%test = readtable('Info_filewise.xlsx')
[~,files] = xlsread('Info_filewise');
[~,txt,~] = xlsread('/mnt/homes/home024/jschipp/Surprise_Drug/Manual_Comp_Rej.xlsx');
txt = txt(:,2:3);

for i = 358:361
%for i = 1:length(txt)   
    
    filein = files{i};
    
    %filein = 'CTG-2_Surprise_20180807_02.ds';
    ID = [filein(1:5) filein(end-5:end-3)]; % Subject ID + Session number + file number
    if regexp(ID, 'URG_S*') % With this subject the session number was missing when registered
        ID = ['URG-1' filein(end-5:end-3)];
    end
    

    rej = txt(i,:);
    rej = [rej{1,1} rej{1,2}];
    rejComps = [];
    
    remain = rej;
    while ~isempty(remain)
        [token,remain] = strtok(remain, '; ');
        token = str2double(token);
        rejComps = [rejComps; token];
    end
    
    rejComps = unique(rejComps)';
    rejComps(isnan(rejComps)) = [];
    
    
    cd(['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/comp_ICA/' ID(1:3) '/'])
    comp2rej = ['comp2rej_' ID '.mat'];
    save(comp2rej,'rejComps')
    
end



