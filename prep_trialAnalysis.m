function [] = prep_trialAnalysis(idx_file)
%PREP_TRIALANALYSIS prepares data for a trial based analysis
%   Data is split into trials, artifacts get removed
%   Output is cleaned data split into trials of one file (normally 4 blocks)
%   Saved data structure contains also more information like response time,
%   hazard rate, displayed samples etc.

        
    % Load important information about files
    % files contains the complete names of all files that must be processed
    % info_EL_blocks: matrix with 5 columns
    % col 1 - first block to process with regard to behavioral data (normally 1 or 5)
    % col 2 - last block to process with regard to behavioral data (normally 4 or 8)
    % col 3 - number of blocks in file (normally 4)
    % col 4 - run (1 or 2)
    % col 5 - eyelink data (1: usable, 0: use veog instead)
    [~,files] = xlsread('/mnt/homes/home024/jschipp/Surprise_Drug/meg_data/Info_filewise.xlsx');
    info_EL_blocks = xlsread('/mnt/homes/home024/jschipp/Surprise_Drug/meg_data/Info_filewise.xlsx');
    filein = files{idx_file};

    ID = [filein(1:5) filein(end-5:end-3)]; % Subject ID + Session number + file number
    if regexp(ID, 'URG_S*') % With this subject the session number was missing when registered
        ID = ['URG-1' filein(end-5:end-3)];
    end
    
    addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
    addpath '/mnt/homes/home024/jschipp/Surprise_Drug/meg_preprocessing'
    addpath '/mnt/homes/home024/jschipp/Surprise_Drug/meg_data/'
    ft_defaults

    fprintf('\n\n ---------------- \n PROCESSING FILE %s...\n ---------------- \n\n', ID)
    
    %% Construct trial matrix
    
    % =====================================================================
    %   1   Construct trial matrix
    % =====================================================================

    if idx_file < 229 % Older datasets (Caro)
        cd /mnt/homes/home024/pmurphy/meg_data/surpriseD/
    else % Newer datasets (Julia)
        cd /mnt/homes/home024/jschipp/Surprise_Drug/meg_data/
    end

    fprintf('\n ---------- 1 Constructing trial matrix ----------\n')
    
    % Specify cfg
    cfgin = [];
    cfgin.idx = idx_file; 
    cfgin.ID = ID;
    cfgin.subj = ID(1:3);
    cfgin.session = ID(5);
    cfgin.fileNum = ID(7:8);
    cfgin.dataset = filein;
    cfgin.trialdef.prestim = 1; % 1 second offset for TF baseline
    cfgin.trialdef.poststim = 0.5; % add 0.5 s after response cue
    cfgin.trialfun = 'trialbasedfun_surpriseD';
    cfgin = ft_definetrial(cfgin);
    
    % Needed cfgin stored in cfgin.trl
    cfgtrial = cfgin.trl;
    clear cfgin
    cfgtrial.alltrl = cfgtrial.trl;
    cfgtrial = rmfield(cfgtrial,'trl');
    
    % Downsample originally created trial matrix
    cfgtrial.alltrl(:,1:3) = round(cfgtrial.alltrl(:,1:3)/1200*400);
    
    % Downsample the block bounds
    cfgtrial.blockBound_trl = round(cfgtrial.blockBound_trl/1200*400);

    %% Loop through blocks

    firstBl = 1;
    
    if strcmp(ID,'UDK-1_03') % 5 block contained, 4th is empty gets removed
        loop_vec = 1:4;
    else
        loop_vec = firstBl:info_EL_blocks(idx_file,3);
    end
          
    for block = loop_vec

        % Consider the actual block numer for behavioral data (1-8)
        actualblock = (block-1) + info_EL_blocks(idx_file,1);
        
       if strcmp(ID,'UDK-1_03') && block == 4
            actualblock = actualblock+1;
       elseif strcmp(ID,'UOC-1_01') 
            actualblock = actualblock+2;
       elseif strcmp(ID,'UOC-1_03') 
            actualblock = actualblock+4;
        end

        fprintf('\n\n ---------------- \n Loop through BLOCK #%d...\n ---------------- \n\n', block)

        % =====================================================================
        %   2   DEFINE BLOCK AND LOAD RELEVANT CHANNELS 
        % =====================================================================

        fprintf('\n ---------- 2 Defining block and relevant channels ----------\n')

        % Specify cfg
        cfgbl = [];
        cfgbl.dataset = filein;
        cfgbl.block = block;
        cfgbl.ID = ID;
        cfgbl.trialdef.prestim = 10; % Add 10 seconds before and after the block
        cfgbl.trialdef.poststim = 10;
        cfgbl.trialfun = 'trialfun_surpriseD_continuous';
        cfgbl = ft_definetrial(cfgbl);
        cfgbl.continuous = 'yes'; % read in data as continuous
        cfgbl.channel = {'meg','EEG001','EEG002','EEG003','EEG057','EEG058','HLC*','UADC*'};
        data = ft_preprocessing(cfgbl);
        fsample_old = data.fsample;
        sampleinfo_old = data.sampleinfo; 
        
        % =====================================================================
        %   3    REMOVE LINE NOISE
        % =====================================================================

        fprintf('\n ---------- 3 Filtering out line noise (50 Hz) and its harmonics -----------\n')

        cfg             = [];
        cfg.bsfilter    = 'yes';
        cfg.bsfreq      = [49 51; 99 101; 149 151];
        data            = ft_preprocessing(cfg, data);
        
        % =====================================================================
        %   4   HIGH PASS FILTER (cutoff 0.1 Hz)
        % =====================================================================

        fprintf('\n ---------- 4 High pass filtering (cutoff 0.1 Hz) -----------\n')

        cfg          = [];
        cfg.hpfilter = 'yes';
        cfg.hpfreq   = 0.1;  
        cfg.hpfiltord = 3;
        cfg.hpfilttyoe = 'firws';
        data = ft_preprocessing(cfg,data);
        
        % =====================================================================
        %   5    DOWN SAMPLE DATA TO 400 HZ
        % =====================================================================

        fprintf('\n ---------- 5 Resampling block -----------\n')

        cfgres.resample = 'yes';
        cfgres.fsample = 1200;
        cfgres.resamplefs = 400;
        cfg.detrend = 'no';
        data = ft_resampledata(cfgres, data);
        sampleinfo_new = [round(sampleinfo_old(1))./fsample_old.*cfgres.resamplefs...
                          round(sampleinfo_old(1))./fsample_old.*cfgres.resamplefs+length(data.time{1}-1)]; 

        % =====================================================================
        %   6    LOAD ICA weights and IDs of artifactual components
        % =====================================================================

        fprintf('\n ---------- 6 Load ICA weights and IDs of artifactual components -----------\n')
        name = ['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/comp_ICA/' ID(1:3)];
        load([name '/comp_' ID '.mat' ]);
        load([name '/comp2rej_' ID '.mat' ]);
        
        % Remove artifactual components
        cfg = [];
        cfg.component = rejComps;
        data_cl = ft_rejectcomponent(cfg, comp, data);
        
        % Separate current block ! ACTUALBLOCK CORRECT?
        cfgtrial.trl = cfgtrial.alltrl(cfgtrial.alltrl(:,5)==actualblock,:);
        
        % Substract block onset from all sample numbers (except for offset) in the trial matrix
        cfgtrial.trl(:,1:2) = cfgtrial.trl(:,1:2) - cfgtrial.blockBound_trl(block,1);
        
        % Segment block's data into trials
        trials = ft_redefinetrial(cfgtrial, data_cl);
        
        % Load previously created artifact matrices (aligned to block onset!)
        name = ['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/preICA_artifactMatrices/' ID(1:3)];
        if strcmp(ID,'UDK-1_03') && block == 4 % block 4 was not started on purpose so in this script block for corresponds to block 5 in meg
            % but to behavior block 5, Since the 4th block meg block was
            % only removed in this script so it corresponds still to
            % artifact matrix block 5 ... confusing but corrected!
            artifacts = ['Artifacts_' ID '_Block_' num2str(5)];
        else
            artifacts = ['Artifacts_' ID '_Block_' num2str(block)];
        end
        load([name '/' artifacts]);
        
        % Mark trials that overlap with an artifact
        trls2remove = [];
        for j = 1:length(trials.sampleinfo)
            % Mark trials with head movements
            for i = 1:size(artifact_headM,1)
                if trials.sampleinfo(j,1) <= artifact_headM(i,1) && trials.sampleinfo(j,2) >= artifact_headM(i,2)
                    trls2remove = [trls2remove j];
                end
            end
            % Mark trials with jumps
            for i = 1:size(artifact_Jump,1)
                if trials.sampleinfo(j,1) <= artifact_Jump(i,1) && trials.sampleinfo(j,2) >= artifact_Jump(i,2)
                    trls2remove = [trls2remove j];
                end
            end
            % Mark trials with muscle artifacts
            for i = 1:size(artifact_Muscle,1)
                if trials.sampleinfo(j,1) <= artifact_Muscle(i,1) && trials.sampleinfo(j,2) >= artifact_Muscle(i,2)
                    trls2remove = [trls2remove j];
                end
            end
            % Mark trials with saccades
            for i = 1:size(artifact_saccade,1)
                if trials.sampleinfo(j,1) <= artifact_saccade(i,1) && trials.sampleinfo(j,2) >= artifact_saccade(i,2)
                    trls2remove = [trls2remove j];
                end
            end
            % Mark trials with blinks
            for i = 1:size(artifact_Ygaze,1)
                if trials.sampleinfo(j,1) <= artifact_Ygaze(i,1) && trials.sampleinfo(j,2) >= artifact_Ygaze(i,2)
                    trls2remove = [trls2remove j];
                end
            end
            % Mark trials with metal artefacts
            for i = 1:size(artifact_metal,1)
                if trials.sampleinfo(j,1) <= artifact_metal(i,1) && trials.sampleinfo(j,2) >= artifact_metal(i,2)
                    trls2remove = [trls2remove j];
                end
            end
        end
        
        % Some trials needed to be removed because of multiple artifacts
        trls2remove = unique(trls2remove); % Remove douplicates
        numtrial = 1:length(trials.sampleinfo);
        trls2keep = numtrial(~ismember(numtrial,trls2remove));
        
        cfg = [];
        cfg.trials = trls2keep;
        % Remove these trials
        trials = ft_selectdata(cfg, trials);
        
        % Concatenate blocks
        if block == 1 
            old_trials = trials;
        else
            cfg = [];
            old_trials = ft_appenddata(cfg, old_trials, trials);
        end
    end
    
    % Define path/folder to save cleaned data
    if strcmp(ID(1:5),'JBK-1')
        ID(1:5) = 'JPK-1';
    end
    mat_name = ['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/data_trials/' ID(1:3) '/'];
    
        if 7==exist(mat_name,'dir')
            cd(mat_name)
        else
            mkdir(mat_name)
            cd(mat_name)
        end
    
    % Rename
    all_trials_cl = old_trials;
    
    all_trials_cl.trialInfoLabel = cfgtrial.trialInfoLabel(:,4:11);
    all_trials_cl.drug = cfgtrial.drug;
    all_trials_cl.HR = cfgtrial.HR;
    all_trials_cl.blockBounds = cfgtrial.blockBound_trl;
    
    clear old_trials;
    
    % Save concatenated data file
    datastore = ['data_clean_postICA_' ID];
    save(datastore,'all_trials_cl','-v7.3');

end

