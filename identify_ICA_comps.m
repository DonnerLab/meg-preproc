function [] = identify_ICA_comps(idx_file)
% IDENTIFY ICA COMPONENTS FOR BLINKS, SACCADES AND HEART BEAT
%   The cleaned data (head movements, muscle artifacts and squid jumps are removed)
%   is first segmented around artifacts of interest, the component
%   time-series are calculated, the frequency decomposition and coherence
%   is calculated and finally the timelocked time-series of the components
%   with the highest coherence values (top 25) and their topographies are
%   plotted.

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
    info_EL_blocks = xlsread('Info_filewise');
    filein = files{idx_file};

    %filein = 'CTG-2_Surprise_20180807_02.ds';
    ID = [filein(1:5) filein(end-5:end-3)]; % Subject ID + Session number + file number
    if regexp(ID, 'URG_S*') % With this subject the session number was missing when registered
        ID = ['URG-1' filein(end-5:end-3)];
    end
    
    % Folder containing data 
    cd(['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/data_clean/' ID(1:3)]) 
    load(sprintf('all_artifacts_%s.mat', ID))   % Load artifacts 
    
    if isempty(artifact_saccade)
        % Only perform ICA for blinks and ecg
        artifact_type = {'blink', 'ecg'};
        chanRej = {'EEG057','EEG058'};
        fprintf('\n +++ NOTE: NO Saccades +++\n')
    else
        % Perform ICA for VEOG, YGAZE, ECG
        artifact_type = {'blink', 'saccade', 'ecg'};
        chanRej = {'EEG057','UADC004','EEG058'};
    end
    
    blinkwin = [-1 2];
    
    
    %% Loop through artifacts (blinks, saccades, ecg)
    for n = 1:length(artifact_type)
        fprintf(['\n ========== PROCESSING ' artifact_type{n} ' =========\n'])
       
        cd(['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/data_clean/' ID(1:3)])    
        load(sprintf('data_clean_%s.mat',ID))       % Load cleaned data
        artifact_seg = [];
        
        % =====================================================================
        %   1   SEGMENT DATA AROUND ARTIFACTS OF INTEREST
        % =====================================================================

        fprintf(['\n ---------- 1 Segmenting data around (' artifact_type{n} ') artifacts ----------\n'])

        if strcmp(artifact_type{n}, 'ecg')
            % Segment data around artifacts of interest (ft_artifact_ecg)
            
            % Define relevant channel
            cfg = [];
            cfg.channel = chanRej{n};
            if strcmp(ID,'TMW-1_04')
                cfg.latency = [0 175]; 
            elseif strcmp(ID, 'FLW-1_01')
                cfg.latency = [725 2725];
            elseif strcmp(ID, 'CTG-2_03')
                cfg.latency = [750 2460];
            end
            art = ft_selectdata(cfg, data);
            art.label{:} = 'ECG';
            
            % Get artifact matrices Mx2 (start and end sample)
            cfg                       = [];
            cfg.dataset               = data;
            cfg.continuous            = 'yes';
            cfg.channel               = chanRej{n};
            cfg.artfctdef.ecg.inspect = chanRej{n};
            cfg.artfctdef.ecg.feedback = 'no';
            cfg.artfctdef.ecg.pretim  = 0.3;
            cfg.artfctdef.ecg.psttim  = 0.3;
            [art, artifact_ecg]        = ft_artifact_ecg(cfg, art);
            
            % Segment data around artifacts of interest
            cfg            = [];
            cfg.trl        = [artifact_ecg,zeros(size(artifact_ecg,1),1)];
            
            if strcmp(ID,'FLW-1_01') | strcmp(ID,'CTG-2_03')
                window_ecg = [3e5 10e5];
            elseif strcmp(ID,'TMW-1_04')
                window_ecg = [1 0.7e5];
            end
                
            if exist('window_ecg')
                trl2k = [];

                for i = 1:length(cfg.trl)
                    if cfg.trl(i,1)>window_ecg(1) && cfg.trl(i,2)<window_ecg(2)
                        trl2k = [trl2k; i];
                    end
                end
                cfg.trl = cfg.trl(trl2k,:);
            end
            
%             cfg.channel    = {'MEG'};
%             cfg.continuous = 'yes';
%             data_art       = ft_preprocessing(cfg,data); % Only MEG channels
%             data_art       = ft_redefinetrial(cfg,data_art); % Segmenting MEG channel
%             
%             % Segment ecg channel
%             cfg.channel    = chanRej{n};
%             art            = ft_preprocessing(cfg,data); % Only ECG channel
%             art            = ft_redefinetrial(cfg,art); % Segmenting ECG channel
%             art.channel{:} = chanRej{n};
%             
        elseif strcmp(artifact_type{n}, 'blink') % Adjust window size (same for all artifacts)
            artifact_seg = [artifact_blink(:,1)+blinkwin(1)*data.fsample...
                artifact_blink(:,1)+blinkwin(2)*data.fsample];
            
            if artifact_seg(1,1) < 0
                s = [1 1201];
                artifact_seg(1,:) = s;
            end
            artifact_blink = artifact_seg;
            
            cfg.trl        = [artifact_blink, zeros(size(artifact_blink,1),1)];
        elseif strcmp(artifact_type{n}, 'saccade') % Adjust window size (same for all artifacts)
            artifact_seg = [artifact_saccade(:,1)+blinkwin(1)*data.fsample...
                artifact_saccade(:,1)+blinkwin(2)*data.fsample];
            
            if artifact_seg(1,1) < 0
                s = [1 1201];
                artifact_seg(1,:) = s;
            end
            artifact_saccade = artifact_seg;
            
            cfg.trl = [artifact_saccade, zeros(size(artifact_saccade,1),1)];
        end
            
        cfg.channel    = {'MEG'};
        cfg.continuous = 'yes';
        if strcmp(ID,'TMW-1_04')
                cfg.latency = [0 175]; 
            elseif strcmp(ID, 'FLW-1_01')
                cfg.latency = [725 2725];
            elseif strcmp(ID, 'CTG-2_03')
                cfg.latency = [750 2460];
            end
        data_art       = ft_preprocessing(cfg,data);
        data_art       = ft_redefinetrial(cfg,data_art);
        cfg.channel    = {chanRej{n}};
        art            = ft_preprocessing(cfg,data);
        art            = ft_redefinetrial(cfg,art);
        art.channel{:} = chanRej{n};
        
        % The first and/or the last segments could potentially contain NaN
        % because of the window (-1s and 2s after the blink/saccade)
        % remove those artifacts/segments in order to avoid issues 
        
        % check if nans are in a segemnt/trial
        seg2rej = ones(size(data_art.trial,2),1);
        for ii = 1:size(data_art.trial,2)
            if any(isnan(data_art.trial{1,ii}(:)))
                seg2rej(ii) = 0;
            end
        end
        
        cfg = [];
        cfg.trials = logical(seg2rej');
        data_art = ft_selectdata(cfg, data_art);
        art = ft_selectdata(cfg, art);
        
        
        % =====================================================================
        %   2   COMPUTE COMPONENT TIME-SERIES 
        % =====================================================================
        % data_art - structure with segmented MEG channels; 
        % art - structure with segmented relevant channel (ecg/veog/ygaze)
        
        fprintf(['\n ---------- 2 Computing components time-series (' artifact_type{n} ')----------\n'])
        
        cd(['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/comp_ICA/' ID(1:3)])      
        load(sprintf('comp_%s.mat',ID))       % ICA weightings

        % Decompose the ECG-locked datasegments into components, using the previously found (un)mixing matrix
        cfg           = [];
        cfg.unmixing  = comp.unmixing;
        cfg.topolabel = comp.topolabel;
        comp_art      = ft_componentanalysis(cfg, data_art); 

        % Append the artifact channel to the data structure
        % NOTE: Relevant channel is the first one!
        comp_art      = ft_appenddata([], art, comp_art);
        
        mat_name = ['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/comp_ICA/' ID(1:3) '/Matrices'];
        if 7==exist(mat_name,'dir')
            cd(mat_name)
        else
            mkdir(mat_name)
            cd(mat_name)
        end
        
        decomp = [artifact_type{n} '_decomp_' ID '.mat'];
        save(decomp,'comp_art','-v7.3')
        
        % Average the components timelocked to the QRS-complex
        cfg           = [];
        cfg.keeptrials = 'yes';
        timelock      = ft_timelockanalysis(cfg, comp_art); % MEG channels + relevant channel included! 
        timelock.trialvar = var(timelock.trial,0,3); % w = 0
        timelock.avg = squeeze(nanmean(timelock.trial)); % Otherwise partly NaNs are included
        
        timelock_art = [artifact_type{n} '_timelock_' ID '.mat'];
        save(timelock_art,'timelock','-v7.3')
        
        
        % =====================================================================
        %   3   COMPUTE FREQUENCY DECOMPOSITION OF ALL COMPONENTS
        % =====================================================================

        fprintf(['\n ---------- 3 Computing frequency decomposition (' artifact_type{n} ')----------\n'])
        
        % Compute a frequency decomposition of all components and relevant channel 
        cfg            = [];
        cfg.method     = 'mtmfft';
        cfg.output     = 'fourier';
        cfg.foilim     = [0 100];
        cfg.taper      = 'hanning';
        cfg.pad        = 'maxperlen';
        freq           = ft_freqanalysis(cfg, comp_art);
        
        freq_info = [artifact_type{n} '_freq_' ID '.mat'];
        save(freq_info,'freq','-v7.3')
        
        % Compute coherence between all components and relevant channel
        cfg            = [];
        cfg.channelcmb = {'all' chanRej{n}};
        cfg.jackknife  = 'no';
        cfg.method     = 'coh';
        fdcomp         = ft_connectivityanalysis(cfg, freq);
        
        fd_comp = [artifact_type{n} '_fdcomp_' ID '.mat'];
        save(fd_comp,'fdcomp','-v7.3')
        
        % =====================================================================
        %   4   PLOT COHERENCE 
        % =====================================================================
        
        plot_ICAcomps(filein, artifact_type{n})
        
    end

end

