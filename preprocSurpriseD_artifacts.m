function [] = preprocSurpriseD_artifacts(idx_file)
% PREPROCESSING surpriseD meg data (continuous approach)
    % Identify and save times of strong head movements, blinks, saccades, squid
    % jumps and muscle artifacts - relative to block onset!

    % Block bounds +/- 10 seconds - be careful with cfg.latency because of
    % the offset!!!
    
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
    [~,files] = xlsread('Info_filewise');
    info_EL_blocks = xlsread('Info_filewise');
    filein = files{idx_file};
    
    %idx_file = find(strcmp(files, filein));
    
    % One of these files contains 4 blocks of meg data
    ID = [filein(1:5) filein(end-5:end-3)]; % Subject ID + Session number + file number
    if regexp(ID, 'URG_S*') % With this subject the session number was missing when registered
        ID = ['URG-1' filein(end-5:end-3)];
    end

    addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
    addpath '/mnt/homes/home024/jschipp/Surprise_Drug/meg_preprocessing'
    addpath '/mnt/homes/home024/jschipp/Surprise_Drug/meg_data/'
    addpath '/mnt/homes/home024/chernandez/meg_data/'
    ft_defaults


    %% Loop through blocks
    %   Filein contains normally 4 blocks of data. Loop trough each block separately
    %   Load data of relevant channels, downsampling to 400Hz, filter out 50Hz
    %   and harmonics (line noise), transform eyelink data from voltage into
    %   pixel values, identify times of blink, saccade, jump and muscle
    %   artifacts and store these matrices
    
    firstBl = 1;
    
    for block = firstBl:info_EL_blocks(idx_file,3)

        % folder with the data
        if idx_file < 230 % Older datasets (Caro)
             cd /mnt/homes/home024/pmurphy/meg_data/surpriseD/
        else % Newer datasets (Julia)
             cd /mnt/homes/home024/jschipp/Surprise_Drug/meg_data/  
        end
      
        fprintf('\n\n ---------------- \n PROCESSING FILE %s, BLOCK #%d...\n ---------------- \n\n', ID, block)

        % =====================================================================
        %   1   DEFINE BLOCK AND LOAD RELEVANT CHANNELS
        % =====================================================================

        fprintf('\n ---------- 1 Defining block and relevant channels ----------\n')

        % Specify cfg
        cfgin = [];
        cfgin.dataset = filein;
        cfgin.block = block;
        cfgin.ID = ID;
        cfgin.trialdef.prestim = 10; % Add 10 seconds before and after the block
        cfgin.trialdef.poststim = 10;
        cfgin.trialfun = 'trialfun_surpriseD_continuous';
        cfgin = ft_definetrial(cfgin);
        cfgin.continuous = 'yes'; % read in data as continuous
        cfgin.channel = {'meg','EEG001','EEG002','EEG003','EEG057','EEG058','HLC*','UADC*'};
        %cfgin.channel = 'meg';
        data = ft_preprocessing(cfgin);

        % Initial power spectrum
        cfgfreq = [];
        cfgfreq.method = 'mtmfft';
        cfgfreq.output = 'pow';
        cfgfreq.taper = 'hanning';
        cfgfreq.channel = 'MEG';
        cfgfreq.foi = 1:130;
        cfgfreq.keeptrials = 'no';
        freq = ft_freqanalysis(cfgfreq, data);

        cnt = 1;
        figure('vis', 'off'), clf
        subplot(2,4,cnt); cnt = cnt + 1;
        loglog(freq.freq, freq.powspctrm, 'linewidth', 0.1); hold on;
        loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
        axis tight; axis square; box off;
        set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'FontSize', 8);
        title('Power spectrum raw');

        clear freq
        
        % =====================================================================
        %   INTERMEDIATE STEP! - PART OF THE ARTEFACT DETECTION
        %   needs to be done before downsampling & filtering
        % =====================================================================
        
        %   a) SACCADES BASED ON EYELINK DATA
        fprintf('\n ---------- a) Identifying saccades -----------\n')
        
        % Transform eyelink data into pixels
        ranges = [5 -5]; % Voltage range
        screen_x = [0 1920];
        screen_y= [0 1080];

        % Find eyelink channels
        % -002 = Xgaze, -004 = Ygaze, -003 = pupil
        ch_mapping(1) = find(strcmp(data.label,'UADC002'));
        ch_mapping(2) = find(strcmp(data.label,'UADC004'));
        ch_mapping(3) = find(strcmp(data.label,'UADC003'));
        
        [data.trial{1}(ch_mapping(1),:), data.trial{1}(ch_mapping(2),:), ~] = ...
            eye_voltage2gaze(data.trial{1}(ch_mapping,:), ranges, screen_x, screen_y, ch_mapping);

        ppd = estimate_pixels_per_degree();
        xcenter = screen_x(2)/2;
        ycenter = screen_y(2)/2;

        x = data.trial{1}(ch_mapping(1),:);
        y = data.trial{1}(ch_mapping(2),:);

        hz = data.fsample; % sampling frequency after resampling (400 Hz)
        threshold = 30; % taken from Niklas' script
        acc_thresh = 2000; % taken from Niklas' script
        amp_thresh = 1.5; % Saccadic amplitude threshold - meant to exclude all but most extreme microsaccades

        % Detect saccades with velocity acceleration approach
        % Returns an Nx2 matrix with start and end samples of each saccade artifact
        artifact_saccade = check_saccade_vel_acc(x, y, hz, threshold, acc_thresh, amp_thresh, ppd);

        % =====================================================================
        %   2    DOWN SAMPLE DATA TO 400 HZ
        % =====================================================================

        fprintf('\n ---------- 2 Resampling block -----------\n')

        cfgres.resample = 'yes';
        cfgres.fsample = 1200;
        cfgres.resamplefs = 400;
        cfg.detrend = 'no';
        data = ft_resampledata(cfgres, data);
        
        % Downsample artifact matrix saccade
        artifact_saccade = round(artifact_saccade/cfgres.fsample*cfgres.resamplefs);

        % =====================================================================
        %   3    REMOVE LINE NOISE
        % =====================================================================

        fprintf('\n ---------- 3 Filtering out line noise (50 Hz) and its harmonics -----------\n')

        cfg             = [];
        cfg.bsfilter    = 'yes';
        cfg.bsfreq      = [49 51; 99 101; 149 151];
        data            = ft_preprocessing(cfg, data);

        % Plot power spectrum after removing line noise
        freq = ft_freqanalysis(cfgfreq, data);

        subplot(2,4,cnt); cnt = cnt + 1;
        loglog(freq.freq, freq.powspctrm, 'linewidth', 0.1); hold on;
        loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
        axis tight; axis square; box off;
        set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'xticklabel', []);
        title('No line noise');

        clear freq

        % =====================================================================
        %   4   DETECT HEAD MOVEMENTS
        % =====================================================================
        % Detect strong head movements with help of head localization data
        fprintf('\n ---------- 4 Identifying changes in head position -----------\n')

        % Compute head location
        cc_rel = computeHeadRotationRest(data);

        % Plot head motion
        subplot(2,4,cnt); cnt = cnt + 1;
        plot(data.time{1},cc_rel); set(gca, 'tickdir', 'out', 'FontSize', 8);
        axis tight; box off;
        title('Head motion')

        % Find outliers
        [~, idx] = deleteoutliers(cc_rel);
        [t,~]    = ind2sub(size(cc_rel),idx);

        % Only take those where the deviation is more than 6 mm
        t = t(any(abs(cc_rel(t, :)) > 6, 2));

        clear cc_rel

        % create Mx2 matrix of sample vector "t" with end start & end sample of motion
        artifact_headM = [];
        if length(t) > 1
            for i = 1:length(t)-1
                if t(i)+1 == t(i+1)
                    artifact_headM(i,1)=t(i);
                    artifact_headM(i,2)=t(i+1);
                else
                    artifact_headM(i,1)=t(i);
                    artifact_headM(i,2)=t(i)+1;
                end
            end
        end

        % Merge consecutive head motions into 1 if they're X ms together
        coalesce = 1;
        if size(artifact_headM,1)>1
            cmsmp = artifact_headM(1,:);
            for b = 1:size(artifact_headM,1)-1
                if artifact_headM(b+1,1) - cmsmp(end,2) < coalesce*data.fsample
                    cmsmp(end,2) = artifact_headM(b+1,2);
                else
                    cmsmp(end+1,:) = artifact_headM(b+1,:);
                end
            end
            artifact_headM = cmsmp; clear cmsmp x y
        end

        % =====================================================================
        %   5   IDENTIFY ARTIFACTS
        % =====================================================================

        fprintf('\n ---------- 5 Identifying artifacts -----------\n')

        % _____________________________________________________________________
        %   b)  BLINK ARTIFACTS
        fprintf('\n ---------- b) Identifying blinks -----------\n')

        % If eyelink data is good (defined in "eyelink_ok') use eyelink data,
        % otherwise use veog channel for blink detection
        if info_EL_blocks(idx_file,5) == 1
            % --- GETTING BLINKS BASED ON Y-GAZE  (EYELINK DATA)
            cfg = [];
            cfg.continuous = 'yes';
            cfg.artfctdef.zvalue.channel = 'UADC004'; % Ygaze
            cfg.artfctdef.zvalue.trlpadding = 0;
            cfg.artfctdef.zvalue.fltpadding = 0;
            cfg.artfctdef.zvalue.artpadding = 0.05;

            % Algorithmic parameters
            cfg.artfctdef.zvalue.bpfilter = 'no';
            cfg.artfctdef.zvalue.hilbert = 'no';

            % Set cutoff
            ygaze_cutoff = 200;  % cutoff in absolute pixels to be used for blink detection (anything below is blink)
            cfg.artfctdef.zvalue.cutoff = -(ygaze_cutoff-mean(data.trial{1}(ch_mapping(2),:)))/std(data.trial{1}(ch_mapping(2),:));
            cfg.artfctdef.zvalue.interactive = 'no';
            data.trial{1}(ch_mapping(2),:) = data.trial{1}(ch_mapping(2),:).*-1;
            [~, artifact_Ygaze] = ft_artifact_zvalue(cfg, data);
            data.trial{1}(ch_mapping(2),:) = data.trial{1}(ch_mapping(2),:).*-1;
        else
            % --- GETTING BLINKS BASED ON VEOG
            cfg.artfctdef.zvalue.channel = 'EEG057'; % VEOG
            cfg.artfctdef.zvalue.trlpadding = 0;
            cfg.artfctdef.zvalue.fltpadding = 0;
            cfg.artfctdef.zvalue.artpadding =  0.05;

            % Algorithmic parameters
            cfg.artfctdef.zvalue.bpfilter = 'yes';
            cfg.artfctdef.zvalue.bpfilttype = 'but';
            cfg.artfctdef.zvalue.bpfreq = [1 15];
            cfg.artfctdef.zvalue.bpfiltord = 4;
            cfg.artfctdef.zvalue.hilbert = 'yes';

            % Set cutoff
            cfg.artfctdef.zvalue.cutoff = 2;
            cfg.artfctdef.zvalue.interactive = 'no';
            [~, artifact_Ygaze] = ft_artifact_zvalue(cfg, data);
        end

        cfg = [];
        cfg.artfctdef.reject = 'partial';
        cfg.artfctdef.eog.artifact = artifact_Ygaze;

        subplot(2,4,cnt); cnt = cnt + 1;
        plot(data.time{1},data.trial{1}(ch_mapping(2),:))
        axis tight; axis square; box off;
        title('Y-gaze')
        
         %   a)  Addition: Saccade artifact
         
        % Remove blink artifacts from saccade artifact matrix
        for n_sacc = 1:size(artifact_saccade,1)
            for n_blink = 1:size(artifact_Ygaze,1)
                 if (artifact_saccade(n_sacc,1) > artifact_Ygaze(n_blink,1)) && (artifact_saccade(n_sacc,2) < artifact_Ygaze(n_blink,2))
                     artifact_saccade(n_sacc,:) = 0;
                 end
            end
        end

        artifact_saccade(artifact_saccade == 0) = [];
        artifact_saccade = reshape (artifact_saccade, [], 2); 
        
        % Merge consecutive saccade artifacts into 1 if they're X ms together
        coalesce = 0.2;
        if size(artifact_saccade,1)>1
            cmsmp = artifact_saccade(1,:);
            for b = 1:size(artifact_saccade,1)-1
                if artifact_saccade(b+1,1) - cmsmp(end,2) < coalesce*data.fsample
                    cmsmp(end,2) = artifact_saccade(b+1,2);
                else
                    cmsmp(end+1,:) = artifact_saccade(b+1,:);
                end
            end
            artifact_saccade = cmsmp; clear cmsmp x y
        end
        
        % _____________________________________________________________________
        %   c)  SQUID JUMPS
        % Compute the power spectrum of all trials and fit a line on the loglog-
        % transformed power spectrum. Jumps cause broad range increase in the power
        % spectrum so trials containing jumps can be selected by detecting outliers
        % in the intercepts of the fitted lines (using Grubb?s test for outliers).
        fprintf('\n ---------- c) Identifying squid jumps -----------\n')

        artifact_Jump = findSquidJumps(data,ID);

        % Merge consecutive jump artifacts into 1 if they're X ms together
        coalesce = 1;
        if size(artifact_Jump,1)>1
            cmsmp = artifact_Jump(1,:);
            for b = 1:size(artifact_Jump,1)-1
                if artifact_Jump(b+1,1) - cmsmp(end,2) < coalesce*data.fsample
                    cmsmp(end,2) = artifact_Jump(b+1,2);
                else
                    cmsmp(end+1,:) = artifact_Jump(b+1,:);
                end
            end
            artifact_Jump = cmsmp; clear cmsmp
        end

        if isempty(artifact_Jump)
            title(sprintf('No jumps'));
        end

        % Calculate power spectra only for clean epochs
%         if isempty(artifact_Jump)
%             freq = ft_freqanalysis(cfgfreq, data);
% 
%             subplot(2,4,cnt); cnt = cnt + 1;
%             loglog(freq.freq, freq.powspctrm, 'linewidth', 0.1); hold on;
%             loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
%             clear freq
%         else
%             powC={};
%             cfg=[]; cfg.channel = {'MEG'};
%             for b = 1:size(artifact_Jump,1)
%                 if b==1 && artifact_Jump(1,1)-data.fsample>1  % if first data segment is clean
%                     cfg.latency = [min(data.time{1}) (artifact_Jump(1,1)-1)./data.fsample];
%                     cdata = ft_selectdata(cfg,data);
%                     cfreq = ft_freqanalysis(cfgfreq, cdata);
%                     powC{end+1,1} = cfreq.powspctrm; powC{end,2} = cfreq.freq;
%                 end
%                 if b~=size(artifact_Jump,1) % following clean period
%                     cfg.latency = ([artifact_Jump(b,2)+1 artifact_Jump(b+1,1)-1])./data.fsample;
%                     cdata = ft_selectdata(cfg,data);
%                     cfreq = ft_freqanalysis(cfgfreq, cdata);
%                     powC{end+1,1} = cfreq.powspctrm; powC{end,2} = cfreq.freq;
%                 elseif b==size(artifact_Jump,1) && artifact_Jump(end,2)<length(data.time{1})-data.fsample % if >1s clean @ end
%                     cfg.latency = [(artifact_Jump(b,2)+1)./data.fsample max(data.time{1})];
%                     cdata = ft_selectdata(cfg,data);
%                     cfreq = ft_freqanalysis(cfgfreq, cdata);
%                     powC{end+1,1} = cfreq.powspctrm; powC{end,2} = cfreq.freq;
%                 end
%             end
%             clear cdata cfreq
% 
%             % plot power spectra separately for each clean period
%             subplot(2,4,cnt);  cnt = cnt + 1; hold on
%             for b = 1:size(powC,1), loglog(powC{b,2}, powC{b,1}, 'linewidth', 0.1); end  % per channel
%             for b = 1:size(powC,1), loglog(powC{b,2}, mean(powC{b,1},1), 'k', 'linewidth', 1); end  % channel-averaged
%             clear powC
%         end
        
%         axis tight; axis square; box off;
%         set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'FontSize', 8, 'xscale', 'log', 'yscale', 'log');
%         title('No line noise');


        % _____________________________________________________________________
        %   d) MUSCLE ARTIFACTS
        fprintf('\n ---------- c) Identifying muscle artefacts -----------\n')

        cfg                              = [];
        cfg.continuous                   = 'yes'; % data has been epoched

        % channel selection, cutoff and padding
        cfg.artfctdef.zvalue.channel     = {'MEG'}; % make sure there are no NaNs
        cfg.artfctdef.zvalue.trlpadding  = 0;
        cfg.artfctdef.zvalue.fltpadding  = 0;
        cfg.artfctdef.zvalue.artpadding  = 0.1;
        cfg.artfctdef.zvalue.interactive = 'no';

        % Algorithmic parameters
        cfg.artfctdef.zvalue.bpfilter    = 'yes';
        cfg.artfctdef.zvalue.bpfreq      = [110 140];
        cfg.artfctdef.zvalue.bpfiltord   = 9;
        cfg.artfctdef.zvalue.bpfilttype  = 'but';
        cfg.artfctdef.zvalue.hilbert     = 'yes';
        cfg.artfctdef.zvalue.boxcar      = 0.2;
        cfg.artfctdef.zvalue.cutoff      = 20;  % THIS MAY NEED TO BE CHECKED AND ADJUSTED FOR DIFFERENT PARTICIPANTS
        [~, artifact_Muscle]             = ft_artifact_zvalue(cfg, data);

        % Merge consecutive muscle artifacts into 1 if they're Xms together
        coalesce = 1;
        if size(artifact_Muscle,1)>1
            cmsmp = artifact_Muscle(1,:);
            for b = 1:size(artifact_Muscle,1)-1
                if artifact_Muscle(b+1,1) - cmsmp(end,2) < coalesce*data.fsample
                    cmsmp(end,2) = artifact_Muscle(b+1,2);
                else
                    cmsmp(end+1,:) = artifact_Muscle(b+1,:);
                end
            end
            artifact_Muscle = cmsmp; clear cmsmp
        end

        % Calculate proportion of total block time that's marked as artifactual
        if ~isempty(artifact_Muscle)
            perc_muscle = sum(artifact_Muscle(:,2)-artifact_Muscle(:,1)+1)/length(data.time{1}).*100;
        else perc_muscle = 0;
        end
        
        % _____________________________________________________________________
        %   e) CAR/OTHER ARTEFACTS
        fprintf('\n ---------- d) Identifying car artefacts -----------\n')
        
        metalThreshold = 2*10^-11; % THIS MAY NEED TO BE CHECKED/ADJUSTED
        % Deviding data in epochs max and min value per channel, mean
        % across all channels 
        
        epochL = data.fsample*2;
        
        for e = 1:floor(size(data.trial{1,1},2)/epochL)
            if e == 1
                data_epochs(:,:,e)= data.trial{1,1}(:,1:epochL);
            else
                data_epochs(:,:,e)= data.trial{1,1}(:,(e-1)*epochL+1:epochL*e);
            end
        end 
        
        maxEpoch = squeeze(max(data_epochs,[],2)); % Find maximum in each epoch
        minEpoch = squeeze(min(data_epochs,[],2)); % FInd minimum in each epoch
        
        megCH = strfind(data.label, 'M'); % Only consider MEG channels 
        megNum = length(megCH(~cellfun('isempty',megCH)));
  
        diffEpoch = maxEpoch(1:megNum,:)-minEpoch(1:megNum,:); % Calculate difference across epoch (per channel, per epoch)
        
        % avgdiff = mean(diffEpoch);
        maxdiff = max(diffEpoch);
        
%         figure(10)
%         plot(maxdiff), hold on, plot(avgdiff);
%         title(['Finding metal artifacts (max-min in epochs) ' ID ' Block ' num2str(block)])
%         legend('channel with maximal difference','average difference across channels');
%         % Save figure
%         figpath = '/mnt/homes/home024/jschipp/Surprise_Drug/meg_preprocessing/test_figures/';
%         figMetal = [figpath 'Metal_' ID '_Block_' num2str(block)];
%         saveas(gca,figMetal,'fig')
        
        idx_metal = find(maxdiff > metalThreshold);
        bad_epochs = zeros(floor(size(data.trial{1,1},2)/epochL),1)';
        bad_epochs(idx_metal) = 1;
   
        
        % If there's an artefact, insert into metal_artefact matrx
        artifact_metal = [];
        for n = 1:length(idx_metal) % Calculate samples from artefact
            artifact_metal(n,1) = idx_metal(n)*epochL-epochL+1;
            artifact_metal(n,2) = idx_metal(n)*epochL;
        end
        
        % Merge consecutive metal artifacts into 1 if they're Xms together
        coalesce = 1;
        if size(artifact_metal,1)>1
            cmsmp = artifact_metal(1,:);
            for b = 1:size(artifact_metal,1)-1
                if artifact_metal(b+1,1) - cmsmp(end,2) < coalesce*data.fsample
                    cmsmp(end,2) = artifact_metal(b+1,2);
                else
                    cmsmp(end+1,:) = artifact_metal(b+1,:);
                end
            end
            artifact_metal = cmsmp; clear cmsmp
        end
    
        % Calculate channel-averaged power spectra for clean/artifactual epochs
%         if ~isempty(artifact_Muscle)
%             powC={}; powM={};
%             cfg=[]; cfg.channel = {'MEG'}; cfgfreq.foi = 5:150;
%             for b = 1:size(artifact_Muscle,1)
%                 if b==1 && artifact_Muscle(1,1)-data.fsample>1  % if first data segment is clean
%                     cfg.latency = [min(data.time{1}) (artifact_Muscle(1,1)-1)./data.fsample];
%                     cdata = ft_selectdata(cfg,data);
%                     cfreq = ft_freqanalysis(cfgfreq, cdata);
%                     powC{end+1,1} = mean(cfreq.powspctrm,1); powC{end,2} = cfreq.freq;
%                 end
%                 cfg.latency = artifact_Muscle(b,:)./data.fsample; % muscle period
%                 cdata = ft_selectdata(cfg,data);
%                 cfreq = ft_freqanalysis(cfgfreq, cdata);
%                 powM{end+1,1} = mean(cfreq.powspctrm,1); powM{end,2} = cfreq.freq;
%                 if b~=size(artifact_Muscle,1) % following clean period
%                     cfg.latency = ([artifact_Muscle(b,2)+1 artifact_Muscle(b+1,1)-1])./data.fsample;
%                     cdata = ft_selectdata(cfg,data);
%                     cfreq = ft_freqanalysis(cfgfreq, cdata);
%                     powC{end+1,1} = mean(cfreq.powspctrm,1); powC{end,2} = cfreq.freq;
%                 elseif b==size(artifact_Muscle,1) && artifact_Muscle(end,2)<length(data.time{1})-data.fsample % if >1s clean @ end
%                     cfg.latency = [(artifact_Muscle(b,2)+1)./data.fsample max(data.time{1})];
%                     cdata = ft_selectdata(cfg,data);
%                     cfreq = ft_freqanalysis(cfgfreq, cdata);
%                     powC{end+1,1} = mean(cfreq.powspctrm,1); powC{end,2} = cfreq.freq;
%                 end
%             end
%             clear cdata cfreq
% 
%             % Plot power spectra separately for clean and muscle periods
%             subplot(2,4,cnt); hold on
%             for b = 1:size(powC,1), plot(powC{b,2}, powC{b,1},'k','linewidth', 1); end
%             for b = 1:size(powM,1), plot(powM{b,2}, powM{b,1},'r','linewidth', 1); end
%             axis tight; axis square; box off; %ylim(ylims);
%             title(['Muscle, n=',num2str(size(artifact_Muscle,1)),', p=',num2str(round(perc_muscle,1))])
%             set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'yscale','log');
%         end


        %%
        % =====================================================================
        %   6  SAVE ARTIFACT MATRICES AND OVERVIEW FIGURE
        % =====================================================================
        fprintf('\n ---------- 6 Saving artifact matrices and figure -----------\n')

        cd('/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/preICA_artifactMatrices/')
        name = ['/mnt/homes/home024/jschipp/Surprise_Drug/meg_analysis/preICA_artifactMatrices/' ID(1:3)];

        if 7==exist(name,'dir')
            cd(name)
        else
            mkdir(name)
            cd(name)
        end

        % Save the artifacts
        artifacts = ['Artifacts_' ID '_Block_' num2str(block)];
        save(artifacts,'artifact_headM','artifact_Ygaze','artifact_saccade','artifact_Muscle','artifact_Jump','artifact_metal')

        % Save the invisible figure
        figOverview = ['Overview_' ID '_Block_' num2str(block)];
        saveas(gca,figOverview,'png')

        % Cleaning up before next block
        fprintf('\n ---------- Cleaning up -----------\n')
        clear data artifact_headM artifact_Ygaze artifact_saccade artifact_Muscle artifact_Jump artifact_metal

        close all
    end
end

