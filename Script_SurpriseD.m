%% SCRIPT FOR PROCESSING AND ANALYSING MEG DATA
% +++ Steps:
% 1. Preprocessing: Artifact matrices for ICA 
% 2. Clean data (muscle, jumps, saccades) & run ICA in order to identify artifacts (blinks, heart beat)
% 3. Plot 25 components with the highest coherence values
% 4. Identify the components that need to be removed manually
% 5. Prepare the trial based analysis // Output: clean trial data per file
% 6. Concatenated trials across run.

% +++ NOTE:
% Steps 1-5 per file
% Step 4 manually!
% Step 6 per subject
%%

clear all
close all

addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
addpath '/mnt/homes/home024/jschipp/Surprise_Drug/meg_preprocessing'
addpath '/mnt/homes/home024/jschippSurprise_Drug//meg_data/'
ft_defaults

% File to process (includes normally 4 blocks of data)
filein = 'CTG-2_Surprise_20180807_01.ds';
subjSess = filein(1:5);

% =====================================================================
%   1    PREPROCESSING: ARTIFACT MATRICES (Preparation ICA)
% =====================================================================

fprintf('\n +++ 1 PREPROCESSING: ARTIFACT MATRICES (Preparation ICA) +++ \n')

preprocSurpriseD_artifacts(filein);

% =====================================================================
%   2    CLEAN DATA, ALIGN ARTIFACTS, PERFORM ICA
% =====================================================================

fprintf('\n +++ 2 CLEAN DATA, ALIGN ARTIFACTS, PERFORM ICA +++ \n')

ICA_SurpriseD(filein);

% =====================================================================
%   3    PLOT 25 COMPONENTS WITH TOP COHERENCE VALUES
% =====================================================================

fprintf('\n +++ PLOT 25 COMPONENTS WITH TOP COHERENCE VALUES +++ \n')

identify_ICA_comps(filein);



%% !!! Manually !!!

% =====================================================================
%   4    IDENTIFY BLINK & ECG COMPONENTS
% =====================================================================

import_comps2rej();



%%
% =====================================================================
%   5    PREPARE TRIALBASED ANALYSIS - RESULT: CLEAN TRIAL DATA
% =====================================================================

fprintf('\n +++ 5    PREPARE TRIALBASED ANALYSIS +++ \n')

prep_trialAnalysis(filein);

