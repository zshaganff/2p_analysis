%%% Analysis Pipeline for the optogenetic stimulation project
% MWHF & RB
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
% MWHF: Assuming you have done the preprocessing using suite2p
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
clear all
close all
addpath(genpath('..'))
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
% MWHF: below are some parameters that need to be manual changed or filled
% out before you run the analysis
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
% mouseID = '3592';

% mouse3579_11222024
% mouse3592_12042024
% mouse3579_12062024 
% mouse3592_12102024
% mouse3580_12102024 % working one
% mouse4330_12122024

% mouse3580_12032024  
% mouse3579_11292024
% mouse4338_01082025
% mouse4301_01092025
% mouse4330_01092025
% mouse4301_01132025
% dbTest
% mouse4338_01152025_combined
% mouse4338_01152025_plane0
% mouse4338_01152025_plane1
% mouse4338_01152025_plane2
% mouse4330_01152025_combined
% mouse4330_01152025_plane0
% mouse4330_01152025_plane1
% mouse4330_01152025_plane2
% MingTest
% test2
% mouse3580_01172025
% mouse4330_01202025_combined
% mouse4330_01202025_plane0
% mouse4330_01202025_plane1
% mouse4330_01202025_plane2
% mouse4338_01202025_combined
% mouse4338_01202025_plane0
% mouse4338_01202025_plane1
% mouse4338_01202025_plane2
% mouse4332_01232025_plane0
% mouse4332_01232025_plane1
% mouse4332_01232025_plane2
% mouse4334_01242025 %%%%
% mouse4330_01272025_plane0
% mouse4330_01272025_plane1
% mouse4330_01272025_plane2
% mouse4301_01282025
% mouse4332_01302025_plane0
% mouse4332_01302025_plane1
% mouse4332_01302025_plane2
% mouse4334_01302025
% mouse4338_02032025_plane0
% mouse4338_02032025_plane1
% mouse4338_02032025_plane2
% mouse4301_02042025 %%%%% Ming used to debug
 % mouse4334_02042025
% mouse4330_02032025_plane0
% mouse4330_02032025_plane1
% mouse4330_02032025_plane2
% mouse4338_02102025
% mouse4338_02172025
% mouse4330_02172025
% mouse4338_02172025_PrePost
% mouse4330_02172025_PrePost
% mouse4332_03042025
% mouse4330_02172025_PrePost
% mouse4332_03042025
% mouse4330_03042025
% mouse4330_03122025
% mouse4332_03122025
% mouse4388_03202025
% mouse4330_03262025
% mouse4330_03122025
% mouse4332_03262025
% mouse4332_03122025
% mouse4388_03262025
% mouse3592_04022025
% % mouse4332_04032025
% mouse3592_04152025
% mouse4401_03252025
% mouse4478_05022025
 % mouse4478_05092025
 mouse4484_07012025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. run printcells to see ROI, markpoints and raw tracesZ:\SunLab_sharedfiles\Project-2photon\Data\4301\02042025
    % printcells3(datadir,showROINum,saveplot,singlescan)
    printcells4_RB(datadir,showROINum,ROInum,saveplot,ch1scan,singlescan) % RB - overlays stim points onto ch1.png in 2nd subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2. This is the formal analysis to see neuronals traces and try to separate any stimulation effect
close all
oriVar = oriMerge(oriName); % MWHF: to extract grating orientation for each block/condition (nconditions by ntrials Matrix)
%%% Main analysis script
tuningWin = 3; % in seconds
plotDiff = 0;
% load speed data
behavioralState = 3;% 1 - running, 0 - still, 3 - all;
if behavioralState < 3
    % speeds = load(['../Data/',mouseID,'/',Tserieslist_all{1}{1},'/', 'treadmill/treadmill_data.mat']);
    % speeds = load(['../Data/',mouseID,'/',Tserieslist_all{1}{1},'/',plane,'/','treadmill/treadmill_data.mat']);% RB 03/07
    speeddata = [datadir '/treadmill/treadmill_data.mat'];
    speeds = load(speeddata);
else
    speeds = [];
end

% main analysis
[~,oriTuningByTrial,calcium_traces,laserstartframes,CondOfInterest,moving_trial] = main_analysis_v11_Ming(mouseID,ROInum, ...
    Tserieslist_all,plotDiff,ncond,0,savedata, ...
    datadir,angles,tuningWin,blockTrials, ...
    CondOfInterest,blockNum,bline,doInferA,trialDur,speeds,oriVar,behavioralState);

%% 
close all
% TuningAnalysis(oriTuningByTrial,singlescan,datadir,ROInum,angles,'min')
% TuningAnalysis_v3_Ming(oriTuningByTrial,singlescan,datadir,ROInum,angles,'min',CondOfInterest)
[oriTuning] = TuningAnalysis_v3_Ming_copy(oriTuningByTrial,singlescan,datadir,ROInum,angles,'min',CondOfInterest);
% saveas(gcf,'/Volumes/144.82.131.60/SunLab_sharedfolder/SunLab_sharedfiles/Project-2photon/Analysis/figurex','png')

