function [oriTuning,oriTuningByTrial,calcium_traces,laserstartframes,CondOfInterest,moving_trial] = main_analysis_v11_Ming( ...
    ID,select,Tserieslist_all,plotDiff,ncond,saveplot,savedata,savedir, ...
    oris,tuningWin,blockTrials,CondOfInterest,blockNum,bline,doInferA,trialDur,speeds,varargin)
% input:
%   nTrials: trials per condition
% global neuron
%% run printcells2.m to extract cell numbers for cell traces
%% Load 'Fall.mat' data
% clear all
% Tserieslist_all = {{'10232024'}};powerLabel = 'maxP';
%%% timing parameters set up
% postim = 19.2; % end of previous photo-stim to the beginning of the second initial delay
% postim = 3.9; % end of previous photo-stim to the beginning of the second initial delay
% initdelay = 4;% initial delay
laserdur = .2;% total duration of photo-stim period
% bline = 4; % in seconds before stimulation
% stim_num = 20;
stim_trace = [];
% for iID = 1:length(miceall)
Tserieslist = Tserieslist_all{1};
% oris = 0:45:315;
colorsCond = [[247 130 150]/255;[194 29 36]/255];
photostimframe = [];
trial = [];
delays = 2; % in seconds
visttimdur = 2+delays;
bline = bline-delays;

for t = 1:length(Tserieslist)
    F = [];
    trialVar = varargin{1};
    if ~isempty(varargin{2})
        bs = varargin{2};% 1 - running, 0-still, 3-all
    end
    datadir = [savedir '/suite2p/Fall.mat'];
    load(datadir,'F','Fneu','iscell')
    % freq = 0.033345982*3; % imaging frequency, s per frame RB 03/06
    freq = 0.033345982; % RB 03/06
    % freq = 0.033474236;
    % ISI = postim+initdelay+stimdur; %interstimulus interval,  seconds
    %% index out cells from 'F', 'Fneu' and 'spks' data
    % Define index condition
    F_1=double(F(find(iscell(:,1)==1),:)); % indexes out data based on above conditions
    Fneu_1=double(Fneu(find(iscell(:,1)==1),:));
    % computing calcium traces
    calcium_traces = F_1-0.7*Fneu_1;
    
    
%     spks_1=spks(cells,:);
    % MWHF: extract block of interests
    
    % behavioralState = 

    %% old code to sort into 
    % if (length(F)/length(blockTrials)*4)>26391
    %     weirdExtra = 26;
    % else
    %     weirdExtra = 0;
    % end
    % blocktemp = [0 cumsum( floor((blockTrials*6+blockpad)/freq) + weirdExtra )];
    trialtemp = [0 cumsum(blockTrials)];


    % identify block indices 
    averagewindow = 1.5;
    behavState = zeros(1,trialtemp(end));
    if bs < 3
        speedthresh = .5;
        [movebase,movestim] = MovingTrial(freq,averagewindow,(speeds),blockTrials,4,speedthresh);
        
        behavState(unique([movebase,movestim])) = 1;
        % behavState(movestim) = 1;
    end

    % pad calcium traces if it is not the integer times of blocks
    % calcium_traces = [calcium_traces,zeros(size(calcium_traces,1),blocktemp(end)-size(calcium_traces,2))];
    % F_1 = [F_1,zeros(size(F_1,1),blocktemp(end)-size(F_1,2))];
    for i = 2:length(trialtemp) % 1 by 8 (all blocks)
        % identifying the end and starting frame for each block
        % block_frames{i-1} = [blocktemp(i-1)+1:blocktemp(i)];% index for each block in the whole time series
        block_trial{i-1} = [trialtemp(i-1)+1:trialtemp(i)];% index for each block's trial in the whole trial sequence
        % block_start(i-1) = blocktemp(i-1)+1; % starting frame for each block in the whole time series
        % block_end(i-1) = blocktemp(i);% end frame for each block in the whole time series
    end
    % check if the last frame is larger than the length of time series
    % if blocktemp(end) > size(calcium_traces,2)
    %     block_frames{end} = block_frames{end} - (blocktemp(end)-size(calcium_traces,2));
    %     block_start(end) = block_start(end) - (blocktemp(end)-size(calcium_traces,2));
    %     block_end(end) = block_end(end) - (blocktemp(end)-size(calcium_traces,2));
    % end
    
    % sorting and dividing into blocks - then focusing only on blocks 
    icondindex = [];
    for i = 1:ncond
        icondindex = [icondindex,find(blockNum==i)];
    end
    % sorting into AABB or AABBCCDD order
    % block_frames = block_frames(icondindex);
    block_trial = block_trial(icondindex);
    % block_start = block_start(icondindex);
    % block_end = block_end(icondindex);
    block_Num = blockNum(icondindex);
    block_blockTrials = blockTrials(icondindex);
    % dividing calcium traces into blocks
    block_cal=[];block_F_1 = [];photostimframe=[];
    for i = 1:ncond
        % block_cal{i} = calcium_traces(:,cell2mat(block_frames((2*i-1):2*i)));% AB or ABCD blocks
        % block_F_1{i} = F_1(:,cell2mat(block_frames((2*i-1):2*i))); % AB or ABCD block
        block_ori{i} = trialVar(cell2mat(block_trial((2*i-1):2*i)));
        cond_trial{i} = cell2mat(block_trial((2*i-1):2*i));
        moving_trial{i} = behavState(cell2mat(block_trial((2*i-1):2*i)));
        block_NumTrials{i} = sum(block_blockTrials((2*i-1):2*i));
        % if ~doInferA
        %     photostimframe{i} = findStimFrame(block_F_1{i});% use B as anchor point
        % else
        %     if i > 1
        %         photostimframe{i} = findStimFrame(block_F_1{i});% use B as anchor point
        %         if i == ncond
        %             photostimframe{1} = photostimframe{2};% use B as anchor point
        %         end
        %     end
        % end
    end
   
    % Calculate calcium changes & select out specific neurons
    % if blockTrials(1) == 32; frameoffset = 1;else frameoffset = 0;end
    % Identify which neurons from 2P software are labelled in suite2p from printcells2.m image
    % neurons = reshape(calcium_traces(select,:),length(select),[],nTrials,ncond);
    % neurons_all = reshape(calcium_traces(select,1:(floor(size(calcium_traces,2)/length(trialVar))*length(trialVar))),length(select),[],length(trialVar));
   
    % extract number of frames for sitmulaiton window
    [~,xindx] = findStimFrame(F_1); % RB 02/20
    LaserDur = xindx(1); % RB 02/20

    % next cut the whole sequences into different trial
    [trial,laserstartframes] = getTrials(calcium_traces,F_1,cond_trial,bline,freq,trialDur,'time'); % RB added 01/22
    % then resort the trials by conditions
    
    % find index for trials for different conditions
    trial = getConditions(trial,block_ori,block_NumTrials);
%%
    % reorder the trials based on the experimental order
    % next find the grating trials only and feed to the original analysis
    % pipepline

    % extract condition of interest for further analysis
    neurons_ori = [];
    for i = 1:length(CondOfInterest)
        % combined gray and ori trial, do not separately plot them
        % find grating trials
        neurons_ori(:,:,:,CondOfInterest(i)) = trial{CondOfInterest(i)}.trial(select,:,:);
        oriVar{CondOfInterest(i)} = block_ori{CondOfInterest(i)};
        %%%%%%%%%%% for grating trials

        % % find grating trials
        % neurons_ori(:,:,:,CondOfInterest(i)) = trial{CondOfInterest(i)}.trial(select,:,trial{CondOfInterest(i)}.index.ori);
        % oriVar{CondOfInterest(i)} = block_ori{CondOfInterest(i)}(block_ori{CondOfInterest(i)}<1000);

        % %%%%%%%%%% for gray trials
        % neurons_gray(:,:,:,CondOfInterest(i))= trial{CondOfInterest(i)}.trial(select,:,trial{CondOfInterest(i)}.index.gray);
        % grayVar{CondOfInterest(i)} = block_ori{CondOfInterest(i)}(block_ori{CondOfInterest(i)}=1000);

        trial_time{CondOfInterest(i)} = (1:(size(neurons_ori(:,:,:,CondOfInterest(i)),2))) * freq - (trialDur - visttimdur);
        % here the trial time is aligned to visual stim onset, which is
        % consistent across all conditions. bline here only means onset of
        % laser

    end
    oriTuning = [];GrayTuning = [];

    if CondOfInterest(1) == CondOfInterest(2)
        % pre post
        preTemp = neurons_ori(:,:,1:size(neurons_ori,3)/2,CondOfInterest(1));
        postTemp = neurons_ori(:,:,size(neurons_ori,3)/2+1:size(neurons_ori,3),CondOfInterest(1));
        oriPreTemp = oriVar{CondOfInterest(1)}(1:size(neurons_ori,3)/2);
        oriPostTemp = oriVar{CondOfInterest(1)}(size(neurons_ori,3)/2+1:size(neurons_ori,3));
        movingPreTemp = moving_trial{CondOfInterest(1)}(1:size(neurons_ori,3)/2);
        movingPostTemp = moving_trial{CondOfInterest(1)}(size(neurons_ori,3)/2+1:size(neurons_ori,3));
        moving_trial{1} = movingPreTemp;
        moving_trial{2} = movingPostTemp;
        
        neurons_ori = [];
        neurons_ori(:,:,:,1) = preTemp;
        neurons_ori(:,:,:,2) = postTemp;
        trial_time_temp{1} = trial_time{CondOfInterest(1)};
        trial_time_temp{2} = trial_time{CondOfInterest(1)};
        trial_time = trial_time_temp;
        oriVar{1} = oriPreTemp;
        oriVar{2} = oriPostTemp;
        bline = repmat(bline(CondOfInterest(1)),1,2);
        
        

        CondOfInterest = [1,2];% 1-pre 2-post
        
    % else % RB added 03/07
    %     % Pre and during (1 session each)
    %     preTemp = neurons_ori(:,:,1:size(neurons_ori,3)/2,CondOfInterest(1));
    %     durTemp = neurons_ori(:,:,1:size(neurons_ori,3)/2,CondOfInterest(2));
    %     oriPreTemp = oriVar{CondOfInterest(1)}(1:size(neurons_ori,3)/2);
    %     oriDurTemp = oriVar{CondOfInterest(2)}(1:size(neurons_ori,3)/2);
    %     movingPreTemp = moving_trial{CondOfInterest(1)}(1:size(neurons_ori,3)/2);
    %     movingDurTemp = moving_trial{CondOfInterest(2)}(1:size(neurons_ori,3)/2);
    %     moving_trial{1} = movingPreTemp;
    %     moving_trial{2} = movingDurTemp;
    % 
    %     neurons_ori = [];
    %     neurons_ori(:,:,:,1) = preTemp;
    %     neurons_ori(:,:,:,2) = durTemp;
    %     trial_time_temp{1} = trial_time{CondOfInterest(1)};
    %     trial_time_temp{2} = trial_time{CondOfInterest(1)};
    %     trial_time = trial_time_temp;
    %     oriVar{1} = oriPreTemp;
    %     oriVar{2} = oriDurTemp;
    %     bline = repmat(bline(CondOfInterest(1)),1,2);
    % 
    %     CondOfInterest = [1,2];% 1-pre 2-during

        


    end


    % [GrayTuning,GrayTuningByTrial] = OriPlot(trial_time,neurons_gray,select,1000,freq,grayVar,bline,visttimdur,tuningWin,saveplot,savedata,Tserieslist,t,ID,ncond);

    
    %% plot the neuron traces
    %%%%%%%%%%%%% raw calcium traces
    calciumPlot(select,calcium_traces,freq,saveplot,savedir,Tserieslist,t,ID,20)
    
    %% plot by orientations 
    [oriTuning,oriTuningByTrial] = OriPlot(trial_time,neurons_ori,select,oris,freq,oriVar,bline,visttimdur,trialDur,tuningWin,saveplot,savedata,Tserieslist,t,ID,CondOfInterest,LaserDur,bs,moving_trial);

    %% plot difference
    if plotDiff & length(CondOfInterest)>1
        figure('Name','SingleNeuron - Difference');
        % set(gcf,'Position',[200 200 2*300/4*length(oris) 1.5*1400/18*length(select)]);
        set(gcf,'Units','Normalized','Position',[0 0 0.6 .6])
        t = tiledlayout(length(select),length(oris), 'TileSpacing', 'tight', 'Padding', 'tight');
        % offset = 1;
        contrastCond = 2;
        OriDiffPlot(trial_time,neurons_ori,select,oris,freq,oriVar,bline,visttimdur,saveplot,contrastCond,trialDur)
        
    end

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% utility functions below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calciumPlot(select,calcium_traces,freq,saveplot,savedir,Tserieslist,t,ID,figNum)
figure(figNum)
set(gcf,'Name',['Calcium traces'])
hold on;
% neurons_plot = zeros(numel(select),length((neurons(1,:)-min(neurons(1,:)))/max(neurons(1,:))*1.1-1));
for i=1:numel(select)
     neurons_plot(i,:)= (calcium_traces(select(i),:)-min(calcium_traces(select(i),:)))/max(calcium_traces(select(i),:))*1.1-i;
     legendname{i} = ['Neuron ',num2str(select(i))];
     colors=rand(i,3);
 end
time_axis = (0:(size(calcium_traces,2)-1))*freq;
colororder(colors);
plot(time_axis,neurons_plot');
% xline(4:8:1024); % RB - vis stim onset
xlabel('Time (s)')
title('Calcium traces')
legend(legendname,'location','bestoutside','FontSize',12);
set(gcf,'Position',[200 200 1200 600]);
%% Save figure
if saveplot
    % saveas(gcf,['Calcium traces for visual stim laser100_' Tserieslist{t} '_Set' num2str(neuset) '.png']);
    saveas(gcf,[savedir, '/Raw_traces_' Tserieslist{t} '_' ID '.png']);
end

function [oriTuning,oriTuningByTrial] = OriPlot(trial_time,neurons,select,oris,freq,oriVar,bline,visttimdur,trialDur,tuningWin,saveplot,savedata,Tserieslist,t,ID,CondOfInterest,LaserDur,bs,moving_trial)
oriTuningByTrial = [];
trial_clrs = [.7,.7,.7,1;.5,.5,.5,.8;.4,.4,.4,.8;.3,.3,.3,.8];
mean_clrs = [[255 187 197]/255,.7;[202 26 33]/255,.7;[202 26 33]/255,.5;[202 26 33]/255,.3];

% light pink 1st condition, brown for second % RB 03/11
% mean_clrs = [[255 187 197]/255,.7; 205/255, 92/255, 0/255, 0.7]; 

figure(100);hold on %%% Ming 05/2025
whichtoplot = [1 2]; % flag for individual vs. mean traces
% set(gcf,'Name',['SingleNeuron',num2str(figNum)]);
set(gcf,'Name',['SingleNeuron - Cond 1&2']);
t2 = tiledlayout(length(select),length(oris), 'TileSpacing', 'tight', 'Padding', 'tight');
% if trial.ntrials(1) == 32; offset = 1;else offset = 1;end
offset=0;
plotted=0;
for iv = whichtoplot
    for icond = CondOfInterest(end:-1:1) 
        
        % Ming - setting the start of frames for computing tuning curve,
        % this will automatically adjust for different delays of laser so
        % that the start is right after laser, in next round adjust the
        % very short delay to be the length of laser stimulation 
        % tuningStart = (floor(max([bline(icond),2])/freq)+1+4)+ floor(offset/freq);
        tuningStart = (floor(max([bline(icond),2])/freq)+1+4)+ floor(offset/freq); 
        baselineEnd = (floor( min([bline(icond),trialDur-visttimdur]) /freq));

        switch bs
            case 1 % running
                trial_select = moving_trial{icond} == 1;
            case 0 % still
                trial_select = moving_trial{icond} == 0;
            case 3 % all
                trial_select = ones(1,length(moving_trial{icond}));
        end

        % update 02/26/25 Ming
        % LaserDur = 
        tunWin_indx = (floor( (trialDur - visttimdur) /freq)):((floor( (trialDur - visttimdur) /freq)) + floor(tuningWin/freq));% 3 = tuning window
        % xline(trial_time{icond}(vstim_indx(1)),'g');
        lstim_indx = (floor( (bline(icond)) /freq)) : (floor( (bline(icond)) /freq) + LaserDur);
        [C,IA,IB] = intersect(tunWin_indx,lstim_indx);
        tunWin_indx(IA) = [];

        % open a figure
        % figNum = 2;% change this to icond if you want to plot separately

        % Ming 05/2025
        numNeurons = size(neurons,1);numConds = length(CondOfInterest);
        % pre-allocate
        p   = nan(numNeurons, numel(oris)-1, numConds);
        h   = nan(numNeurons, numel(oris)-1, numConds);
        for ii=1:length(select)
            figure(100);% Ming 05/2025
            % subplot(ceil(length(select)/1),1,ii); % change depending on size of 'select'
            for jj = 1:length(oris)
                nexttile(jj+(ii-1)*length(oris));hold on
                if plotted < inf
                    % plotted < length(select)*length(oris)
                    % visualize baseline region
                    xregion(trial_time{icond}(1),trial_time{icond}( baselineEnd ),'FaceColor',[.1 .1 .1],'FaceAlpha',.1)
                    % visualize regions used for computing tuning curves
                    % xregion(trial_time{icond}((floor( (trialDur-visttimdur)/freq ))+offset),trial_time{icond}(floor( trialDur/freq - (visttimdur-tuningWin)/freq)+offset),'FaceColor',[0 1 0],'FaceAlpha',.1)
                    
                    % xregion(trial_time{icond}( tuningStart ),trial_time{icond}( floor((trialDur)/freq - (visttimdur-tuningWin) * 30)+ floor(offset/freq) ),'FaceColor',[0 1 0],'FaceAlpha',.1)

                    % xregion(trial_time{icond}( tuningStart ),trial_time{icond}( tuningStart + floor(tuningWin/freq) ),'FaceColor',[0 1 0],'FaceAlpha',.1)

                    % update 2/27/25 ming - segmented plot to ignore the
                    % laser during stimulation
                    xregion([trial_time{icond}( tunWin_indx(1) );trial_time{icond}( tunWin_indx(end) )],'FaceColor',[0 1 0],'FaceAlpha',.1)
                    xregion([trial_time{icond}( lstim_indx(1) );trial_time{icond}( lstim_indx(end) )],'FaceColor',[1 1 1],'FaceAlpha',.8)
                    
                    plotted = plotted+1;
                    % if icond == 2
                    %     xregion(trial_time{icond}(floor((bline(icond)+stimdur)/freq)),trial_time{icond}(floor( (bline(icond)+2*stimdur)/freq )),'FaceColor',[0 1 0])
                    % end
                end
                if iv == 1
                    % indivdiual traces
                    plot(trial_time{icond}, squeeze(neurons(ii,:,oriVar{icond}==oris(jj) & trial_select,icond) - ...
                        repmat(mean(neurons(ii,1:( baselineEnd ),oriVar{icond}==oris(jj) & trial_select,icond),2),1,size(neurons,2),1,1)),'color',trial_clrs(icond,:),'LineWidth',.5)
                        % [0.8 0.8 0.8]-.2)
                elseif iv == 2
                    % mean response 
                    plot(trial_time{icond}, mean(neurons(ii,:,oriVar{icond}==oris(jj)& trial_select,icond) - ...
                        repmat(mean(neurons(ii,1:baselineEnd,oriVar{icond}==oris(jj)& trial_select,icond),2),1,size(neurons,2),1,1),3),'Color',mean_clrs(icond,:),'LineWidth',2)
                end
                xlim([-(trialDur-visttimdur)-.2 visttimdur+.2])
                ylim([min(squeeze(neurons(ii,:,:,:) - ...
                    repmat(mean(neurons(ii,1:baselineEnd,:,:),2),1,size(neurons,2),1,1)),'','all'),...
                    max(squeeze(neurons(ii,:,:,:) - ...
                    repmat(mean(neurons(ii,1:baselineEnd,:,:),2),1,size(neurons,2),1,1)),'','all')])
                %ylim ([-25 inf]) %change depending on y axis and inh or exc neurons
        %         set(gca,'xtick',[])
                % set(gca,'xtick',[-1,0:10:30])
                xlabel('Time (s)');
                ylabel('Calcium changes');
        %                 title(legendname(1,ii))
                % axis off
                xline(trial_time{icond}(floor( (trialDur - visttimdur) /freq)))% visual stim onset
                xline(trial_time{icond}(floor( (bline(icond)) /freq)),'r') % laser onset
                
                yline(0)
                % find the response to the current orientation#
                meanTraces = mean(neurons(ii,:,oriVar{icond}==oris(jj)& trial_select,icond) - ...
                    repmat(mean(neurons(ii,1:baselineEnd,oriVar{icond}==oris(jj)& trial_select,icond),2),1,size(neurons,2),1,1),3);
                
                Traces = neurons(ii,:,oriVar{icond}==oris(jj)& trial_select,icond) - ...
                    repmat(mean(neurons(ii,1:baselineEnd,oriVar{icond}==oris(jj)& trial_select,icond),2),1,size(neurons,2),1,1);
                
                % pad traces with nan for missing trials
                theoretical_trials = size(neurons,3)/length(oris);
                if size(Traces,3)<theoretical_trials
                    Traces(:,:,size(Traces,3)+1:theoretical_trials) = nan(size(Traces,1),size(Traces,2),theoretical_trials - size(Traces,3) );
                end
                % oriTuning(ii,jj,icond) = mean( meanTraces( (floor(bline(icond)/freq)+1+4+ floor(offset/freq)) : (floor((trialDur)/freq - (visttimdur-tuningWin) * 30)+ floor(offset/freq)) ) );
                % oriTuningByTrial(ii,jj,:,icond) = mean( Traces(:, (floor(max([2,bline(icond)])/freq)+1+4 + floor(offset/freq)):(floor((trialDur)/freq - (visttimdur-tuningWin) * 30)+ floor(offset/freq)),1:8 ),2 );
                
                % oriTuning(ii,jj,icond) = mean( meanTraces( tuningStart : tuningStart + floor(tuningWin/freq) ) );
                % oriTuningByTrial(ii,jj,:,icond) = mean( Traces(:, tuningStart : tuningStart + floor(tuningWin/freq) ,1:8 ),2 );

                % update 02/26/25 Ming
                oriTuning(ii,jj,icond) = mean( meanTraces( tunWin_indx ) );
                oriTuningByTrial(ii,jj,:,icond) = mean( Traces(:, tunWin_indx ,: ),2 );
                
                if ii == 1
                    if oris(jj)<1000
                        title(sprintf('Orientation: %d',oris(jj)))
                    else
                        title('Gray Stim')
                    end
                end
            end
            %%% Ming 05/2025: add ranksum test
            % baseline distribution
            if iv==1 
                neurons_baseSub = neurons(:,:,:,:) - ...
                    repmat(mean(neurons(:,1:baselineEnd,:,:),2),1,size(neurons,2),1,1);

                masks = arrayfun(@(o) (oriVar{icond}==o) & trial_select, oris(1:8), ...
                     'UniformOutput', false);
                % extract the two time‐windows once
                baseIdx = 1:baselineEnd;
                % stimIdx = tunWin_indx(1):size(neurons,2);
                stimIdx = tunWin_indx;
      
                % arrayfun over jjj = 1:numel(oris)
                % [pRow, hRow] = arrayfun(@(j) ...
                %     ranksum( ...
                %         reshape(neurons(ii, baseIdx,   masks{j}, icond),[],1), ...
                %         reshape(neurons_baseSub(ii, stimIdx,    masks{j}, icond),[],1), ...
                %         'tail','left' ...
                %     ), 1:numel(oris(1:8)));
                [pRow, hRow] = arrayfun(@(j) ...
                    ranksum( ...
                        reshape(neurons(ii, baseIdx,   :, icond),[],1), ...
                        reshape(neurons_baseSub(ii, stimIdx,    masks{j}, icond),[],1), ...
                        'tail','left' ...
                    ), 1:numel(oris(1:8)));
        
                % store
                p(ii, :, icond) = pRow; 
                h(ii, :, icond) = hRow;
                if length(pRow == min(pRow))>1
                    nopref = oris(pRow == min(pRow));
                    oriPref(ii,icond) = -1;
                    
                else
                    oriPref(ii,icond) = oris(pRow == min(pRow));
                end
                isplot = 1;
                if isplot && icond==1
                    nbins=20;
                    figure(102)
                    set(gcf,'Units','normalized','Position',[0 1-.2*length(select) 1 .2*length(select)])
                    for jjj = 1:8
                        subplot(length(select),9,jjj+(ii-1)*9);hold on
                        % [N1,edges1] = histcounts(reshape(neurons(ii, baseIdx,   masks{jjj}, icond),[],1),nbins,'Normalization', 'probability');
                        [N1,edges1] = histcounts(reshape(neurons(ii, baseIdx,   :, icond),[],1),nbins,'Normalization', 'probability');
                        stairs(edges1,[N1 0],'Color',trial_clrs(icond,:),'LineWidth',3);hold on
                        
                        [N2,edges2] = histcounts(reshape(neurons_baseSub(ii, stimIdx,   masks{jjj}, icond),[],1),nbins,'Normalization', 'probability');
                        stairs(edges2,[N2 0],'Color',mean_clrs(icond,:),'LineWidth',3);hold on
                        currMax = ylim;
                        ylim([0,1.1*max([currMax(2),N1,N2])]);
                        
                        % histogram(base_dist(:),nbins,'Normalization','probability','FaceColor',[.6 .6 .6]*(1/icond));hold on
                        % histogram(stim_dist(:),nbins,'FaceColor',[1 0 0]*(1/icond),'Normalization','probability');
                        text(30,.8*max([N1,N2]),sprintf('p: %.2f',p(ii, jjj, icond)*9))
                        if jjj==8
                            subplot(length(select),9,jjj+1 +(ii-1)*9);hold on
                            plot(oriTuning(ii,1:8,icond),'Color',mean_clrs(icond,:),'LineWidth',3)
                        end
                    end
    
                end
                respNeu = sum(p(:,:,icond)*8<0.05,2)>0;
                if ii == length(select)
                    display(' ')
                    display(sprintf('for Condition %d',icond))
                    display('P values compared to baseline are: ')
                    display(p(:,:,icond)*9)
                    display(' ')
                    display('Significance: ')
                    display(p(:,:,icond)*8<0.05)
                    % display(['For Condition ', num2str(icond) ...
                    %     ', visual responsive neurons are: ', num2str(select(sum(p(:,:,icond)*9<0.05,2)>0))])

                    display(sprintf(['For Condition %d, visual responsive neurons are: ',repmat('%d(%d°) ',1,length(respNeu))],icond,reshape([select(respNeu)',oriPref(respNeu,icond)]',1,[])))
                    display(' ')
                end
            end
        end
        % set(gcf,'Position',[200 200 2*300/4*length(oris) 1.5*1400/18*length(select)]);
        set(figure(100),'Units','normalized','Position',[0 + (length(oris)>1)*.1 0 .8/8*length(oris) .8])
        
    end
end


%% Save figure and workspace
if saveplot
    saveas(gcf,[savedir '/Mean_traces_' Tserieslist{t} '_' ID '.png']);
end
% close gcf
if savedata
    save(['./results_rhys/Calcium_traces/processed_data_' Tserieslist{t} '_' ID '.mat'])
end
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% function to find the anchor frame for each photostimulation
%%%%%%%%%% points
function [startframe,xindx] = findStimFrame(F)
[~,neu_anchor] = max(std(double(F),[],2));% find the best trace to identify photostim points
stimFrames = find((abs(double(F(neu_anchor,:))-min(double(F(neu_anchor,:)))))<1);
% stimFrames = find(double(F(neu_anchor,:))<.1);
[~,xindx] = find((stimFrames(2:end)-stimFrames(1:end-1))>10);
startframe = [stimFrames(1),stimFrames(xindx+1)];% starting frame of each photostimulation point

% startframe = unique(startframe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% function to find the start and end of each trial based on 
%%%%%%%%%% the anchor point of the photostimulation 
function [trial,laserstartframes] = getTrials(calcium_traces,F_1,cond_trial,bline,freq,trialDur,method) % RB added 01/22
if strcmp(method,'shutter')
    [laserstartframes,~] = findStimFrame(F_1);
    for icond = 1:length(cond_trial)
        photostimframe{icond} = laserstartframes(cond_trial{icond});
        for i = 1:length(photostimframe{icond} ) % looping through trials
            trial{icond}.start(i) = photostimframe{icond}(i)-floor((bline(icond))/freq); % frames before photostimulation
            trial{icond}.end(i) = photostimframe{icond}(i) + ceil((trialDur-bline(icond))/freq);% frames after photostimulation
            trial{icond}.trial(:,:,i) = calcium_traces(:,trial{icond}.start(i):trial{icond}.end(i));
        end
    end
elseif strcmp(method,'time') % reshape to extract trials w/o laser 
    [laserstartframes,~] = findStimFrame(F_1);
    nblocks = length(cond_trial)*2;
    ntrials = length(cond_trial{1})/2;
    % try using reshape to cut into trials
    % y1 = calcium_traces(10,:)';
    % blocks = reshape(y1,[],nblocks);

    x = reshape(calcium_traces,size(calcium_traces,1),[],nblocks);% n neurons by time by block
    % x = reshape(F_1,size(F_1,1),[],nblocks);% n neurons by time by block

    % plot(x(10,:,2)-calcium_traces(10,size(x,2)+1:2*size(x,2)));hold on
    % plot(calcium_traces(1,1:size(x,2)),'--')

    theoreticalLength = (ntrials*6+4)/freq;
    shrinkageFactor = length(x)/theoreticalLength;
    laserpoints = floor(((1:((6/freq)):length(x))+4/freq)*shrinkageFactor);
    trialStart = floor(((1:((6/freq)):length(x))+2/freq)*shrinkageFactor);
    trials_all = [];
    for ib = 1:nblocks
        % trials_all = [trials_all;reshape(x(:,:,ib),size(x,1),[],36)];
        for itr = 1:ntrials
            indxtemp = trialStart(itr):trialStart(itr+1)-1;
            if length(indxtemp)<180
                trials_all(:,:,(ib-1)*ntrials+itr) = x(:,trialStart(itr):trialStart(itr+1),ib);
            else
                trials_all(:,:,(ib-1)*ntrials+itr) = x(:,trialStart(itr):trialStart(itr+1)-1,ib);
            end
        end
    end
    % plot(squeeze(trials_all(2,:,cond_trial{2})))

    if 1-1e-3<=shrinkageFactor &  shrinkageFactor<= 1+1e-3
        warning(['Actual test is the same as expected: ',num2str(length(x)) '/' num2str(round(theoreticalLength)) 'frames'])
    elseif shrinkageFactor < 1-1e-3
        warning(['Actual test is shorted than expected: ',num2str(length(x)) '/' num2str(round(theoreticalLength)),'frames'])
    elseif shrinkageFactor > 1+1e-3
        warning(['Actual test is longer than expected: ',num2str(length(x)) '/' num2str(round(theoreticalLength)),'frames'])
    end

    for icond = 1:length(cond_trial)
        trial{icond}.trial = trials_all(:,:,cond_trial{icond});
    end
    % plot(squeeze(trial{3}.trial(1,:,:)))

    % plot(mean(squeeze(trial{3}.trial(1,:,:)),2))

    isplot =  0;
    if isplot
        figure
        set(gcf,'units','normalized','position',[0 .8 1 .2])
        plot(F_1(1,1:size(x,2)))
        % plot(x(:,:,1)','r');hold on
        xline(laserpoints)
    end
        
        



    % initialDelay = 2;% check this for your experiment, this should be fixed at 2s in normal cases 
    % StartFrames = initialDelay/freq + (1:round(trialDur/freq-1):length(calcium_traces)-4/freq);
    % padindex = zeros(1,length(StartFrames));padindex(37:36:length(padindex)) = round(2/freq);
    % StartFrames = StartFrames + padindex;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% function to find the trials index (grating vs. gray) for each
%%%%%%%%%% condition
function trial = getConditions(trial,block_ori,block_NumTrials)
for icond = 1:length(trial) % looping through conditions
    % number of total trials for each condition
    trial{icond}.ntrials = block_NumTrials{icond};
    % index for ori trials
    trial{icond}.index.ori = find(block_ori{icond}<1000);
    % index for gray screen trials
    trial{icond}.index.gray = find(block_ori{icond}==1000);
    
end

%%%
function OriDiffPlot(trial_time,neurons_ori,select,oris,freq,oriVar,bline,visttimdur,saveplot,whichCond,trialDur)
offset = 0;
for ii=1:length(select)
    % subplot(ceil(length(select)/1),1,ii); % change depending on size of 'select'

    for jj = 1:length(oris)
        meanpre = [];meanpost=[];meandiff=[];
        nexttile(jj+(ii-1)*length(oris));hold on

        % neurons: ncells by time points by trials by condition
        diffs = (neurons_ori(ii,:,oriVar{whichCond}==oris(jj),whichCond) - ...
            repmat(mean(neurons_ori(ii,1:(round( min([bline(whichCond),trialDur-visttimdur])/freq )),oriVar{whichCond}==oris(jj),2),whichCond),1,size(neurons_ori,2),1,1)) - ...
        mean(neurons_ori(ii,:,oriVar{1}==oris(jj),1) - ...
            repmat(mean(neurons_ori(ii,1:(round(min([bline(1),trialDur-visttimdur])/freq)),oriVar{1}==oris(jj),1),2),1,size(neurons_ori,2),1,1),3);

        diffall = (neurons_ori(ii,:,:,whichCond) - ...
            repmat(mean(neurons_ori(ii,1:(round(min([bline(whichCond),trialDur-visttimdur])/freq)),:,whichCond),2),1,size(neurons_ori,2),1,1)) - ...
        mean(neurons_ori(ii,:,:,1) - ...
            repmat(mean(neurons_ori(ii,1:(round(min([bline(1),trialDur-visttimdur])/freq)),:,1),2),1,size(neurons_ori,2),1,1),3);

        plot(trial_time{whichCond}, squeeze(diffs),'Color',[0.8 0.8 0.8]-.2)

        % plot vis
        % meanpre = mean(neurons(ii,:,oriVar{1}==oris(jj),1) - ...
        %     repmat(mean(neurons(ii,1:(round(bline/freq)),oriVar{1}==oris(jj),1),2),1,size(neurons,2),1,1),3);
        % plot(trial_time,meanpre,'Color',colorsCond(1,:))

        % plot vis+opto
        % meanpost = mean(neurons(ii,:,oriVar{2}==oris(jj),2) - ...
        %     repmat(mean(neurons(ii,1:(round(bline/freq)),oriVar{2}==oris(jj),2),2),1,size(neurons,2),1,1),3);
        % plot(trial_time,meanpost,'Color',colorsCond(2,:))

        % plot difference
        meandiff = mean(neurons_ori(ii,:,oriVar{whichCond}==oris(jj),whichCond) - ...
            repmat(mean(neurons_ori(ii,1:(round(min([bline(whichCond),trialDur-visttimdur])/freq)),oriVar{whichCond}==oris(jj),whichCond),2),1,size(neurons_ori,2),1,1),3) - ...
            mean(neurons_ori(ii,:,oriVar{1}==oris(jj),1) - ...
            repmat(mean(neurons_ori(ii,1:(round(min([bline(1),trialDur-visttimdur])/freq)),oriVar{1}==oris(jj),1),2),1,size(neurons_ori,2),1,1),3);
        plot(trial_time{whichCond}, meandiff,'Color',[0 0 1 .5])
        xlim([-(trialDur-visttimdur)-.1 visttimdur+.1])
        % ylim([ ...
        %     min((neurons(ii,:,:,2) - repmat(mean(neurons(ii,1:(round(bline/freq)),:,2),2),1,size(neurons,2),1,1)) - ...
        %         (neurons(ii,:,:,1) - repmat(mean(neurons(ii,1:(round(bline/freq)),:,1),2),1,size(neurons,2),1,1)),'','all'),...
        %     max((neurons(ii,:,:,2) - repmat(mean(neurons(ii,1:(round(bline/freq)),:,2),2),1,size(neurons,2),1,1)) - ...
        %         (neurons(ii,:,:,1) - repmat(mean(neurons(ii,1:(round(bline/freq)),:,1),2),1,size(neurons,2),1,1)),'','all')])
        ylim([min(diffall(:)), max(diffall(:))])



        % ylim([min([meanpre;meanpost;meandiff],'','all'),...
        %     max([meanpre;meanpost;meandiff],'','all')])
        
        %ylim ([-25 inf]) %change depending on y axis and inh or exc neurons
%         set(gca,'xtick',[])
        % set(gca,'xtick',[-1,0:10:30])
        xlabel('Time (s)');
        ylabel('Calcium changes');
%                 title(legendname(1,ii))
        axis off
        % xline(trial_time(round((bline)/freq)+offset))
        xline(trial_time{whichCond}(floor( (trialDur - visttimdur) /freq)+offset))% visual stim onset
        xline(trial_time{whichCond}(floor( (bline(whichCond)) /freq)+offset)) % laser onset
        yline(0)
        if ii == 1
            title(sprintf('Orientation: %d',oris(jj)))
            if jj==1
                % legend({'Vis','Opto+vis','Diff'})
            end
        end
    end
end

function [movebase,movestim] = MovingTrial(samplefreq,averagewindow,speeds,blockTrials,padframenum,speedthre)
allmean = [];
for i = 1:length(speeds.speed_data)
    blockframenum {i} = length(speeds.speed_data{i}); % = 6599 each block
end
% speedthre = 1;
orisequence = cell2mat(speeds.speed_data');
averagecore = round(averagewindow/samplefreq);% =45 frames
% absorisequence = abs(orisequence);
orisequence = abs(orisequence);
allmean = movmean(orisequence, averagecore, 'Endpoints','shrink')/samplefreq;
idx = round(linspace(1,length(allmean)+1,length(blockTrials)+1));
for i = 1:length(blockTrials)
    blockglidmean{i} = allmean(idx(i):idx(i+1)-1);
end
for i = 1:length(blockTrials) % cut off the 4s padding from each block 
    blockglidmean{i} = blockglidmean{i}(1:blockframenum{i}-padframenum);
end
sequence = []; % splicing blocks with padding cutted
for i = 1:length(blockglidmean)
    sequence = [sequence, blockglidmean{i}];
end
Trialnum = sum(blockTrials); % = 4*36 =144 trials
framepertrial = length(sequence)/Trialnum; % = 180 frames/trial
for i = 1:Trialnum % cut the whole sequences & store data of each frame in a cell element
    trialsequence{i} = sequence((i-1)*framepertrial + 1:i*framepertrial);
end
stimmovmean = cellfun(@(x) x(121:end), trialsequence, 'UniformOutput', false);
basemovmean = cellfun(@(x) x(1:120), trialsequence, 'UniformOutput', false);% anonymous function indexing each row vector in cell
% find moving trials and trials with speed threshold
nonzerostim = find(cellfun(@(x) any(x(:) ~= 0), stimmovmean));
nonzerobase = find(cellfun(@(x) any(x(:) ~= 0), basemovmean));
movestim = find(cellfun(@(x) any(x > speedthre), stimmovmean));
movebase = find(cellfun(@(x) any(x > speedthre), basemovmean));

