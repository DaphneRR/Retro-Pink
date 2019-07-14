%% Analysis -- retropink data ---------------------------------------------
% Written by Daphné Rimsky-Robert -> drr@tuta.io

clear
close all
clc

datapath = 'add repository path here';
cd(datapath)
count = 1;

% data was pre-processed in matlab: trials with RTs under 150ms were 
% removed, trials where timings were not respected also (hardware failure)
% data from all included participants were pooled into a single csv file
% for easy access and use in R and Matlab.

% read pre-processed data from csv starting below the header
data = csvread('RP4.csv',1,0);
sj_ID = unique(data(:,1));



for iSj = 1:length(sj_ID)
    
    % reformat data for matlab and readability
    subjDat = data(:,1) == sj_ID(iSj);
    d.rt = data(subjDat,2);
    d.accuracy = data(subjDat,3);
    d.congr = data(subjDat,4);
    d.soa = data(subjDat,5);
    d.targLoc = data(subjDat,6);
    d.cueLoc = data(subjDat,7);
    d.response = ~isnan(d.rt);
    
    % make relevant logical index vectors
    Hit = d.response == 1 & ~isnan(d.targLoc);
    FA = d.response == 1 & isnan(d.targLoc);
    goTrials = ~isnan(d.targLoc);
    noCue = isnan(d.cueLoc);
    
    %% calculate pFA for cue present trials
    trials = ~noCue & isnan(d.targLoc);
    pFA = double(FA(trials));
    pFA_mean = nnz(pFA)/length(pFA);
    
    if pFA_mean == 0
        pFA_mean = 1/(2*length(pFA));
    elseif pFA_mean == 1
        pFA_mean = 1 - 1/(2*length(pFA));
    end
    pFA_all(count) = pFA_mean;
    
    %% get data for cue-absent trials

    trials = noCue;
    pHit = double(Hit(trials));
    pHit_mean = nnz(pHit)/length(pHit);
    
    if pHit_mean == 0
        pHit_mean = 1/(2*length(pHit));
    elseif pHit_mean == 1
        pHit_mean = 1 - 1/(2*length(pHit));
    end
    
    trials = noCue;
    pFA = double(FA(trials));
    pFA_mean = nnz(pFA)/length(pFA);
    
    if pFA_mean == 0
        pFA_mean = 1/(2*length(pFA));
    elseif pFA_mean == 1
        pFA_mean = 1 - 1/(2*length(pFA));
    end
    
    pHit_noCue(count) = pHit_mean;
    pFA_noCue(count) = pFA_mean;
    
    [dprime_noCue(count), beta_noCue(count), c_noCue(count)] = dprime(pHit_mean, pFA_mean);
        
    trials = noCue & goTrials == 1 & d.response == 1;
    RT_noCue(count) = nanmedian(d.rt(trials));
    
    trials = noCue & goTrials == 0 & d.response == 1;
    RT_noCue_FA(count) = nanmedian(d.rt(trials));
    
    trials =  ~noCue & goTrials == 0 & d.response == 1;
    RT_FA(count) = nanmedian(d.rt(trials));

    %% get data per SOA
    for i = 1:4
        % get percent correct ---------------------------------------------
        all_Soa = [-550, -100, 200, 500];
        mySoa = all_Soa(i);

        % RT analysis and histograms
        trials = d.soa == mySoa & goTrials == 1 & d.response == 1;
        RT_all(i,count) = nanmedian(d.rt(trials));
        
        trials = d.soa == mySoa  & d.congr == 1 & goTrials == 1 & d.response == 1;
        RT_C(i,count) = nanmedian(d.rt(trials));
        
        trials = d.soa == mySoa  & d.congr == 0 & goTrials == 1 & d.response == 1;
        RT_IC(i,count) = nanmedian(d.rt(trials));

        % calculate dprime and accuracy for all trials
        trials = d.soa == mySoa & goTrials == 1;
        pHit = double(Hit(trials));
        pHit_mean = nnz(pHit)/length(pHit);
        if pHit_mean == 0
            pHit_mean = 1/(2*length(pHit));
        elseif pHit_mean == 1
            pHit_mean = 1 - 1/(2*length(pHit));
        end
        pHit_all(i,1,count) = pHit_mean;
        [dprime_all(i,count), beta_all(i, count), c_all(i,count)] = dprime(pHit_all(i,1,count), pFA_all(count));
        
                
        % calculate dprime and accuracy for d.congr == 1
        trials = d.soa == mySoa & d.congr == 1 & goTrials == 1;
        pHit = double(Hit(trials));
        pHit_mean = nnz(pHit)/length(pHit);
        if pHit_mean == 0
            pHit_mean = 1/(2*length(pHit));
        elseif pHit_mean == 1
            pHit_mean = 1 - 1/(2*length(pHit));
        end
        pHit_all(i,2,count) = pHit_mean;
        [dprime_C(i,count), beta_C(i, count), c_C(i,count)] = dprime(pHit_all(i,2,count), pFA_all(count));
        
        
        % calculate dprime and accuracy for d.congr == 0
        trials = d.soa == mySoa & d.congr == 0 & goTrials == 1;
        pHit = double(Hit(trials));
        pHit_mean = nnz(pHit)/length(pHit);
        if pHit_mean == 0
            pHit_mean = 1/(2*length(pHit));
        elseif pHit_mean == 1
            pHit_mean = 1 - 1/(2*length(pHit));
        end
        pHit_all(i,3,count) = pHit_mean;
        [dprime_IC(i,count), beta_IC(i, count), c_IC(i,count)] = dprime(pHit_all(i,3,count), pFA_all(count));
    end
    
    count = count + 1;
    
    
    
    
end
%% means for plotting
dprime_C_all = nanmean(dprime_C,2);
dprime_IC_all = nanmean(dprime_IC,2);
dprime_noCue_all = nanmean(dprime_noCue,2);

RT_C_all = nanmean(RT_C,2);
RT_IC_all = nanmean(RT_IC,2);
RT_noCue_all = nanmean(RT_noCue,2);
dprime_nocue_sem = std(dprime_noCue)/sqrt(length(dprime_noCue));
RT_nocue_sem = std(RT_noCue)/sqrt(length(RT_noCue));
%% compute SED -------------------------------------------------------------
for j = 1:4
    RT_SED(j) = nansed(RT_C(j,:),RT_IC(j,:));
    dprime_SED(j) = nansed(dprime_C(j,:), dprime_IC(j,:));
end
%% statistical tests


detectionAnova = [dprime_C' dprime_IC'];
[efs,F,cdfs,p,eps,dfs,b,y2,sig]=repanova(detectionAnova,[2,4]);


RTAnova = [RT_C' RT_IC'];
[efs,F,cdfs,p,eps,dfs,b,y2,sig]=repanova(RTAnova,[2,4]);

% confidence intervals
nboot = 10000;
bootfun = @(x,y) nanmean(x - y);

ci_soa = bootci(nboot,{bootfun, dprime_all', repmat(dprime_noCue,4,1,1)'}, 'alpha', 0.05/4)
ci_congr = bootci(nboot,{bootfun, dprime_C',dprime_IC'}, 'alpha', 0.05/4)
ci = bootci(nboot,{bootfun, RT_all', repmat(RT_noCue,4,1,1)'}, 'alpha', 0.05/4)
ci = bootci(nboot,{bootfun, RT_C',RT_IC'}, 'alpha', 0.05/4)
%% figures 

figure
hold on
title('Detection (d'')')
% plot target appearance
AREA = area([0 53], [4 4], 'FaceColor', [.4 .4 .4], 'EdgeColor', 'none');
% plot BL
plot([-600 -150 150 450], [dprime_noCue_all dprime_noCue_all dprime_noCue_all dprime_noCue_all], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([-450 -450], [dprime_noCue_all - dprime_nocue_sem, dprime_noCue_all + dprime_nocue_sem])
% plot full lines
errorbar([-600 -150], dprime_C_all(1:2), dprime_SED(1:2), '.-', 'Color', [.4 .6 .9], 'MarkerSize', 45, 'LineWidth', 1.5);
errorbar([-600 -150], dprime_IC_all(1:2), dprime_SED(1:2), '.-', 'Color', [.9 0 0], 'MarkerSize', 45, 'LineWidth', 1.5);
errorbar([150 450], dprime_C_all(3:4), dprime_SED(3:4), '.-', 'Color', [.4 .6 .9], 'MarkerSize', 45, 'LineWidth', 1.5);
errorbar([150 450], dprime_IC_all(3:4), dprime_SED(3:4), '.-', 'Color', [.9 0 0], 'MarkerSize', 45, 'LineWidth', 1.5);
% plot dotted lines
plot([-150 150], dprime_C_all(2:3), '.--', 'Color', [.4 .6 .9], 'LineWidth', 3)
plot([-150 150], dprime_IC_all(2:3), '.--', 'Color', [.9 0 0], 'LineWidth', 3)
% adjust linewidth without changing error bars
plot([-600 -150], dprime_C_all(1:2),'.-', 'Color', [.4 .6 .9], 'LineWidth', 3);
plot([-600 -150], dprime_IC_all(1:2), '.-', 'Color', [.9 0 0], 'LineWidth', 3)
plot([150 450], dprime_C_all(3:4),'.-', 'Color', [.4 .6 .9], 'LineWidth', 3)
plot([150 450], dprime_IC_all(3:4),'.-', 'Color', [.9 0 0], 'LineWidth', 3)

ylim([2.1 2.9])
xlim([-610 460])
set(gca, 'XTick', [-600 -150 150 450],'FontSize', 25, 'YTick', [2 2.2 2.4 2.6 2.8])
AREA.Face.ColorType = 'truecoloralpha';
drawnow; pause(0.05);  % This needs to be done for transparency to work
AREA.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3
set(gcf,'Position', [752 284 963 814])



% RT for hits
figure
hold on
AREA = area([0 53], [700 700], 'FaceColor', [.4 .4 .4], 'EdgeColor', 'none');
% plot BL
plot([-600 -150 150 450], [RT_noCue_all RT_noCue_all RT_noCue_all RT_noCue_all], '--', 'Color', [.6 .6 .6], 'LineWidth', 2)
plot([-450 -450], [RT_noCue_all - RT_nocue_sem, RT_noCue_all + RT_nocue_sem])
% plot full lines
errorbar([-600 -150], RT_C_all(1:2), RT_SED(1:2), '.-', 'Color', [.4 .6 .9], 'MarkerSize', 45, 'LineWidth', 1.5)
errorbar([-600 -150], RT_IC_all(1:2), RT_SED(1:2), '.-', 'Color', [.9 0 0], 'MarkerSize', 45, 'LineWidth', 1.5);
errorbar([150 450], RT_C_all(3:4), RT_SED(3:4), '.-', 'Color', [.4 .6 .9], 'MarkerSize', 45, 'LineWidth', 1.5)
errorbar([150 450], RT_IC_all(3:4), RT_SED(3:4), '.-', 'Color', [.9 0 0], 'MarkerSize', 45, 'LineWidth', 1.5);
% plot dotted lines
plot([-150 150], RT_C_all(2:3), '.--', 'Color', [.4 .6 .9], 'LineWidth', 3)
plot([-150 150], RT_IC_all(2:3), '.--', 'Color', [.9 0 0], 'LineWidth', 3)
% adjust linewidth without changing error bars
plot([-600 -150], RT_C_all(1:2),'.-', 'Color', [.4 .6 .9], 'LineWidth', 3);
plot([-600 -150], RT_IC_all(1:2), '.-', 'Color', [.9 0 0], 'LineWidth', 3)
plot([150 450], RT_C_all(3:4),'.-', 'Color', [.4 .6 .9], 'LineWidth', 3)
plot([150 450], RT_IC_all(3:4),'.-', 'Color', [.9 0 0], 'LineWidth', 3)
ylim([460 580])
xlim([-610 460])
set(gca, 'XTick', [-600 -150 150 450],'FontSize', 25);%, 'YTick', [480 520 560 600 640] )
AREA.Face.ColorType = 'truecoloralpha';
drawnow; pause(0.05);  % This needs to be done for transparency to work
AREA.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3
set(gcf,'Position', [752 284 963 814])


