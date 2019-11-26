function cp_groupdemo
% Demo on synthetic 1D group data.
% 1/ create data for Ns subjects
% 2/ smooth them with classic Gaussian or tissue-weighted Gaussian
% 3/ build explicit masks
% 
% Then show some results...
% 
% Q: should I consider the T-SPOON approach ?
%__________________________________________________________________________
% Copyright (C) 2019 GIGA Institute

% Written by C. Phillips, 2019.
% Cyclotron Research Centre, University of Liege, Belgium

%% Parameters
% Number of subjects
Ns = 20;
r_jitter = 1;
P_GmWmCsf = cell(Ns,1);
pP_GmWmCsf = cell(3,1);
gP_GmWmCsf = cell(Ns,1);
ggP_GmWmCsf = cell(3,1);
twsP_signal = cell(Ns,1);

%% Do the processing
% Deal with the Ns subjects, one at a time.
for ii=Ns:-1:1
    % Create the signal + tissue probs
    [P_signal(ii,:), P_GmWmCsf{ii}] = cp_create_data(r_jitter);
    % Smooth the signals, Gaussian & tissue-weighted
    data_ii = struct('P_signal',P_signal(ii,:),'P_GmWmCsf',P_GmWmCsf{ii});
    [gP_signal(ii,:),gP_GmWmCsf{ii},twsP_signal{ii}] = ...
        cp_smooth_data(data_ii,8);
    % Reorganize 
    %   * smoothed tissue classes (-> expl mask) and 
    %   * tissue-weighted smoothed signals
    for jj=1:3
        ggP_GmWmCsf{jj}(ii,:) = gP_GmWmCsf{ii}(jj,:);
        pP_GmWmCsf{jj}(ii,:) = P_GmWmCsf{ii}(jj,:);
    end
    for jj=1:2
        ttwsP_signal{jj}(ii,:) = twsP_signal{ii}(jj,:);
    end
end
% Create explicit mask
[exMask] = cp_explmask(ggP_GmWmCsf);

% Create clean signal
[cP_signal, cP_GmWmCsf, T_names] = cp_create_data(0,[0 0 0]); %#ok<*ASGLU>

%% Check how the mean smoothed signal matches the original signal
% Measure Root Mean Square Error, overall and over each segment based on 
% explicit mask, for the G-smoothed and TW-smoothed signals, w.r.t. the
% true noise-free signal

% Deal with GM
%-------------
% All GM in explicit mask
l_GMall = find(exMask(1,:));
RMSE_GMall = sqrt([...
    sum((cP_signal(l_GMall) - mean(P_signal(:,l_GMall))).^2) / ...
        numel(l_GMall) ;
    sum((cP_signal(l_GMall) - mean(gP_signal(:,l_GMall))).^2) / ...
        numel(l_GMall) ;
    sum((cP_signal(l_GMall) - mean(ttwsP_signal{1}(:,l_GMall))).^2) / ...
        numel(l_GMall)]);
% By segment
l_end_segm = [0 find(diff(l_GMall)>1) numel(l_GMall)];
Nsegm_GM = numel(l_end_segm)-1; % number of segmented
Nel_GMSegm = diff(l_end_segm); % number of voxels/segm
RMSE_GMsegm = zeros(2,Nsegm_GM);
for ii=1:Nsegm_GM
    l_ii = l_GMall(l_end_segm(ii)+1):l_GMall(l_end_segm(ii+1));
    RMSE_GMsegm(1,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(P_signal(:,l_ii))).^2) / ...
            Nel_GMSegm(ii) );
    RMSE_GMsegm(2,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(gP_signal(:,l_ii))).^2) / ...
            Nel_GMSegm(ii) );
    RMSE_GMsegm(3,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(ttwsP_signal{1}(:,l_ii))).^2) / ...
            Nel_GMSegm(ii) );
end
    
% Deal with WM
%-------------
% All WM in explicit mask
l_WMall = find(exMask(2,:));
RMSE_WMall = sqrt([...
    sum((cP_signal(l_WMall) - mean(P_signal(:,l_WMall))).^2) / ...
        numel(l_WMall) ;
    sum((cP_signal(l_WMall) - mean(gP_signal(:,l_WMall))).^2) / ...
        numel(l_WMall) ;
    sum((cP_signal(l_WMall) - mean(ttwsP_signal{2}(:,l_WMall))).^2) / ...
        numel(l_WMall)]);
% By segment
l_end_segm = [0 find(diff(l_WMall)>1) numel(l_WMall)];
Nsegm_WM =numel(l_end_segm)-1; % number of segmented
Nel_WMSegm = diff(l_end_segm); % number of voxels/segm
RMSE_WMsegm = zeros(2,Nsegm_WM);
for ii=1:Nsegm_WM
    l_ii = l_WMall(l_end_segm(ii)+1):l_WMall(l_end_segm(ii+1));
    RMSE_WMsegm(1,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(P_signal(:,l_ii))).^2) / ...
        Nel_WMSegm(ii) );
    RMSE_WMsegm(2,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(gP_signal(:,l_ii))).^2) / ...
        Nel_WMSegm(ii) );
    RMSE_WMsegm(3,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(ttwsP_signal{2}(:,l_ii))).^2) / ...
        Nel_WMSegm(ii) );
end

% Plot values
% figure,
% plot(Nel_GMSegm,RMSE_GMsegm,'bo',Nel_WMSegm,RMSE_WMsegm,'ro')

% Print out some numbers
fprintf('RMSE over the explicit mask for \n')
fprintf('\tGM signal, no-sm %2.2f, G-sm %2.2f and TW-sm %2.2f\n',RMSE_GMall(1),RMSE_GMall(2),RMSE_GMall(3))
fprintf('\tWM signal, no-sm %2.2f, G-sm %2.2f and TW-sm %2.2f\n',RMSE_WMall(1),RMSE_WMall(2),RMSE_WMall(3))
fprintf('\n')
fprintf('RMSE ratio\n')
fprintf('\tGM, no-sm/TW-m %2.2f and G-sm/TW-m %2.2f\n',RMSE_GMall(1)/RMSE_GMall(3),RMSE_GMall(2)/RMSE_GMall(3))
fprintf('\tWM, no-sm/TW-m %2.2f and G-sm/TW-m %2.2f\n',RMSE_WMall(1)/RMSE_WMall(3),RMSE_WMall(2)/RMSE_GMall(3))

%% Plot things

% Original signals, noisy and mean, + same but G-smoothed 
% Plot things, eiter into a single figure or multiple ones.
% -> use sub-function for each (sub)plot !

% 1/ Plot everything in a single figure
figure,
subplot(2,2,1)
plot_Psignal(P_signal,cP_signal)
subplot(2,2,2)
plot_gPsignal(gP_signal,cP_signal)
subplot(2,2,3)
plot_twsPsignal(ttwsP_signal,exMask,cP_signal)
subplot(4,2,6)
plot_pPGmWmCsf(pP_GmWmCsf)
subplot(4,2,8)
plot_ggPGmWmCsf(ggP_GmWmCsf,exMask)
set(gcf,'Position',[500 150 1600 1200])
saveas(gcf,'TissueW_smoothing_demo.png');

% 2/ Plot in different figures
figure,
plot_Psignal(P_signal,cP_signal)
saveas(gcf,'demo_OriginalSignal.png');

figure,
plot_gPsignal(gP_signal,cP_signal)
saveas(gcf,'demo_GsmoothedSignal.png');

figure,
plot_twsPsignal(ttwsP_signal,exMask,cP_signal)
saveas(gcf,'demo_TWsmoothedSignal.png');

figure,
subplot(2,1,1)
plot_pPGmWmCsf(pP_GmWmCsf)
subplot(2,1,2)
plot_ggPGmWmCsf(ggP_GmWmCsf,exMask)
saveas(gcf,'demo_TissueProb.png');

end
%% PLOTTING FUNCTIONS

function plot_Psignal(P_signal,cP_signal)

plot(P_signal','LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(mean(P_signal),'LineWidth',2,'Color',[.1 .1 .1])
plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
title('Noisy signals, its mean (-), and true signal (--)')

end

function plot_gPsignal(gP_signal,cP_signal)

plot(gP_signal','LineWidth',.3,'Color',[1 .8 .8])
hold on
plot(mean(gP_signal),'LineWidth',2,'Color',[.1 .1 .1])
plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
title('Smoothed noisy signals, its mean (-), and true signal (--)')

end

function plot_twsPsignal(ttwsP_signal,exMask,cP_signal)

Ns = size(ttwsP_signal{1},1);
plot(ttwsP_signal{1}'.*(exMask(1,:)'*ones(1,Ns)),'LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(ttwsP_signal{2}'.*(exMask(2,:)'*ones(1,Ns)),'LineWidth',.3,'Color',[1 .8 .8])
mM_ttws1P_signal = mean(ttwsP_signal{1}).*exMask(1,:); 
mM_ttws1P_signal(mM_ttws1P_signal==0) = NaN;
mM_ttws2P_signal = mean(ttwsP_signal{2}).*exMask(2,:); 
mM_ttws2P_signal(mM_ttws2P_signal==0) = NaN;
plot(mM_ttws1P_signal,'LineWidth',2,'Color',[0 0 1],'LineStyle','-')
plot(mM_ttws2P_signal,'LineWidth',2,'Color',[1 0 0],'LineStyle','-')
% plot(mean(ttwsP_signal{1}).*exMask(1,:),'LineWidth',2,'Color',[.1 .1 .1])
% plot(mean(ttwsP_signal{2}).*exMask(2,:),'LineWidth',2,'Color',[.1 .1 .1])
plot(cP_signal,'LineWidth',2,'Color',[.2 .2 .2],'LineStyle','--')
title('Tissue-w. smoothed noisy signals, its mean (-), and true signal (--)')

end

function plot_pPGmWmCsf(pP_GmWmCsf)

plot(pP_GmWmCsf{1}','LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(mean(pP_GmWmCsf{1}),'LineWidth',2,'Color',[0 0 1])
plot(pP_GmWmCsf{2}','LineWidth',.3,'Color',[1 .8 .8])
plot(mean(pP_GmWmCsf{2}),'LineWidth',2,'Color',[1 0 0])
title('Noisy tissue probabilities and their mean (-)')

end

function plot_ggPGmWmCsf(ggP_GmWmCsf,exMask)

plot(ggP_GmWmCsf{1}','LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(mean(ggP_GmWmCsf{1}),'LineWidth',2,'Color',[0 0 1])
plot(ggP_GmWmCsf{2}','LineWidth',.3,'Color',[1 .8 .8])
plot(mean(ggP_GmWmCsf{2}),'LineWidth',2,'Color',[1 0 0])
plot(exMask(1,:),'LineWidth',2,'Color',[0 0 .5],'LineStyle','--')
plot(exMask(2,:),'LineWidth',2,'Color',[.5 0 0],'LineStyle','--')
title('Smoothed noisy tissue prob, their mean (-), and explicit mask (--)')

end

%% SOME OLD STUFF TO KEEP AT HAND...

% figure,
% 
% subplot(2,2,1)
% plot(P_signal','LineWidth',.3,'Color',[.8 .8 1])
% hold on
% plot(mean(P_signal),'LineWidth',2,'Color',[.1 .1 .1])
% plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
% title('Noisy signals, its mean (-), and true signal (--)')
% 
% subplot(2,2,2)
% plot(gP_signal','LineWidth',.3,'Color',[1 .8 .8])
% hold on
% plot(mean(gP_signal),'LineWidth',2,'Color',[.1 .1 .1])
% plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
% title('Smoothed noisy signals, its mean (-), and true signal (--)')
% 
% subplot(2,2,3)
% hold on
% plot(ttwsP_signal{1}'.*(exMask(1,:)'*ones(1,Ns)),'LineWidth',.3,'Color',[.8 .8 1])
% plot(ttwsP_signal{2}'.*(exMask(2,:)'*ones(1,Ns)),'LineWidth',.3,'Color',[1 .8 .8])
% mM_ttws1P_signal = mean(ttwsP_signal{1}).*exMask(1,:); 
% mM_ttws1P_signal(mM_ttws1P_signal==0) = NaN;
% mM_ttws2P_signal = mean(ttwsP_signal{2}).*exMask(2,:); 
% mM_ttws2P_signal(mM_ttws2P_signal==0) = NaN;
% plot(mM_ttws1P_signal,'LineWidth',2,'Color',[0 0 1],'LineStyle','-')
% plot(mM_ttws2P_signal,'LineWidth',2,'Color',[1 0 0],'LineStyle','-')
% % plot(mean(ttwsP_signal{1}).*exMask(1,:),'LineWidth',2,'Color',[.1 .1 .1])
% % plot(mean(ttwsP_signal{2}).*exMask(2,:),'LineWidth',2,'Color',[.1 .1 .1])
% plot(cP_signal,'LineWidth',2,'Color',[.2 .2 .2],'LineStyle','--')
% title('Tissue-w. smoothed noisy signals, its mean (-), and true signal (--)')
% 
% subplot(4,2,6)
% plot(pP_GmWmCsf{1}','LineWidth',.3,'Color',[.8 .8 1])
% hold on
% plot(mean(pP_GmWmCsf{1}),'LineWidth',2,'Color',[0 0 1])
% plot(pP_GmWmCsf{2}','LineWidth',.3,'Color',[1 .8 .8])
% plot(mean(pP_GmWmCsf{2}),'LineWidth',2,'Color',[1 0 0])
% title('Noisy tissue probabilities and their mean (-)')
% 
% subplot(4,2,8)
% plot(ggP_GmWmCsf{1}','LineWidth',.3,'Color',[.8 .8 1])
% hold on
% plot(mean(ggP_GmWmCsf{1}),'LineWidth',2,'Color',[0 0 1])
% plot(ggP_GmWmCsf{2}','LineWidth',.3,'Color',[1 .8 .8])
% plot(mean(ggP_GmWmCsf{2}),'LineWidth',2,'Color',[1 0 0])
% plot(exMask(1,:),'LineWidth',2,'Color',[.1 .1 .1],'LineStyle','--')
% plot(exMask(2,:),'LineWidth',2,'Color',[.1 .1 .1],'LineStyle','--')
% title('Smoothed noisy tissue prob, their mean (-), and explicit mask (--)')
% 
% set(gcf,'Position',[500 150 1600 1200])

% 
% % Original signals, noisy and mean, + same but G-smoothed 
% figure,
% subplot(3,1,1)
% plot(P_signal','LineWidth',.3,'Color',[.8 .8 1])
% hold on
% plot(mean(P_signal),'LineWidth',2,'Color',[.1 .1 .1])
% plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
% title('Noisy signals, its mean (-), and true signal (--)')
% subplot(3,1,2)
% plot(gP_signal','LineWidth',.3,'Color',[1 .8 .8])
% hold on
% plot(mean(gP_signal),'LineWidth',2,'Color',[.1 .1 .1])
% plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
% title('Smoothed noisy signals, its mean (-), and true signal (--)')
% subplot(3,1,3)
% hold on
% plot(ttwsP_signal{1}'.*(exMask(1,:)'*ones(1,Ns)),'LineWidth',.3)
% plot(ttwsP_signal{2}'.*(exMask(2,:)'*ones(1,Ns)),'LineWidth',.3)
% mM_ttws1P_signal = mean(ttwsP_signal{1}).*exMask(1,:); 
% mM_ttws1P_signal(mM_ttws1P_signal==0) = NaN;
% mM_ttws2P_signal = mean(ttwsP_signal{2}).*exMask(2,:); 
% mM_ttws2P_signal(mM_ttws2P_signal==0) = NaN;
% plot(mM_ttws1P_signal,'LineWidth',2,'Color',[.1 .1 .1],'LineStyle','-')
% plot(mM_ttws2P_signal,'LineWidth',2,'Color',[.1 .1 .1],'LineStyle','-')
% % plot(mean(ttwsP_signal{1}).*exMask(1,:),'LineWidth',2,'Color',[.1 .1 .1])
% % plot(mean(ttwsP_signal{2}).*exMask(2,:),'LineWidth',2,'Color',[.1 .1 .1])
% plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
% title('TW-smoothed noisy signals, its mean (-), and true signal (--)')
% set(gcf,'Position',[1000 150 800 1200])

