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
[cP_signal, cP_GmWmCsf, T_names] = cp_create_data(0,[0 0 0]);

%% Check how the mean smoothed signal matches the original signal
% Measure Root Mean Square Error, over each segment and overall.

%% Plot things

% Original signals, noisy and mean, + same but G-smoothed 
figure,

subplot(2,2,1)
plot(P_signal','LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(mean(P_signal),'LineWidth',2,'Color',[.1 .1 .1])
plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
title('Noisy signals, its mean (-), and true signal (--)')

subplot(2,2,2)
plot(gP_signal','LineWidth',.3,'Color',[1 .8 .8])
hold on
plot(mean(gP_signal),'LineWidth',2,'Color',[.1 .1 .1])
plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
title('Smoothed noisy signals, its mean (-), and true signal (--)')

subplot(2,2,3)
hold on
plot(ttwsP_signal{1}'.*(exMask(1,:)'*ones(1,Ns)),'LineWidth',.3,'Color',[.8 .8 1])
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

subplot(4,2,6)
plot(pP_GmWmCsf{1}','LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(mean(pP_GmWmCsf{1}),'LineWidth',2,'Color',[0 0 1])
plot(pP_GmWmCsf{2}','LineWidth',.3,'Color',[1 .8 .8])
plot(mean(pP_GmWmCsf{2}),'LineWidth',2,'Color',[1 0 0])
title('Noisy tissue probabilities and their mean (-)')

subplot(4,2,8)
plot(ggP_GmWmCsf{1}','LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(mean(ggP_GmWmCsf{1}),'LineWidth',2,'Color',[0 0 1])
plot(ggP_GmWmCsf{2}','LineWidth',.3,'Color',[1 .8 .8])
plot(mean(ggP_GmWmCsf{2}),'LineWidth',2,'Color',[1 0 0])
plot(exMask(1,:),'LineWidth',2,'Color',[.1 .1 .1],'LineStyle','--')
plot(exMask(2,:),'LineWidth',2,'Color',[.1 .1 .1],'LineStyle','--')
title('Smoothed noisy tissue prob, their mean (-), and explicit mask (--)')

set(gcf,'Position',[500 150 1600 1200])


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

end