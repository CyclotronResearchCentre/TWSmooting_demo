function cp_groupdemo(opt)
% Demo on synthetic 1D group data.
% 0/ Options, parameters & signals
% 1/ Processing
%       - create data for Ns subjects
%       - smooth them with classic Gaussian or tissue-weighted Gaussian
%       - build explicit masks
%       - create clean signal/tissue probs
% 2/ Calculate match between smoothed signal and clean one
% 3/ Plot the results
% 
% INPUT
% opt   : optional parameters
% .Ns           , number of subjects (20, def)
% .r_jitter     , random jitter introducing some anatomical variability (1, def)
% .smoK         , width of smoothing kernel (8, def)
% .plot_all     , (flag) create the big figs with multi-plots (1) or not (0, def)
% .create_data  , (flag) create the data (1) or try to load them (0, def)
% .fn_data      , filename for saved data
% .save_fig     , (flag) save figures into .png files (1) or not (0, def)
%
% Example call
% - start anew and force data generation then save in 'test.mat' file, 
% everything else being default:
%   cp_groupdemo(struct('create_data',true','fn_data','test.mat'))
% - re-use existing data in szpecific file, generate figure and save them
%   cp_groupdemo(struct('fn_data','test.mat','save_fig',true'))
%__________________________________________________________________________
% Copyright (C) 2019 Cyclotron Research Centre

% Written by C. Phillips, 2019.
% GIGA Institute, University of Liege, Belgium

% Q: should I consider the T-SPOON approach ?
% -> for the sake of completeness, yes!

%% 0/ Options
% default options
opt_def = struct(...
    'Ns', 20, ... % Number of subjects
    'r_jitter', 1, ... % random jitter introducing some anatomical variability
    'smoK', 8, ... % width of smoothing kernel
    'plot_all', false, ... % no big figs with multi-plots
    'create_data', false, ... % try 1st to load the data
    'save_fig', false, ... % not saving figures into file
    'fn_data', 'TWSdata_demo.mat');

% check options & defaults
if nargin==0, opt = []; end
[opt] = crc_check_flag(opt_def,opt);

% get options
fl_plot_all    = opt.plot_all;
fl_create_data = opt.create_data;
fl_save_fig    = opt.save_fig;
fn_data        = opt.fn_data;

% check if data file exist
if ~exist(fn_data,'file') && ~fl_create_data
    fprintf('Can''t find the data file, creating some then!\n')
    fl_create_data = 1;
end

% Parameters & Signals
% --------------------
Ns = opt.Ns; % Number of subjects
r_jitter = opt.r_jitter; % spatial jitter
smoK = opt.smoK; % width of smoothing kernel

% empty output cell arrays
P_GmWmCsf    = cell(Ns,1);
pP_GmWmCsf   = cell(3,1);
gsP_GmWmCsf  = cell(Ns,1);
ggsP_GmWmCsf = cell(3,1);
twsP_signal  = cell(Ns,1);
ttwsP_signal = cell(2,1);

tosP_signal  = cell(Ns,1);
ttosP_signal = cell(2,1);

% Signals
% -------
% Ns / Np     : number of subjects / number of 'pixels'
% P_signal    : signal profiles of Ns subjects, [Ns x Np] array
% P_GmWmCsf   : tissue probabilites profiles of Ns subjects, 
%               {Ns x 1} cell array of [3 x Np] array
% pP_GmWmCsf  : same as P_GmWmCsf but reorganized per tissue class, 
%               as a {3 x 1} cell array of [Ns x Np] array
% gsP_signal  : Gaussian smoothed signals of Ns subjects, [Ns x Np] array
% gsP_GmWmCsf : Gaussian smoothed tissue GM, WM & CSF probabilities of 
%               Ns subjects, {Ns x 1} cell array of [3 x Np] array
% ggsP_GmWmCsf: same as gsP_GmWmCsf but reorganized per tissue class, 
%               as a {3 x 1} cell array of [Ns x Np] array
% twsP_signal : tissue-weighted smoothed signal for GM & WM of 
%               Ns subjects, {Ns x 1} cell array of [2 x Np] array
% ttwsP_signal: same as twsP_signal but reorganized per tissue class,
%               as a {1 x 2} cell array of [Ns x Np] array
% exMask      : group-wise explicit mask, 
%               using majority and >20% for GM & WM, a [2 x Np] array

% tosP_signal : TSPOON smoothed signal for GM & WM of 
%               Ns subjects, {Ns x 1} cell array of [2 x Np] array
% ttosP_signal: same as tosP_signal but reorganized per tissue class,
%               as a {1 x 2} cell array of [Ns x Np] array
% gsP_iexMask  : Gaussian smoothed individual explicit mask for GM & WM of 
%               Ns subjects, {Ns x 1} cell array of [2 x Np] array

%% 1/ Do the processing
% Deal with the Ns subjects, one at a time.
if fl_create_data
    % Create the true clean signal + tissue probabilities
    [cP_signal, cP_GmWmCsf, T_names] = cp_create_data(0,[0 0 0]); %#ok<*ASGLU>

    for ii=Ns:-1:1
        % Create the noisy signal + tissue probabilities, per subject
        [P_signal(ii,:), P_GmWmCsf{ii}] = cp_create_data(r_jitter);
        
        % Smooth the signals, Gaussian & tissue-weighted
        data_ii = struct('P_signal',P_signal(ii,:),'P_GmWmCsf',P_GmWmCsf{ii});
        [gsP_signal(ii,:), gsP_GmWmCsf{ii}, twsP_signal{ii}, ~, tosP_signal{ii}] = ...
            cp_smooth_data(data_ii, smoK);
        
        % Reorganize signals per type rather than per subject
        %   * smoothed tissue classes (-> expl mask) and
        %   * tissue-weighted & TSPOON smoothed signals
        for jj=1:3
            ggsP_GmWmCsf{jj}(ii,:) = gsP_GmWmCsf{ii}(jj,:);
            pP_GmWmCsf{jj}(ii,:) = P_GmWmCsf{ii}(jj,:);
        end
        for jj=1:2
            ttwsP_signal{jj}(ii,:) = twsP_signal{ii}(jj,:);
            ttosP_signal{jj}(ii,:) = tosP_signal{ii}(jj,:);
        end
    end
    % Create explicit mask at group level
    [exMask] = cp_explmask(ggsP_GmWmCsf);
    
    % Save everything useful: options with parameter, tissue probabilities 
    % & names, signal, explicit mask.
    save(fn_data, 'opt', '*_GmWmCsf', '*_signal', 'T_names', 'exMask')
    
else
    load(fn_data)
end

%% 2/ Check how the mean smoothed signal matches the original signal
% Measure Root Mean Square Error, overall and over each segment based on
% explicit mask, for the G-smoothed and TW-smoothed signals, w.r.t. the
% true noise-free signal

% Deal with GM
%-------------
% All GM in explicit mask
l_GMall = find(exMask(1,:));
% Look at difference between true clean signal and
% - mean noisy signal
% - mean G-smoothed noisy signal
% - mean TW-smoothed noisy signal
% - mean TSPOON-smoothed noisy signal

RMSE_GMall = sqrt([...
    sum((cP_signal(l_GMall) - mean(P_signal(:,l_GMall))).^2) / ...
    numel(l_GMall) ;
    sum((cP_signal(l_GMall) - mean(gsP_signal(:,l_GMall))).^2) / ...
    numel(l_GMall) ;
    sum((cP_signal(l_GMall) - mean(ttwsP_signal{1}(:,l_GMall))).^2) / ...
    numel(l_GMall) ; 
    sum((cP_signal(l_GMall) - mean(ttosP_signal{1}(:,l_GMall))).^2) / ...
    numel(l_GMall) ]);
% By segment
l_GMsegmEnd = [0 find(diff(l_GMall)>1) numel(l_GMall)];
Nsegm_GM = numel(l_GMsegmEnd)-1; % number of segmented
Nel_GMSegm = diff(l_GMsegmEnd); % number of voxels/segm
RMSE_GMsegm = zeros(4,Nsegm_GM);
for ii=1:Nsegm_GM
    l_ii = l_GMall(l_GMsegmEnd(ii)+1):l_GMall(l_GMsegmEnd(ii+1));
    RMSE_GMsegm(1,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(P_signal(:,l_ii))).^2) / ...
        Nel_GMSegm(ii) );
    RMSE_GMsegm(2,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(gsP_signal(:,l_ii))).^2) / ...
        Nel_GMSegm(ii) );
    RMSE_GMsegm(3,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(ttwsP_signal{1}(:,l_ii))).^2) / ...
        Nel_GMSegm(ii) );
    RMSE_GMsegm(4,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(ttosP_signal{1}(:,l_ii))).^2) / ...
        Nel_GMSegm(ii) );
end

% Deal with WM
%-------------
% All WM in explicit mask
l_WMall = find(exMask(2,:));
RMSE_WMall = sqrt([...
    sum((cP_signal(l_WMall) - mean(P_signal(:,l_WMall))).^2) / ...
    numel(l_WMall) ;
    sum((cP_signal(l_WMall) - mean(gsP_signal(:,l_WMall))).^2) / ...
    numel(l_WMall) ;
    sum((cP_signal(l_WMall) - mean(ttwsP_signal{2}(:,l_WMall))).^2) / ...
    numel(l_WMall) ;
    sum((cP_signal(l_WMall) - mean(ttosP_signal{2}(:,l_WMall))).^2) / ...
    numel(l_WMall) ]);
% By segment
l_WMsegmEnd = [0 find(diff(l_WMall)>1) numel(l_WMall)];
Nsegm_WM =numel(l_WMsegmEnd)-1; % number of segmented
Nel_WMSegm = diff(l_WMsegmEnd); % number of voxels/segm
RMSE_WMsegm = zeros(4,Nsegm_WM);
for ii=1:Nsegm_WM
    l_ii = l_WMall(l_WMsegmEnd(ii)+1):l_WMall(l_WMsegmEnd(ii+1));
    RMSE_WMsegm(1,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(P_signal(:,l_ii))).^2) / ...
        Nel_WMSegm(ii) );
    RMSE_WMsegm(2,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(gsP_signal(:,l_ii))).^2) / ...
        Nel_WMSegm(ii) );
    RMSE_WMsegm(3,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(ttwsP_signal{2}(:,l_ii))).^2) / ...
        Nel_WMSegm(ii) );
    RMSE_WMsegm(4,ii) = sqrt( ...
        sum((cP_signal(l_ii) - mean(ttosP_signal{2}(:,l_ii))).^2) / ...
        Nel_WMSegm(ii) );
end

% Plot values
% figure,
% plot(Nel_GMSegm,RMSE_GMsegm,'bo',Nel_WMSegm,RMSE_WMsegm,'ro')

% Print out some numbers
fprintf('RMSE over the explicit mask for \n')
fprintf('\tGM signal, no-sm %2.2f, G-sm %2.2f, TW-sm %2.2f and TSP-sm %2.2f\n', ...
    RMSE_GMall(1),RMSE_GMall(2),RMSE_GMall(3),RMSE_GMall(4))
fprintf('\tWM signal, no-sm %2.2f, G-sm %2.2f, TW-sm %2.2f and TSP-sm %2.2f\n', ...
    RMSE_WMall(1),RMSE_WMall(2),RMSE_WMall(3),RMSE_WMall(4))
fprintf('\n')
% fprintf('RMSE ratio\n')
% fprintf('\tGM, no-sm/TW-m %2.2f and G-sm/TW-m %2.2f\n',RMSE_GMall(1)/RMSE_GMall(3),RMSE_GMall(2)/RMSE_GMall(3))
% fprintf('\tWM, no-sm/TW-m %2.2f and G-sm/TW-m %2.2f\n',RMSE_WMall(1)/RMSE_WMall(3),RMSE_WMall(2)/RMSE_GMall(3))

%% 3/ Plot things

% Original signals, noisy and mean, + same but G-smoothed
% Plot things, eiter into a single figure or multiple ones.
% -> use sub-function for each (sub)plot !

% 1/ Plot everything in a single figure
%======================================
if fl_plot_all
    figure,
    subplot(2,2,1)
    plot_Psignal(P_signal,cP_signal)
    subplot(2,2,2)
    plot_gPsignal(gsP_signal,cP_signal)
%     subplot(2,2,3)
%     plot_msktwsPsignal(ttwsP_signal,exMask,cP_signal)
    subplot(4,2,5)
    plot_msktwsPsignal(ttwsP_signal,exMask,cP_signal)
    subplot(4,2,7)
    plot_msktwsPsignal(ttosP_signal,exMask,cP_signal,'TSPOON')
    subplot(4,2,6)
    plot_pPGmWmCsf(pP_GmWmCsf)
    subplot(4,2,8)
    plot_ggPGmWmCsf(ggsP_GmWmCsf,exMask)
    set(gcf,'Position',[500 150 1600 1200])
    if fl_save_fig
        saveas(gcf,'demo_ALL_TissueW_smoothing.png');
    end
end

% 2/ Plot in different figures
%=============================
% Multi-subj noisy signal
% -----------------------
figure, 
plot_Psignal(P_signal,cP_signal)
if fl_save_fig
    saveas(gcf,'demo_OriginalSignal.png');
end

% Gaussian smoothed multi-subj signal
% -----------------------------------
figure,
plot_gPsignal(gsP_signal,cP_signal)
if fl_save_fig
    saveas(gcf,'demo_GsmoothedSignal.png');
end

% Tissue-weighted/TSPOON smoothed multi-subj signal
% -------------------------------------------------
% No masking
figure,
plot_twsPsignal(ttwsP_signal,cP_signal)
if fl_save_fig
    saveas(gcf,'demo_TWsmoothedSignal.png');
end
figure,
plot_twsPsignal(ttosP_signal,cP_signal,'TSPOON')
if fl_save_fig
    saveas(gcf,'demo_TSPsmoothedSignal.png');
end

% With explicit masking
figure,
subplot(2,1,1)
plot_msktwsPsignal(ttwsP_signal,exMask,cP_signal)
subplot(2,1,2)
plot_msktwsPsignal(ttosP_signal,exMask,cP_signal,'TSPOON')
if fl_save_fig
    saveas(gcf,'demo_mskTWsmoothedSignal.png');
end

% Tissue probabilities, original & smoothed
% -----------------------------------------
figure,
subplot(2,1,1)
plot_pPGmWmCsf(pP_GmWmCsf)
subplot(2,1,2)
plot_ggPGmWmCsf(ggsP_GmWmCsf,exMask)
if fl_save_fig
    saveas(gcf,'demo_TissueProb.png');
end

% Signals within the explicit mask, original + G/TW/TSP-smoothed
% --------------------------------------------------------------
figure
subplot(2,1,1)
x_st = 0;
for ii=1:Nsegm_GM
    px = (1:Nel_GMSegm(ii))+x_st;
    x_st = x_st+Nel_GMSegm(ii);
    l_ii = l_GMall(l_GMsegmEnd(ii)+1):l_GMall(l_GMsegmEnd(ii+1));
    plot(px,cP_signal(l_ii),'k--','LineWidth',3)
    hold on
    plot(px,mean(P_signal(:,l_ii)),'b-','LineWidth',1.5)
    plot(px,mean(gsP_signal(:,l_ii)),'color',[.4 .4 1],'LineWidth',1.5)
    plot(px,mean(ttwsP_signal{1}(:,l_ii)),'c-','LineWidth',1.5)
    plot(px,mean(ttosP_signal{1}(:,l_ii)),'color',[.8 .2 1],'LineWidth',1.5)
end
% Add grey lines to seperate segments
y_wd = get(gca,'YLim');
for x_st = [0 cumsum(Nel_GMSegm)]+.5
    plot([x_st x_st],y_wd,'Color',[.5 .5 .5])
end
legend('True signal', ...
    sprintf('Av. noisy signal, RMSE %2.2f',RMSE_GMall(1)), ...
    sprintf('Av. G-smoothed signal, RMSE %2.2f',RMSE_GMall(2)), ...
    sprintf('Av. TW-smoothed signal, RMSE %2.2f',RMSE_GMall(3)), ...
    sprintf('Av. TSP-smoothed signal, RMSE %2.2f',RMSE_GMall(4)), ...
    'Location','NorthWest')
title('GM segments')
subplot(2,1,2)
x_st = 0;
for ii=1:Nsegm_WM
    px = (1:Nel_WMSegm(ii))+x_st;
    x_st = x_st+Nel_WMSegm(ii);
    l_ii = l_WMall(l_WMsegmEnd(ii)+1):l_WMall(l_WMsegmEnd(ii+1));
    plot(px,cP_signal(l_ii),'k--','LineWidth',3)
    hold on
    plot(px,mean(P_signal(:,l_ii)),'r-','LineWidth',1.5)
    plot(px,mean(gsP_signal(:,l_ii)),'color',[1 .8 .2],'LineWidth',1.5)
    plot(px,mean(ttwsP_signal{2}(:,l_ii)),'m-','LineWidth',1.5)
    plot(px,mean(ttosP_signal{2}(:,l_ii)),'color',[1 .5 .6],'LineWidth',1.5)
end
% Add grey lines to seperate segments
y_wd = get(gca,'YLim');
for x_st = [0 cumsum(Nel_WMSegm)]+.5
    plot([x_st x_st],y_wd,'Color',[.5 .5 .5])
end
legend('True signal', ...
    sprintf('Av. noisy signal, RMSE %2.2f',RMSE_WMall(1)), ...
    sprintf('Av. G-smoothed signal, RMSE %2.2f',RMSE_WMall(2)), ...
    sprintf('Av. TW-smoothed signal, RMSE %2.2f',RMSE_WMall(3)), ...
    sprintf('Av. TSP-smoothed signal, RMSE %2.2f',RMSE_WMall(4)), ...
    'Location','SouthWest')
title('WM segments')
set(gcf,'Position',[600 120 500 800])

if fl_save_fig
    saveas(gcf,'demo_RMSE_segments.png');
end

% Same figure but with TW and TSPOON smoothiung only
figure
subplot(2,1,1)
x_st = 0;
for ii=1:Nsegm_GM
    px = (1:Nel_GMSegm(ii))+x_st;
    x_st = x_st+Nel_GMSegm(ii);
    l_ii = l_GMall(l_GMsegmEnd(ii)+1):l_GMall(l_GMsegmEnd(ii+1));
    plot(px,cP_signal(l_ii),'k--','LineWidth',3)
    hold on
    plot(px,mean(ttwsP_signal{1}(:,l_ii)),'c-','LineWidth',1.5)
    plot(px,mean(ttosP_signal{1}(:,l_ii)),'color',[.8 .2 1],'LineWidth',1.5)
end
% Add grey lines to seperate segments
y_wd = get(gca,'YLim');
for x_st = [0 cumsum(Nel_GMSegm)]+.5
    plot([x_st x_st],y_wd,'Color',[.5 .5 .5])
end
legend('True signal', ...
    sprintf('Av. TW-smoothed signal, RMSE %2.2f',RMSE_GMall(3)), ...
    sprintf('Av. TSP-smoothed signal, RMSE %2.2f',RMSE_GMall(4)), ...
    'Location','NorthWest')
title('GM segments')
subplot(2,1,2)
x_st = 0;
for ii=1:Nsegm_WM
    px = (1:Nel_WMSegm(ii))+x_st;
    x_st = x_st+Nel_WMSegm(ii);
    l_ii = l_WMall(l_WMsegmEnd(ii)+1):l_WMall(l_WMsegmEnd(ii+1));
    plot(px,cP_signal(l_ii),'k--','LineWidth',3)
    hold on
    plot(px,mean(ttwsP_signal{2}(:,l_ii)),'m-','LineWidth',1.5)
    plot(px,mean(ttosP_signal{2}(:,l_ii)),'color',[1 .5 .6],'LineWidth',1.5)
end
% Add grey lines to seperate segments
y_wd = get(gca,'YLim');
for x_st = [0 cumsum(Nel_WMSegm)]+.5
    plot([x_st x_st],y_wd,'Color',[.5 .5 .5])
end
legend('True signal', ...
    sprintf('Av. TW-smoothed signal, RMSE %2.2f',RMSE_WMall(3)), ...
    sprintf('Av. TSP-smoothed signal, RMSE %2.2f',RMSE_WMall(4)), ...
    'Location','SouthWest')
title('WM segments')
set(gcf,'Position',[600 120 500 800])

if fl_save_fig
    saveas(gcf,'demo_RMSE_segments_TWTSP.png');
end


end
%% PLOTING SUB-FUNCTIONS

function plot_Psignal(P_signal,cP_signal)

plot(P_signal','LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(mean(P_signal),'LineWidth',2,'Color',[.1 .1 .1])
plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
title('Noisy signals, its mean (-), and true signal (--)')

end

function plot_gPsignal(gsP_signal,cP_signal)

plot(gsP_signal','LineWidth',.3,'Color',[1 .8 .8])
hold on
plot(mean(gsP_signal),'LineWidth',2,'Color',[.1 .1 .1])
plot(cP_signal,'LineWidth',2,'Color',[.5 .5 .5],'LineStyle','--')
title('Smoothed noisy signals, its mean (-), and true signal (--)')

end

function plot_twsPsignal(ttwsP_signal, cP_signal, meth_lab)
% Display the tissue-weighted/TSPOON smoothed signal for GM and WM in
% different subplots

% Define method label, typically TW (tissue weighted)
if nargin<3 || isempty(meth_lab), meth_lab = 'TW'; end

subplot(2,1,1)
plot(ttwsP_signal{1}','LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(mean(ttwsP_signal{1}),'LineWidth',2,'Color',[0 0 1])
plot(cP_signal,'LineWidth',2,'Color',[.2 .2 .2],'LineStyle','--')
title([meth_lab,' smoothed noisy signals for GM and true signal (--)'])
axis([0 200 -20 120])

subplot(2,1,2)
plot(ttwsP_signal{2}','LineWidth',.3,'Color',[1 .8 .8])
hold on
plot(mean(ttwsP_signal{2}),'LineWidth',2,'Color',[1 0 0])
plot(cP_signal,'LineWidth',2,'Color',[.2 .2 .2],'LineStyle','--')
title([meth_lab,' smoothed noisy signals for WM and true signal (--)'])
axis([0 200 -20 120])

end

function plot_msktwsPsignal(ttwsP_signal, exMask, cP_signal, meth_lab)
% Display the tissue-weighted/TSPOON smoothed signal per segment in a 
% single plot

% Define method label, typically TW (tissue weighted)
if nargin<4 || isempty(meth_lab), meth_lab = 'TW'; end

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
title(['Masked ',meth_lab,' smoothed signals, GM (blue) and WM (red), ', ...
    'mean (-), and true signal (--)'])
end

function plot_pPGmWmCsf(pP_GmWmCsf)

plot(pP_GmWmCsf{1}','LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(mean(pP_GmWmCsf{1}),'LineWidth',2,'Color',[0 0 1])
plot(pP_GmWmCsf{2}','LineWidth',.3,'Color',[1 .8 .8])
plot(mean(pP_GmWmCsf{2}),'LineWidth',2,'Color',[1 0 0])
title('Noisy tissue probabilities, GM (blue) and WM (red), and their mean (-)')

end

function plot_ggPGmWmCsf(ggsP_GmWmCsf,exMask)

plot(ggsP_GmWmCsf{1}','LineWidth',.3,'Color',[.8 .8 1])
hold on
plot(mean(ggsP_GmWmCsf{1}),'LineWidth',2,'Color',[0 0 1])
plot(ggsP_GmWmCsf{2}','LineWidth',.3,'Color',[1 .8 .8])
plot(mean(ggsP_GmWmCsf{2}),'LineWidth',2,'Color',[1 0 0])
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
% plot(gsP_signal','LineWidth',.3,'Color',[1 .8 .8])
% hold on
% plot(mean(gsP_signal),'LineWidth',2,'Color',[.1 .1 .1])
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
% plot(ggsP_GmWmCsf{1}','LineWidth',.3,'Color',[.8 .8 1])
% hold on
% plot(mean(ggsP_GmWmCsf{1}),'LineWidth',2,'Color',[0 0 1])
% plot(ggsP_GmWmCsf{2}','LineWidth',.3,'Color',[1 .8 .8])
% plot(mean(ggsP_GmWmCsf{2}),'LineWidth',2,'Color',[1 0 0])
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
% plot(gsP_signal','LineWidth',.3,'Color',[1 .8 .8])
% hold on
% plot(mean(gsP_signal),'LineWidth',2,'Color',[.1 .1 .1])
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

