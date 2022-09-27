function [gsP_signal,gsP_GmWmCsf,twsP_signal,exMask]  = ...
    cp_smooth_data(data,sm_kern,plot_fig)
% Function to smooth the synthetic 1D data, using
% - the standard Gaussian smoothing kernel
% - the tissue weighted smoothing
% also returning the explicit mask as defined for VBQ/hMRI analysis
%
% INPUT:
% sm_kern  : smoothing kernel width [def. 8]
% data     : structure with the input signal
%   .P_signal  : (noisy) signal to smooth (a [1xN] array)
%   .P_GmWmCsf : tissue probabilities for GM, WM & CSF (a [3 x N] array)
%   .T_names   : tissue classes name, only used for plotting
% plot_fig : flag to plot or not some figures [def. 0 = no plot]
%
% OUTPUT:
% gsP_signal   : smoothed signal, with usual Gaussian kernel, [1 x N] array
% gsP_GmWmCsf  : smoothed tissue GM, WM & CSF probabilities, [3 x N] array
% twsP_signal : tissue-weighted smoothed signal for GM & WM, [2 x N] array
% exMask      : explicit mask, using majority and >20% for GM & WM
%               Just for that single subject -> not group level!
%__________________________________________________________________________
% Copyright (C) 2019 Cyclotron Research Centre

% Written by C. Phillips, 2019.
% GIGA Institute, University of Liege, Belgium

%% Some parameters
if nargin<3, plot_fig = 0; end
if nargin<2, sm_kern = []; end % -> use default
if isempty(sm_kern)
    sm_kern = 8; % default smoothing kernel size
end

%% Get the data
P_signal  = data.P_signal;
P_GmWmCsf = data.P_GmWmCsf;
if plot_fig
    T_names   = data.T_names;
end

%% Apply standard smoothing
wg = gausswin(sm_kern);
wg = wg/sum(wg); % normalize

gsP_signal = filtfilt(wg,1,P_signal);
gsP_GmWmCsf = filtfilt(wg,1,P_GmWmCsf')';

%% Plot all profiles
if plot_fig
    figure,
    % display intensity profile
    subplot(2,1,1)
    plot(gsP_signal)
    ylabel('Intensities')
    
    % display tissue probability profile
    for ii=1:3
        subplot(6,1,3+ii)
        plot(gsP_GmWmCsf(ii,:))
        ylabel(T_names{ii})
    end
    
    % Plot orginal and G-smoothed signal
    figure,
    plot(P_signal), hold on
    plot(gsP_signal,'r')
end

%% Appli tissue-weighted smoothing

% Appplying the smoothing as implemented for VBQ,
% assuming the TPMs are like the tissue probability but smoothed with a
% kernel twice the size, for simplicity.
for ii=1:2
    tmp1 = P_signal .* P_GmWmCsf(ii,:) .* ...
        (filtfilt(2*wg,1,P_GmWmCsf(ii,:))>.05); % Like the TPM masking
    twsP_signal(ii,:) = filtfilt(wg,1,tmp1) ./ gsP_GmWmCsf(ii,:) .* ...
        (gsP_GmWmCsf(ii,:)>.05); %#ok<*AGROW> % masking from smoothed tissue
end

if plot_fig
    figure,
    plot(P_signal),
    hold on
    plot(twsP_signal(1,:),'r')
    plot(twsP_signal(2,:),'c')
end

%% Explicit mask
% majority and above 20%
exMask = cp_explmask(gsP_GmWmCsf);

if plot_fig
    figure,
    plot(P_signal),
    hold on
    plot(twsP_signal(1,:).*exMask(1,:),'r')
    plot(twsP_signal(2,:).*exMask(2,:),'c')
end

end