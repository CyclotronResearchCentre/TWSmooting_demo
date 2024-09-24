function [gsP_signal, twsP_signal, tosP_signal] = aj_smooth_1d(data_1D, data_GmWmCsf, sm_kern, plot_fig)
% aj_smooth_1d Perform 3 smoothing methods (Gaussian, tissue-weighted, TSPOON) on 1D data.
% 
% INPUTS:
% data_1D      : 1D data profile to smooth (signal) [1 x N]
% data_GmWmCsf : Tissue probabilities for GM, WM, and CSF [3 x N]
% sm_kern      : Smoothing kernel width (default = 8)
% plot_fig     : Flag to plot comparison of the methods (default = false)
%
% OUTPUTS:
% gsP_signal   : Smoothed signal using Gaussian smoothing
% twsP_signal  : Tissue-weighted smoothed signal (GM and WM)
% tosP_signal  : TSPOON smoothed signal (GM and WM)

%% Default parameters
if nargin < 4, plot_fig = false; end
if nargin < 3, sm_kern = 8; end  % Default kernel size

%% Apply Gaussian smoothing
wg = gausswin(sm_kern);
wg = wg / sum(wg); % Normalize the kernel

gsP_signal = filtfilt(wg, 1, data_1D); % Gaussian smoothing on signal

%% Tissue-weighted smoothing
% Using a wider Gaussian kernel (2x width) for tissue probability maps
wg2 = gausswin(2 * sm_kern);
wg2 = wg2 / sum(wg2); % Normalize the kernel

% Initialize the tissue-weighted signal (for GM and WM)
twsP_signal = zeros(2, length(data_1D));

for ii = 1:2  % Only use GM and WM (ignoring CSF)
    tmp1 = data_1D .* data_GmWmCsf(ii,:) .* (filtfilt(wg2, 1, data_GmWmCsf(ii,:)) > 0.05);
    twsP_signal(ii,:) = filtfilt(wg, 1, tmp1) ./ filtfilt(wg, 1, data_GmWmCsf(ii,:)) .* (data_GmWmCsf(ii,:) > 0.05);
end

%% TSPOON smoothing
% Create explicit mask based on majority vote and 20% threshold for GM and WM
iexMask = double(data_GmWmCsf(1,:) > 0.2 | data_GmWmCsf(2,:) > 0.2);

% Smooth the explicit mask with Gaussian kernel
gsP_iexMask = filtfilt(wg, 1, iexMask);

% Initialize the TSPOON signal
tosP_signal = zeros(2, length(data_1D));

for ii = 1:2  % For GM and WM
    tosP_signal(ii,:) = filtfilt(wg, 1, iexMask .* data_1D) ./ gsP_iexMask;
end

%% Plot comparison of smoothing methods if requested
if plot_fig
    figure;
    subplot(3,1,1); plot(data_1D, 'k'); hold on;
    plot(gsP_signal, 'r'); title('Gaussian Smoothing'); legend('Original', 'Smoothed');
    
    subplot(3,1,2); plot(twsP_signal(1,:), 'r'); hold on;
    plot(twsP_signal(2,:), 'b'); title('Tissue-Weighted Smoothing');
    legend('GM-weighted', 'WM-weighted');
    
    subplot(3,1,3); plot(tosP_signal(1,:), 'r'); hold on;
    plot(tosP_signal(2,:), 'b'); title('TSPOON Smoothing');
    legend('GM-smoothed', 'WM-smoothed');
end

end
