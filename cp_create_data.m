function [P_signal, P_GmWmCsf,T_names] = cp_create_data(r_jitter,T_noise, plot_fig)
% Function to create some synthetic 1D data:
% - some "image" values, supposed to be like T1w or MT
% - some tissue probabilities for GM, WM & CSF
%
% INPUT:
% r_jitter : random jitter to add to the profile indexes, i.e. introducing
%            some pseudo-anatomocal variability. Def. r_jitter=0
% T_noise  : noise variance for the signal. Def. T_noise = [2 2 10]
% plot_fig : flag to plot or not some figures Def. 0 = no plot
%
% OUTPUT:
% P_signal  : signal profile, [1 x N] array
% P_GmWmCsf : tissue probabilites profiles, [3 x N] array
%__________________________________________________________________________
% Copyright (C) 2019 GIGA Institute

% Written by C. Phillips, 2019.
% Cyclotron Research Centre, University of Liege, Belgium

%% Key parameters
if nargin<3, plot_fig = 0; end
if nargin<2, T_noise  = [2 2 10]; end
if nargin<0, r_jitter=0; end


% for GM, WM & CSF respectively
T_names  = {'GM', 'WM', 'CSF'};
T_signal = [50 100 5];
T_Wsegm  = [32 16 24 26 24 8 12 8 12 4 32]; % width of tissue segments
T_probGmWmCsf = [ ...
    1 98  2  2  2 95  5  2  5 98  1 ; ...
    1  1 97  2 97  2 94  4 94  1  1 ; ...
    98  1  1 96  1  3  1 95  1  1 98]; % tissue prob per segment
T_tcnoise = 2; %error on tissue class probability
T_tcthresh = .1; % minimum tissue class prob in % !

T_index = cumsum([1 T_Wsegm]); % T_index(end) = T_index(end)+1;
Nsegm = numel(T_Wsegm);
% Adding anatomical variability, if requested
if r_jitter
    i_jitter = floor(rand(1,numel(T_index))*(r_jitter*2+1) - r_jitter);
    i_jitter(1) = 0; i_jitter(end) = 0; % not moving signal border
    T_index = T_index + i_jitter;
end

%% Creating data and tissue profiles
P_GmWmCsf = zeros(3,max(T_index));
P_signal  = zeros(1,max(T_index));
for ii=1:Nsegm
    l_ind = T_index(ii):T_index(ii+1)-1; % indexes
    for jj=1:3 % 3 tissue probabily profiles
        % exact tissue prob
        P_GmWmCsf(jj,l_ind) = T_probGmWmCsf(jj,ii);
        % signal from tissue prob + noise
        P_signal(l_ind) = P_signal(l_ind) + ...
            T_probGmWmCsf(jj,ii)/100*T_signal(jj) + ...
            T_probGmWmCsf(jj,ii)/100.*randn(1,numel(l_ind))*T_noise(jj) ;
    end
end
% add noise to tissue prob -> need to ensure p>0 and sum==1
P_GmWmCsf = P_GmWmCsf + randn(size(P_GmWmCsf))*T_tcnoise;
% minimum value should not be lower than .1%
P_GmWmCsf(P_GmWmCsf<T_tcthresh) = T_tcthresh;
% sum to 1
P_GmWmCsf = P_GmWmCsf ./ (ones(3,1)*sum(P_GmWmCsf)) ;

%% Plot all profiles
if plot_fig
    figure,
    % display intensity profile
    subplot(2,1,1)
    plot(P_signal)
    ylabel('Intensities')
    
    % display tissue probability profile
    for ii=1:3
        subplot(6,1,3+ii)
        plot(P_GmWmCsf(ii,:))
        ylabel(T_names{ii})
    end
end

end
