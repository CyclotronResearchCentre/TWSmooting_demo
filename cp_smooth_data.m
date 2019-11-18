% function out = cp_smooth_data
% %__________________________________________________________________________
% % Copyright (C) 2019 GIGA Institute
% 
% % Written by C. Phillips, 2019.
% % Cyclotron Research Centre, University of Liege, Belgium

%% Some parameters
sm_kern = 8; % smoothing kernel size

%% Get the data
[P_signal, P_GmWmCsf, T_names] = cp_create_data;

%% Apply standard smoothing
wg = gausswin(sm_kern);
wg = wg/sum(wg); % normalize

fP_signal = filtfilt(wg,1,P_signal);
fP_GmWmCsf = filtfilt(wg,1,P_GmWmCsf')';

%% Plot all profiles
figure,
% display intensity profile
subplot(2,1,1)
plot(fP_signal)
ylabel('Intensities')

% display tissue probability profile
for ii=1:3
    subplot(6,1,3+ii)
    plot(fP_GmWmCsf(ii,:))
    ylabel(T_names{ii})
end

% Plot orginal and G-smoothed signal
figure,
plot(P_signal), hold on
plot(fP_signal,'r')

%% Appli tissue-weighted smoothing

% Appplying the smoothing as implemented for VBQ,
% assuming the TPMs are like the twice smoothed tissue prob for simplicity.
for ii=1:2
    tmp1 = P_signal .* P_GmWmCsf(ii,:) .* ...
        (filtfilt(wg,1,fP_GmWmCsf(ii,:))>.05); % Like the TPM masking
    twsP_signal(ii,:) = filtfilt(wg,1,tmp1) ./ fP_GmWmCsf(ii,:) .* ...
                    (fP_GmWmCsf(ii,:)>.05); % masking from smoothed tissue
end
%#ok<*SAGROW>


figure,
plot(P_signal), 
hold on
plot(twsP_signal(1,:),'r')
plot(twsP_signal(2,:),'c')

%% Explicit mask
% majority and above 20%

exMask = [  fP_GmWmCsf(1,:)>fP_GmWmCsf(2,:) & ...
            fP_GmWmCsf(1,:)>fP_GmWmCsf(3,:) & ...
            fP_GmWmCsf(1,:)>.2 ; ...
            fP_GmWmCsf(2,:)>fP_GmWmCsf(1,:) & ...
            fP_GmWmCsf(2,:)>fP_GmWmCsf(3,:) & ...
            fP_GmWmCsf(2,:)>.2 ] ; 

figure,
plot(P_signal), 
hold on
plot(twsP_signal(1,:).*exMask(1,:),'r')
plot(twsP_signal(2,:).*exMask(2,:),'c')
