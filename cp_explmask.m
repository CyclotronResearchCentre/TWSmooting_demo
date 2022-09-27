function [exMask] = cp_explmask(gP_GmWmCsf)
% Function to create the explicit tissue mask from the (smoothed) tissue
% probabilities. As in M.F. Callaghan paper, to be in GM a voxel must
% - have a higher probability of being GM than  WM or CSF
% - have a probabiliy larger than 20%
% 
% The input could come from a single subject or multiple ones. In the 
% latter case, then the tissue probabilities are averaged across subjects.
% 
% INPUT
% gP_GmWmCsf : array or cell array of (smoothed) tissue probabilities
%   if [3 x Npt] array of GM, MW, CSF tissue probabilities from 1 subject
%   if {3 x 1} cell array of [Ns x Npt] tissue probabilities (GM, WM, CSF) 
%   from Ns subject
% 
% OUTPUT
% exMask   : [2 x Npt] array with explicit GM and WM mask
%__________________________________________________________________________
% Copyright (C) 2019 Cyclotron Research Centre

% Written by C. Phillips, 2019.
% GIGA Institute, University of Liege, Belgium

%% Deal with input
% Single subject = array data -> put in cell
if ~iscell(gP_GmWmCsf)
    tmp = gP_GmWmCsf;
    gP_GmWmCsf = cell(3,1);
    for ii=1:3
        gP_GmWmCsf{ii} = tmp(ii,:);
    end
    need2avg = 0;
else
    need2avg = 1;
end
Npt = size(gP_GmWmCsf{1},2);

%% Explicit mask
% majority and above 20%

% Check if tissue probabilities need to be averaged across subjects
if need2avg
    avgP_GmWmCsf = zeros(3,Npt);
    for ii=1:3
        avgP_GmWmCsf(ii,:) = mean(gP_GmWmCsf{ii});
    end
else
    avgP_GmWmCsf = zeros(3,Npt);
    for ii=1:3
        avgP_GmWmCsf(ii,:) = gP_GmWmCsf{ii};
    end    
end

exMask = [  ...
    avgP_GmWmCsf(1,:)>avgP_GmWmCsf(2,:) & ... % GM>WM
    avgP_GmWmCsf(1,:)>avgP_GmWmCsf(3,:) & ... % GM>CSF
    avgP_GmWmCsf(1,:)>.2 ; ...                % GM>.2
    avgP_GmWmCsf(2,:)>avgP_GmWmCsf(1,:) & ... % WM>GM
    avgP_GmWmCsf(2,:)>avgP_GmWmCsf(3,:) & ... % WM>CSF
    avgP_GmWmCsf(2,:)>.2 ] ;                  % WM>.2

end