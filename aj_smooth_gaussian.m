function [gs_imgaussfilt3_paths, gs_spm_paths] = aj_smooth_gaussian(MPM_paths, pth_out, fwhm, dim)
%--------------------------------------------------------------------------
% Function to apply standard Gaussian smoothing to 1D, 2D or 3D data.
%
% INPUTS
% con_info: A structure array containing image volume information (from spm_vol) 
%           of data to be smoothed
% pth_out:  Directory to save output files
% param:    Smoothing parameters, including the width of smoothing kernel
%           in mm (from default)
% dim:      Dimension of the data (1D, 2D or 3D)
%
% OUTPUTS
% gsP_signal:       Gaussian smoothed signal from imgaussfilt3
% gs_path:          Gaussian smoothed file path from spm_smooth
%
% REFERENCE ???
%
%--------------------------------------------------------------------------
% 1D GAUSSIAN SMOOTHING
% filtfilt: Zero-phase forward and reverse digital IIR filtering.
%
% gsP_signal = filtfilt(wg, 1, Y) filters the data in vector Y with the 
% filter described by vectors wg and 1 to create the filtered data gsP_signal.
% LIMITATION: filtfilt function (double filtering low-pass filtering) 
% requires that the length of the data is greater than three times the 
% length of the filter core.
%
%--------------------------------------------------------------------------
% 2D GAUSSIAN SMOOTHING
% imgaussfilt: 2-D Gaussian filtering of images.
%
% gsP_signal = imgaussfilt(Y, 'FilterSize', param.sm_kern_gaussian)
% filters image Y with a 2D Gaussian smoothing kernel with standard 
% deviation of 0.5 which can be modified by using imgaussfilt(A,SIGMA,...).
% 'FilterSize' is a scalar or 2-element vector, of positive, odd integers 
% that specifies the size of the Gaussian filter.
%
%--------------------------------------------------------------------------
% 3D GAUSSIAN SMOOTHING
% imgaussfilt3: 3-D Gaussian filtering of 3-D images.
%
% gsP_signal = imgaussfilt3(Y, 'FilterSize', param.sm_kern_gaussian)
% filters 3-D image Y with a 3D Gaussian smoothing kernel with standard 
% deviation of 0.5 which can be modified by using imgaussfilt3(A,SIGMA,...).
% 'FilterSize' is a scalar or 3-element vector, of positive, odd integers 
% that specifies the size of the Gaussian filter.
%
% FORMAT spm_smooth(P,Q,s,dtype): 3 dimensional convolution of an image.
% P     - image(s) to be smoothed (or 3D array)
% Q     - filename for smoothed image (or 3D array)
% s     - [sx sy sz] Gaussian filter width {FWHM} in mm (or edges)
% dtype - datatype [Default: 0 == same datatype as P]
%--------------------------------------------------------------------------
% NOTES
% imgaussfilt3 tends to propagate Nan values >< spm_smooth where not finite
% values are fixed to 0.0 in spm_con_vol.c (i.e. lines 22, 23).
% https://github.com/spm/spm/blob/main/src/spm_conv_vol.c
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------

if nargin < 4
    warning('GS issue: not enough inputs.');
end

% if isvector(MPM_paths)
%     nb_MPM = length(MPM_paths);
% else
%     warning('GS issue: con_info is not a vector (i.e. not specific to a subject).\n');
% end

% Generate Gaussian kernel
wg = gausswin(fwhm);
wg = wg / sum(wg);


% Apply standard Gaussian smoothing based on data dimensions
switch dim
    case 1
        con = spm_read_vols(MPM_paths);
        gsP_signal = filtfilt(wg, 1, con);

    case 2
        con = spm_read_vols(MPM_paths);
        gsP_signal = imgaussfilt(con, 'FilterSize', fwhm);

    case 3
        nMPM = size(MPM_paths,1);
        gs_spm_paths = cell(nMPM,1);
        gs_imgaussfilt3_paths = cell(nMPM,1);
        for i = 1:nMPM
            gs_spm_paths{i} = spm_file(MPM_paths(i,:), 'prefix', 'spm_smooth_', 'path', pth_out);
            spm_smooth(MPM_paths(i,:), gs_spm_paths{i}, fwhm);
%             hmri_proc_zero2nan(gs_spm_paths{i});
%             aj_proc_thr2nan(gs_spm_paths{i}, 0.05);

            % test imgaussfilt3, a matlab function => NaN absorption issue
%             gsP_signal = imgaussfilt3(spm_read_vols(con_info{i}), 'FilterSize', param.sm_kern_gaussian);
%             V_info = con_info{i}; % Copy the structure of the original file header
%             output_filename = spm_file(V_info.fname, 'path', pth_out, 'prefix', ['imgaussfilt3', '_']);
%             V_info.fname = output_filename;
%             spm_write_vol(V_info, gsP_signal);
%             gs_imgaussfilt3_paths{i} = V_info.fname;
        end

    otherwise
        warning('Data dimension (dim) is neither 1D, 2D nor 3D.');
end
end