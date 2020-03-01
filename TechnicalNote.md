# Tissue-weighted smoothing evaluation

[TOC]
---

## Introduction

After being warped into a common space, typically a group averaged in line with the MNI space, MR images are smoothed. This is useful to : 1/ remove some high-frequency noise, 2/ reduce the multiple comparison problem, and 3/ reduce the remaining inter-subject anatomical variability.

With quantitative MRI, one is interested in the signal coming from the GM and/or the WM, standard Gaussian smoothing would therefore mix the signal coming from these 2 tissues classes, i.e. introducing some partial volume effect! One would therefore want to smooth the images "within each tissue class". Such an approach was introduced by Drangaski et al. (2011) as a specific way to smooth data in the case of "Voxel-Based Quantification" (VBQ) analysis. 

Another technique is the "*T*issue-*SP*ecific, sm*OO*thing-compe*N*sated" method, aka. T-SPOON, by Eun Lee et al. It is not considered here, at least for the moment...

---
## VBQ methods

Following spatial warping, the last 2 spatial processing steps are smoothing, using the a tissue-weighted approach) and the creation of explicit masks. 

### Tissue-weighted smoothing
For each tissue class of interest, typically GM and WM, the quantitative map is smoothed according to the tissue class *posterior* probability. Tissue-weighted smoothing is thus defined as follow:

<img src="https://latex.codecogs.com/gif.latex?p_j=\frac{g*(w_js_j)}{g*w_j}\:m_{\tiny\mbox{TPM}}\:m_j" />

where:

| Parameter | Meaning  |
| ------ | -------- |
| <img src="https://latex.codecogs.com/gif.latex?p_j" /> | Quantitative map for subject <img src="https://latex.codecogs.com/gif.latex?j" /> after tissue-weighted smoothing |
| <img src="https://latex.codecogs.com/gif.latex?s_j" /> | Participant-specific quantitative map warped to group space by deformation <img src="https://latex.codecogs.com/gif.latex?\Phi_j" /> |
| <img src="https://latex.codecogs.com/gif.latex?\Phi_j" /> | Participant-specific deformation mapping from native to group space |
| <img src="https://latex.codecogs.com/gif.latex?w_j" /> | Participant-specific tissue weights given by <img src="https://latex.codecogs.com/gif.latex?\|J_j\| t_j" /> |
| <img src="https://latex.codecogs.com/gif.latex?\|J_j\|" /> | Jacobian determinants of deformation <img src="https://latex.codecogs.com/gif.latex?\Phi_j" />              |
| <img src="https://latex.codecogs.com/gif.latex?t_j" /> | Participant-specific tissue *posterior* probability map warped by deformation <img src="https://latex.codecogs.com/gif.latex?\Phi_j" /> |
| <img src="https://latex.codecogs.com/gif.latex?g*" /> | Convolution by a Gaussian smoothing kernel, i.e. Gaussian smoothing. |
| <img src="https://latex.codecogs.com/gif.latex?m_{\tiny\mbox{TPM}}" /> | TPM-specific mask identifying voxels with probability > 5%   |
| <img src="https://latex.codecogs.com/gif.latex?m_j" /> | Participant-specific mask defined as <img src="https://latex.codecogs.com/gif.latex?g*w_j>5\%" /> |
| <img src="https://latex.codecogs.com/gif.latex?./." /> | Ratio applied voxel by voxel over the images |

The point of the 2 masks is to ensure only voxels with sufficient 

- *a priori* probability (<img src="https://latex.codecogs.com/gif.latex?m_{\tiny\mbox{TPM}}" />) of being of the tissue class of interest, i.e. keeping voxels that are where they are supposed to be, and 
- *a posteriori* probability (<img src="https://latex.codecogs.com/gif.latex?m_j" />) of being of the tissue class of interest, i.e. masking out voxels that are unlikely to be of the class of interest 

A further explicit mask should be defined at the group level for the statistical analysis.

### Explicit mask

For each tissue class of interest, i.e. GM and WM, an explicit mask is generated from the group-mean smoothed tissue probabilities.  The GM (resp. WM) mask includes voxels having, on average across all subjects, 

- \> 20% chance of being GM (resp. WM), and 
- larger probability of being GM (resp. WM) than CSF or WM (resp. GM), 

This explicit mask ensures that for a voxel will only appear in one tissue class, either GM or WM, the one with the largest probability (on average across the group)

---

## Simulation

To illustrate how classic and tissue-weighted smoothing affect the quantitative signal, we created synthetic data in 1D, for an easy visualization, for the 20 pseudo-subjects. Then the average across subject is displayed.

### Data generation

Here are the characteristics of the simulated data for each of the 20 subjects:

- the profile contains segments of 3 tissue classes, say GM, WM, and CSF, with true intensity of 50, 100, and 5 resp., see the 1st figure here under;
- some "anatomical variability" is introduced by randomly shifting the edge of all the segments of the profile by -1, 0 or +1;
- each "tissue segment" has quite distinct probabilities, i.e.  >94% or <5%, but they are always  >1% and the total at each voxel =100%;
- the signal is constructed by summing over the 3 tissue classes the product of the tissue probability with their corresponding intensity, and adding some random noise. The standard deviation of noise for the GM, WM, and CSF is 2, 2, 10 resp.;
- some randomness is also added to the tissue class profiles but they  still remain >1% and the total at each voxel =100%, see last figure for their profile.

There are therefore different sources of variability in this simulation:

- anatomical, with edges being shifted left or right randomly;
- intensity, with some added noise in the signal;
- tissue class, with some added noise in the probability.

### Data smoothing & averaging

 The width of the Gaussian smoothing kernel is arbitrarily set to 8 voxels here. For each subject, 

- the tissue probabilities are smoothed using the same Gaussian kernel.
- the signal is smoothed using a (standard) Gaussian kernel and using the (VBQ) tissue-weighted method for the GM and WM.

Then these signals, original and smoothed, and tissue probabilities are averaged over subjects. These would thus represent the group mean (smoothed) signals and tissue probabilities.

### Limitations of the simulation

There are 2 obvious limitations with this simulation:

- there is no volume modulation by the Jacobian determinant of the spatial deformation. Here this term is assumed to be 1.
- there is no *a priori* tissue probability maps (TPMs). Here we simply consider the true tissue probabilities and smooth them with a kernel twice as wide as what is used for the signal, i.e. 16 voxels, as they tend to be pretty smooth 

Further simulation could lift these limitations and explicitly test for the effect of the volume modulation and the TPMs.

---

## Evaluations

### Plots

Different versions of  signals and  tissue probabilities profiles are displayed

- true and noisy,
- subject individual and averaged, 
- smoothed with Gaussian kernel and tissue-weighted

### Discrepancy measurement

After averaging over the 20 subjects, one would expect that the resulting mean (smoothed) signals and tissue probabilities are close to the true underlying profile. The discrepancy is measured as the "Root Mean Square Error" (RMSE). This is calculated separately for the signal in the GM and WM, over their respective explicit mask, for the cases where 

- no smoothing is applied; 
- the signal is smoothed with the standard Gaussian and the tissue-weighted method.


<img src="https://latex.codecogs.com/gif.latex?\mbox{RMSE}_{TC}=\sqrt{\frac{\sum_{i\in\mbox{ExplMsk}_{TC}}(tS_i-sS_i)^2}{\vert\mbox{ExplMsk}_{TC}\vert}}" />

where

| Parameter          | Meaning                                      |
| ------------------ | -------------------------------------------- |
| <img src="https://latex.codecogs.com/gif.latex?\mbox{RMSE}_{TC}" /> | Root Mean Square Error for tissue class <img src="https://latex.codecogs.com/gif.latex?TC" /> |
| <img src="https://latex.codecogs.com/gif.latex?\mbox{ExplMsk}_{TC}" /> | Explicit mask for tissue class <img src="https://latex.codecogs.com/gif.latex?TC" /> |
| <img src="https://latex.codecogs.com/gif.latex?tS_i" /> | True signal at voxel <img src="https://latex.codecogs.com/gif.latex?i" /> |
| <img src="https://latex.codecogs.com/gif.latex?sS_i" /> | Averaged (over 20 subjects) signal with/without smoothing at voxel <img src="https://latex.codecogs.com/gif.latex?i" /> |
| <img src="https://latex.codecogs.com/gif.latex?\vert\mbox{ExplMsk}_{TC}\vert" /> | Number of voxels in explicit mask for tissue class <img src="https://latex.codecogs.com/gif.latex?TC" /> |

---

## Results

- Signal from the 20 subjects, thin lines, and the average signal intensity bold black line. The true underlying signal is represented by the dashed-grey line. From left to right, each segment contains the following tissue with width (vx): CSF (24), GM (24), WM(24), CSF (26), WM (24), GM (12), WM (8), CSF (12), WM (12), GM (6), CSF (26). Some segments are quite broader than the smoothing kernel (8) but others are of equivalent size (12 or 8) or even a bit smaller (6)
<img src="demo_OriginalSignal.png" style="zoom: 150%;" />
  One can see the variability in signal intensity and some "fussiness" close to the edges. Even when averaging the noisy signal from the 20 subjects, there remains some variability.
  
- Tissue probabilities, with noise (top) and after Gaussian smoothing (bottom), with the corresponding explicit mask (bottom, dashed line) with, in blue, the GM and, in red, the WM tissue. .
	<img src="demo_TissueProb.png" style="zoom: 150%;" />
  This simply illustrates the fact that the tissue probabilities are also affected by the standard Gaussian smoothing, still the explicit mask ensures that the location of the underlying tissue classes is correctly recovered.

- The signal from the 20 subjects, thin lines, and the mean signal intensity, bold line, smoothed with the Gaussian kernel.
	<img src="demo_GsmoothedSignal.png" style="zoom: 150%;" />
  One can see the mixing effect of standard Gaussian smoothing. The signal close to the edges (between tissue classes) is strongly affected and deviates from the true signal (dashed line).
  
- The signal from the 20 subjects, thin lines, and the mean signal intensity, bold line, smoothed with the tissue-weighted method for the GM (blue top) and WM (red bottom) tissues, plus the true signal (black dashed line). 
	<img src="demo_TWsmoothedSignal.png" style="zoom: 150%;" />
  The signal is fairly flat over the actual tissue segment and matches closely the true signal. Because of the smoothing kernel width, the signal also extends into the neighbouring segments and the value outside the tissue segment is dragged up or down, depending on the adjourning signal.
- Same as the previous figure but after applying the explicit masking for GM and WM and plotting both together, plus the true signal (black dashed line). 
	<img src="demo_mskTWsmoothedSignal.png" style="zoom: 150%;" />
  One can notice that the the signal is relatively homogeneously smoothed within each tissue class and only deviates slightly from the true signal very close to the edges (between tissue classes). See here under for a zoom in on each tissue segment.
  
- Same as the previous figure but zooming in on the true and average signals over the explicit mask segments, for GM (top) and WM (bottom). The RMSE for the signal without smoothing, with Gaussian smoothing, and tissue-weighted (TW) smoothing is also calculated.
	<img src="demo_RMSE_segments.png" style="zoom: 1O0%;" />
  Quite obviously, averaging the noisy signal over 20 subjects irons out the signal over the middle part of the segments. Specifically
  
    - Without smoothing (blue line), the signal deviates on the segment extremities, 1 or 2 voxels at most, because of the anatomical variance introduced. 
    - When applying Gaussian smoothing (red line), the signal is a bit smoother in the middle part of each tissue segment but deviates largely close to the edges. This is due to the "partial volume effect", i.e. signal from different tissue classes are mixed up. WM signal is always dragged down because both adjacent GM and CSF tissues have lower intensities (100 versus 50 and 5), while the GM signal is dragged up by adjacent WM (50 vs 100) or down by adjacent CSF (50 vs5). 
  - On the contrary, the signal after tissue-weighted smoothing (green line) is smooth and remains flat **over the whole tissue segment**.
  
  The RMSE values reflects these observations with a large value for the signal without smoothing (large-ish RMSE), after Gaussian smoothing (even larger RMSE), and after tissue-weighted smoothing (tiny RMSE) . 


## Discussion

With the tissue-weighted smoothing, the averaged smoothed signal within each explicitly mask tissue segment is very close to the original signal. Some deviations are visible close to the edges, especially for the small segments w.r.t. the Gaussian kernel size, e.g. WM segments, of width 8 and 12, and GM segments, of size 6 and 12, on the right side of the profiles. 
With the standard smoothing tissue edges are completely erased and the signal mixed up between compartments.

The RMSE of the averaged noisy signal (over the 20 subjects) without smoothing is a measure of the remaining noise and the aim of the smoothing is to reduce that noise.
The RMSE using the tissue-weighted smoothing and explicit masking is very small, compared to the true signal intensity (50 and 100 for GM and WM respectively) and that of signal without smoothing. On the contrary, with the standard Gaussian smoothing, the RMSE is much larger, about 14 to 20 times larger, and increased compared to without smoothing !

## Conclusion

Tissue-weighted smoothing does make a huge difference, compared to Gaussian smoothing, and the resulting signal over the explicit mask segment is very close to the true underlying signal. 

At least this is the case in this simple simulation.

## References
- Draganski et al. (2011), "Regional specificity of MRI contrast parameter changes in normal ageing revealed by voxel-based quantification (VBQ)."
  https://doi.org/10.1016/j.neuroimage.2011.01.052
- Callaghan et al. (2014), "Widespread age-related differences in the human brain microstructure revealed by quantitative magnetic resonance imaging."
  https://doi.org/10.1016/j.neurobiolaging.2014.02.008
- Eun Lee et al. (2009), "A study of diffusion tensor imaging by tissue-specific, smoothing-compensated voxel-based analysis."
  https://doi.org/10.1016/j.neuroimage.2008.09.041
