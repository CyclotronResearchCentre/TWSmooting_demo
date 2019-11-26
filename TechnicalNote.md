# Tissue-weighted smoothing

## Introduction

After being warped into a common space, typically a group averaged in line with the MNI space, MR images are smoothed. This is useful to : 1/ remove some high-frequency noise, 2/ reduce the multiple comparison problem, and 3/ reduce the remaining inter-subject anatomical variability.

With quantitative MRI, one is interested in the signal coming from the GM and/or the WM, standard Gaussian smoothing would therefore mix the signal coming from these 2 tissues classes, i.e. introducing some partial volume effect! One would therefore want to smooth the images "within each tissue class". Such an approach was introduced by Drangaski et al. (2011) as a specific way to smooth data in the case of "Voxel-Based Quantification" (VBQ) analysis. 

Another technique is the "*t*issue-*s*pecific, sm*oo*thing-compe*n*sated" method, aka. T-SPOON, by Eun Lee et al.

## VBQ tissue-weighted smoothing
For each tissue class of interest, typically GM and WM, the quantitative map is smoothed according to the tissue class *posterior* probability. Tissue-weighted smoothing is thus defined as follow:

$$p_j = [g*(w_j s_j)]./[g*w_j]\:m_{\tiny \mbox{TPM}}\:m_j$$

where:

| Parameter | Meaning  |
| ------ | -------- |
| $p_j$ | Quantitative map for subject $j$ after tissue-weighted smoothing |
| $s_j$ | Participant-specific quantitative map warped to group space by deformation $\Phi_j$ |
| $\Phi_j$ | Participant-specific deformation mapping from native to group space |
| $w_j$ | Participant-specific weights given by $J_j t_j$          |
| $J_j$ | Jacobian determinants of deformation $\Phi_j$              |
| $t_j$ | Participant-specific tissue *posterior* probability map warped by deformation $\Phi_j$ |
| $g*$ | Convolution by a Gaussian smoothing kernel, i.e. Gaussian smoothing. |
| $m_{\tiny \mbox{TPM}}$ | TPM-specific mask identifying voxels with probability > 5%   |
| $m_j$ | Participant-specific mask defined as $g*w_j > 5\%$ |
| $./$ | Ratio applied voxel by voxel over the smoothed images |

The point of the 2 masks is to ensure only voxels with sufficient 

- *a priori* probability ($m_{\tiny \mbox{TPM}}$) of being of the tissue class of interest, i.e. keeping voxels that are where they are supposed to be, and 
- *a posteriori* probability ($m_j$) of being of the tissue class of interest, i.e. masking out voxels that are unlikely to be of the class of interest 

A further explicit mask should be defined at the group level for the statistical analysis.

## Explicit mask

For each tissue class of interest, i.e. GM and WM, an explicit mask is generated from the group-mean smoothed tissue probabilities.  The GM (resp. WM) mask includes voxels having, on average across all subjects, 

- \> 20% chance of being GM (resp. WM), and 
- larger probability of being GM (resp. WM) than CSF or WM (resp. GM), 

This explicit mask ensures that for a voxel will only appear in one tissue class, either GM or WM, the one with the largest probability (on average across the group)

## Simulation

To illustrate how classic and tissue-weighted smoothing affect the quantitative signal, we created synthetic data in 1D, for an easy visualization, for the 20 pseudo-subjects. Then the average across subject is displayed.

### Data generation

Here are the characteristics of the simulated data,:

- the profile contains segments of 3 tissue classes, say GM, WM, and CSF, with true intensity of 50, 100, and 5 resp.;
- some "anatomical variability" is introduced by randomly shifting the edge of all the segments of the profile by -1, 0 or +1;
- each "tissue segment" has quite distinct probabilities, i.e.  $\geq 94\%$ or $\leq5\%$, but we always have $\geq 1\%$ and the total at each voxel $=100\%$;
- the signal is constructed by summing over the 3 tissue classes the product of the tissue probability with their corresponding intensity, and adding some random noise of standard deviation 2, 2, 10 resp.;
- some randomness is also added to the tissue class profiles but we still keep $\geq 1\%$ and the total at each voxel $=100\%$.

There are therefore different sources of variability:

- anatomical, with edges being shifted left or right randomly;
- intensity, with some added noise in the signal;
- tissue class, with some added noise in the probability.

For each subject, 

- the signal is smoothed using a (standard) Gaussian kernel and using the tissue-weighted method for the GM and WM.
- the tissue probabilities are smoothed using the same Gaussian kernel

Then these signals, original and smoothed, and tissue probabilities are averaged.

### Resulting plots

1. Top left, the signal from the 20 subjects, thin lines, and the mean signal intensity bold line.
   One can see the variability in signal intensity and some "fussiness" close to the edges.
2. Top right, the signal from the 20 subjects, thin lines, and the mean signal intensity, bold line, smoothed with the Gaussian kernel.
3. Bottom Left, the signal from the 20 subjects, thin lines, and the mean signal intensity, bold line, smoothed with the tissue-weighted method.

![](TissueW_smoothing_demo.png)

## References
- Draganski et al. (2011), https://doi.org/10.1016/j.neuroimage.2011.01.052
- Callaghan et al. (20XY)
- Eun Lee et al. (2009), "A study of diffusion tensor imaging by tissue-specific, smoothing-compensated voxel-based analysis" https://doi.org/10.1016/j.neuroimage.2008.09.041