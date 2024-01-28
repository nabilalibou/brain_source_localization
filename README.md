# Brain Source Localization
University Project: Localization of brain sources measured by EEG using 3 approaches: Gibbs sampler 
(Markov chain Monte Carlo algorithm), Minimum Norm Estimates (MNE) and Source Imaging based on Structured Sparsity 
(SISSY).

Inverse problems consist of estimating the parameters of a model from incomplete and noisy measurements. 
Here, the aim is to identify the distribution of dipolar sources present in the cortex illustrated by the unknown source
vector s from measurements of the potential generated by it at the level of the scalp corresponding to observation x. 
We thus have: x = **A***s*+*n*  
With **A** the mixing matrix (called forward operator) of size NxD multiplied by the sources amplitudes *s*. 
N is the number of sensors and D the number of dipoles. 
With N << D and *n* the additive white Gaussian noise of unknown variance corrupting the signal. 
Using different simulation hypotheses, the aim is to reproduce the 3D brain with active/non-active areas as close as 
possible to the given result (finding the best model for the observed data is called the forward problem). 
To achieve this, knowing the ground truth ([Fig. 1](#fig1)), our task will be to implement 3 approaches:
- The Gibbs sampler is a Markov chain Monte Carlo (MCMC) algorithm.
- Minimum Norm Estimates (MNE) based on Tikhonov (L2 norm) regularization.
- Source Imaging based on Structured Sparsity (SISSY) based on sparse Total Variation (sparse TV) regularization or 
TV-L1 regularization.

These 3 algorithms will be tested and their performance compared using the "dipole localization error" (DLE) evaluation 
criterion.

<br/>
<a id="fig1"> </a>
<p align="center">
<img src="doc\source_original.png" width="382" height="274">
</p>

<center>
<em>
Figure 1: Visualization of the 2 active sources.
</em>
</center>

The problem is under-determined as the number of parameters to estimate is greater than the number of measurements. 
In such settings, several different source configurations can explain the experimental data and additional constraints 
are needed to provide a sound solution. To tackle this problem, we need to use a priori knowledge on the characteristics 
of realistic brain sources configuration and localization. The solutions will be regularized/penalized according to
these hypotheses:
The first assumptions we make about the data, and which will be used in the estimation of s are :
- Time correlation is discarded: *s* (19k x 200) is reduced to a particular *t* and we then consider *s* (19k x 1).
- Sparsity assumption: as shown in Fig. 1, the surface of the cortex contains a lot of blue and therefore a lot of 
inactive dipoles corresponding to a source value of 0.
- Spatial correlation hypothesis: there are 2 active zones, grouped together (this hypothesis is not exploited by the 
Gibbs sampler but will be by the 2 others).
- Using the vector *s* (19k x 1), we can map the 3D figure (surface of the brain with a localization criterion and a 
gray level criterion).

# Gibbs Sampler

The Markov chain Monte Carlo method is one of a class of sampling methods based on probability distributions. 
The idea is to generate a sequence of samples where each sample depends uniquely on the preceding sample (we would like 
to ensure that if the preceding sample is accepted, then the second sample generated is close to the one just accepted).


<br/>
<a id="fig2"> </a>
<p align="center">
<img src="doc\snr_impact_gibbs.png" width="1102" height="252">
</p>

<center>
<em>
Figure 2: Influence of SNR on the Gibbs sampler.
</em>
</center>


<br/>
<a id="fig3"> </a>
<p align="center">
<img src="doc\MNE.png" width="618" height="540">
</p>

<center>
<em>
Figure 3: MNE algorithm implementation.
</em>
</center>


<br/>
<a id="fig4"> </a>
<p align="center">
<img src="doc\lambda_determination.png" width="719" height="279">
</p>

<center>
<em>
Figure 4: Lambda determination.
</em>
</center>


<br/>
<a id="fig5"> </a>
<p align="center">
<img src="doc\sissy_lambda_impact.png" width="450" height="620">
</p>

<center>
<em>
Figure 5: Lambda influence.
</em>
</center>


<br/>
<a id="fig6"> </a>
<p align="center">
<img src="doc\sissy_alpha_impact.png" width="577" height="594">
</p>

<center>
<em>
Figure 6: Alpha influence.
</em>
</center>


<br/>
<a id="fig7"> </a>
<p align="center">
<img src="doc\lambda_determination2.png" width="568" height="328">
</p>

<center>
<em>
Figure 7: Lambda determination.
</em>
</center>


<br/>
<a id="fig8"> </a>
<p align="center">
<img src="doc\comparison_3_algo_with_noise.png" width="597" height="600">
</p>

<center>
<em>
Figure 8: Analysis of the 3 algorithms with Gaussian noise.
</em>
</center>


<br/>
<a id="fig9"> </a>
<p align="center">
<img src="doc\comparison_3_algo_with_corr_noise.png" width="609" height="592">
</p>

<center>
<em>
Figure 9: Analysis of the 3 algorithms with spatially correlated Gaussian noise.
</em>
</center>