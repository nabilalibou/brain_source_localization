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
<div align="center">Figure 1: Visualization of the 2 active sources.</div>
<br/>

The problem is under-determined as the number of parameters to estimate is greater than the number of measurements. 
In such settings, several different source configurations can explain the experimental data and additional constraints 
are needed to provide a sound solution. To tackle this problem, we need to use a priori knowledge on the characteristics 
of realistic brain sources configuration and localization (we regularize/penalize the solution according to hypotheses).  

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

This principle led to the creation of the Metropolis-Hastings algorithm, whose general idea is to sample elements 
according to a simple conditional law that will be kept with a certain probability (see the algorithm 
in the statement). The goal of this method is therefore to converge gradually towards the correct distribution. 
This is the most general among the MCMC family. A variant of this algorithm is the Gibbs sampler.

**Implementation of the algorithm for the Gibbs sampler**

The samples will asymptotically approach the target distributions by updating the error. 
At each iteration, we calculate the parameters for a "direction", update them, and thus create estimators. 
Once a certain number of iterations have been completed, these estimators will be injected back into the initial 
equation *x = As + n* and we will model the associated cortex. The result as a function of the SNR after sampling 100 
times can be seen in [Fig. 2](#fig2):

<br/>
<a id="fig2"> </a>
<p align="center">
<img src="doc\snr_impact_gibbs.png" width="1102" height="252">
</p>
<div align="center">Figure 2: Influence of SNR on the Gibbs sampler.</div>
<br/>

The Gibbs sampler is a good way to estimate the parameters of a target distribution and gives correct results and is 
relatively close to the truth, but it is very sensitive to noise and cannot be used if the noise is too dominant.

# Tikhonov regularization and TV-L1 regularization

When the inverse problem is ill-posed (in the sense of Hadamard) and A cannot be inverted because there are more 
unknowns than equations (N<<D), we need to regularize the problem. The regularization parameter, which corresponds to a 
Lagrange multiplier, allows us to manage the balance between the term associated with the model (error minimization) 
and the constraint term.  

A very popular regularization method is Tikhonov regularization of type L2 norm. Often, a term multiplying the source s 
in the L2 norm called the "Tikhonov matrix", must be carefully chosen for the problem considered. 
T corresponds to the mech function which, when applied to s, corresponds to the gradient of the source and then 
gives the spatial difference. This is done by calculating the gradient between each neighbor constituting the cortex 
mesh and will allow to include parsimony in the problem.

It should be noted that this regularization allows for an analytical solution that we will try to find using the MNE 
algorithm. The second type of regularization that we will use is the type of TV-L1 norm regularization.  
In contrast to the previous method, this one reflects the parsimonious character of the solution that we want to find, 
hence its name Sparse a priori.
In fact, a measure of the parsimony of a vector *s* is given by the L0 semi-norm. However, it has the disadvantage of 
leading to problems that are very difficult to calculate with standard optimization algorithms. 
That is why we prefer to proceed to its continuous and convex relaxation which is the L1 norm.
The disadvantage of this regularization is that the cost function is non-differentiable. This regularization then does 
not allow for a fixed solution and must therefore go through an optimization algorithm with a stopping criterion.
Here we will use a proximal optimization method (to a solution of type proxy mal) which is the ADMM (alternating 
direction method of multipliers). The algorithm allowing to implement it is the SISSY algorithm. T in the L1 norm is a 
linear operator that implements the gradient on the surface of the cortex mesh. The operator multiplied by s then gives 
a term translating the spatial difference between neighboring sources.

## Implementing the Minimum Norm Estimates algorithm (MNE)

We implement the algorithm and display ([Fig. 3](#fig3)) the results as a function of parameter and SNR variations.

<br/>
<a id="fig3"> </a>
<p align="center">
<img src="doc\MNE.png" width="618" height="540">
</p>
<div align="center">Figure 3: MNE algorithm implementation.</div>
<br/>

The higher the signal-to-noise ratio, the closer the solution is to the truth. In addition, the higher the 
regularization parameter, the more the constraint is taken into account, and therefore the source foci are more and 
more evident. However, although this algorithm gives interesting results, it is observed that when the parameter is 
large and the L2 norm constraint is taken into account, the solution foci are zones that extend to the edges. 
The regularization parameter should therefore not be too high.  
The problem with Tikhonov regularization is that the constraint translates the smooth character of the components of 
the vectors on the surface (dipole gain). The solution is then an extended patch with the amplitude decreasing at the 
edges, which is not physiologically realistic. The L2 norm constraint allows for smoothness and not parsimony, as we 
can see.

*Criteria for lambda ([Fig. 4](#fig4))*

<br/>
<a id="fig4"> </a>
<p align="center">
<img src="doc\lambda_determination.png" width="719" height="279">
</p>
<div align="center">Figure 4: Lambda determination.</div>
<br/>

**L-curve criterion**: This corresponds to the graph of the residual norm (the construction error) as a function of 
the solution norm varying with lambda. The curve then takes a L-shape where the "abrupt" parts correspond to a dominant 
reconstruction error while the "flat" parts correspond to a dominant regularization error. The optimal value of lambda is 
then located at the fold of the L and corresponds to approximately 0.1.

**Discrepancy principle**: This criterion exploits prior information on the quality of the measurements: the 
reconstruction error is due to noise. This heuristic will allow us to find the good lambda that satisfies the prior. 
The value of lambda such that the power of the reconstruction error is approximately equal to that 
of the noise is 300. The value of lambda is too high, probably due to an error on our part rather than the prior which would 
be incorrect.

**Generalized cross-validation**: In order to choose a good value of lambda, this heuristic is an analogy to cross-validation in 
classification and starts from the principle that this good should also allow to correctly predict the observations 
excluded in the calculation x = As + n. The regularization parameter is then chosen such that it minimizes the function 
characterizing the generalized prediction error and is equal to 30.
This is the result closest to the value of the regularization parameter deduced previously.

## Implementing the Source Imaging based on Structured Sparsity (SISSY) algorithm

The second regularization approach for the optimization problem uses the TV-L1 norm. 
This method is based on dual ascent, which consists of simplifying the primary optimization problem by several 
sub-optimization problems for which the objective is to have analytical solutions. 
Treating each term separately and iteratively has given the name alternating method of ADMM.
We then need to add terms translating these constraints and integrate them into the cost function as a penalty term. 
In order to take these constraints into account, we add two other Lagrange multipliers
(we minimize the difference between the terms). Finally, in order to increase the robustness of the method and 
accelerate its convergence, we introduce the multipliers of the augmented Lagrangians. 

**Lambda influence [Fig. 5](#fig5)**

<br/>
<a id="fig5"> </a>
<p align="center">
<img src="doc\sissy_lambda_impact.png" width="450" height="620">
</p>
<div align="center">Figure 5: Lambda influence.</div>
<br/>

Unlike the first regularization, the solution is more sparse and the source is more concentrated. 
However, it dissolves over time as the regularization parameter increases after a certain value.

**Alpha influence [Fig. 6](#fig6) (lambda set at 10)**

<br/>
<a id="fig6"> </a>
<p align="center">
<img src="doc\sissy_alpha_impact.png" width="577" height="594">
</p>
<div align="center">Figure 6: Alpha influence.</div>
<br/>

We see that regularization with the L1 norm reflects the parsimonious character of the solution well. 
The higher alpha, the smaller the source area becomes. It seems that the fact of not applying any operator to s in the L1 
norm (unlike the L2 norm where we apply the operator T) leads to favoring solutions whose norms are small and amounts 
to looking for an epileptic point rather than a focus (which is not physiologically accurate). 
Despite this hypothesis, we can say that in our case, the regularization of type TV-L1 norm gives results closer to 
reality than the regularization of type L2 norm. The good choices of parameters are therefore lambda = 10 and alpha = 0.1.

<br/>
<a id="fig7"> </a>
<p align="center">
<img src="doc\lambda_determination2.png" width="568" height="328">
</p>
<div align="center">Figure 7: Lambda determination.</div>
<br/>

[Fig. 7](#fig7) shows the measure of the sparsity of a vector is given by the L0 semi-norm. This norm is not taken into account in 
regularization because it leads to optimization problems that are very costly in terms of computational time. 
This is why a "convex relaxation" is preferred using the l1 norm (although the solution is no longer analytical). 
However, this heuristic restores the fact that the constraint that is actually to be imposed is based on the L0 norm 
and not l1. Thus, the appropriate regularization parameter is the one for which the L0-type constraint is minimal.

# Performance analysis and comparison of algorithms studied

In order to compare our 3 algorithms, Gibbs sampler, MNE and SISSY in detail, we will use the "dipole localization 
error (DLE)" criterion. This evaluation criterion allows us to compare the active dipoles of the original source 
configuration (the ground truth) with the active dipoles of the inverse solution estimated by our algorithms. 
In order to obtain a performance of the algorithms that does not depend on the noise, we will start by calculating the 
DLE for a certain number of noise realizations. By repeating several times with different noise vectors drawn randomly 
from a given distribution, we will be able to evaluate the DLE using the estimated statistical parameters (mean, 
variance, etc.). The comparison of the algorithms will be done qualitatively and then quantitatively.

## Qualitative comparison

<br/>
<a id="fig8"> </a>
<p align="center">
<img src="doc\comparison_3_algo_with_noise.png" width="597" height="600">
</p>
<div align="center">Figure 8: Analysis of the 3 algorithms with Gaussian noise.</div>
<br/>

[Fig. 8](#fig8) we can observe that the Gibbs sampler is very sensitive to noise. In fact, its convergence speed tends 
to infinity as the SNR decreases. In conclusion, if the SNR is not sufficiently high, the results of Gibbs are not usable.

MNE displays solutions that are too extended and sources that are too dispersed. (L2 constraint). Despite this, its 
performance is more acceptable than for Gibbs additionally it is much more stable with respect to noise. 
Its performance is not excellent but its robustness is interesting.

Finally, Sissy has without a doubt the best performance when the Gaussian noise is low. In addition, although its 
performance is very sensitive to noise, the consideration of the parsimony of this type of regularization means that 
even with a low SNR, the results are more physiologically realistic.


<br/>
<a id="fig9"> </a>
<p align="center">
<img src="doc\comparison_3_algo_with_corr_noise.png" width="609" height="592">
</p>
<div align="center">Figure 9: Analysis of the 3 algorithms with spatially correlated Gaussian noise.</div>
<br/>

[Fig. 9](#fig9) corresponds to a type of noise models Gaussian activity of dipoles that do not correspond to the sources, 
called background activity. Gibbs is less sensitive to this type of noise, it converges more easily to correct estimators. 
The results are more accurate than for Gaussian noise at the sensor level. Again, these results are more accurate but 
are only usable if the SNR is high as this method is not very robust. MNE shows better performance than Gibbs. 
The results are constant with respect to the SNR but are limited by this fact. This method is very robust but the
corresponding results are less close to reality than for the Sissy algorithm. 
Sissy is more robust when the noise is background activity at the dipole level. Its performance is much higher than 
that of the other two algorithms and is very close to reality from a SNR of 1.