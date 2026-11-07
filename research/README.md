# GSCA OOB correction experiment

This research module evaluates corrections for the unnormalized predictive contribution \(\Delta\mathrm{EPE}_p(N)\). It keeps the public `GSCA_Reg` interface unchanged. The methods below are experimental: a computational check is not evidence of predictive-null calibration.

## Target and replacement

For each loading or structural edge, fit the complete GSCA model once. Replace its source component by its linear projection onto the **other specified parents of the same equation**, including an intercept. Do not refit the reduced outcome equation. The nuisance projection is conditional on the fitted component weights.

Observed-indicator losses are expressed in raw indicator units. Structural losses concern the fitted endogenous component, standardized to unit variance under reference X moments. Consequently, this is an algorithm-dependent component target; it is not risk for a fixed unobserved latent variable. This choice fixes outcome units across fits and introduces no EPE denominator. An externally supplied reference covariance gives an oracle experiment; the default full-sample X covariance gives a plug-in approximation. The latter has additional nuisance estimation and dependence, and must be assessed separately.

A fixed first-indicator sign convention is used in every fit, with explicit failure near a zero anchor. Lambda is prespecified and scalar. Selecting lambda from the full data and freezing it during OOB evaluation is deliberately not used: that would let the held-out observations inform tuning. A study of tuned GSCA requires selection to be repeated inside each training sample.

## Candidate correction

Let \(e_{1,bi}\) and \(e_{-p,bi}\) be the paired full and removed residuals for observation \(i\) when it is OOB in draw \(b\). With \(K_i\) such draws,

\[
\widehat{\Delta\mathrm{EPE}}_p^{\mathrm{raw}}
=\frac1N\sum_i\overline{(e_{-p,bi}^2-e_{1,bi}^2)},\qquad
\widehat C_p=\frac1N\sum_i\{s^2(e_{1,bi})-s^2(e_{-p,bi})\},
\]
\[
\widehat{\Delta\mathrm{EPE}}_p^{\mathrm{VC}}
=\widehat{\Delta\mathrm{EPE}}_p^{\mathrm{raw}}+\widehat C_p.
\]

The sample variances use divisor \(K_i-1\). Each original observation has equal weight; `rawReplicateMean` is also returned to identify the separate effect of aggregation. **The correction may be negative.**

This construction satisfies the exact identity

\[
\widehat{\Delta\mathrm{EPE}}_p^{\mathrm{VC}}
=\frac1N\sum_i(\bar e_{-p,i}^{\,2}-\bar e_{1,i}^{\,2})
+\frac1N\sum_i\frac{s^2(e_{1,bi})-s^2(e_{-p,bi})}{K_i}.
\]

It is a bagging/finite-bootstrap identity. It is **not a theorem that this estimates the original unbagged finite-N target without bias**. It removes a variance contribution but need not remove second-order bias of learned weights or penalized coefficients. The simulation compares it with independent-training truth rather than assuming the OLS \(+C\) result transfers.

If the outcome \(t_b\) changes with the fitted weights, using prediction variances alone is incorrect:

\[
s^2(t_b-f_b)-s^2(t_b-g_b)
=s_f^2-s_g^2-2\widehat{\operatorname{Cov}}(t_b,f_b-g_b).
\]

The implementation uses residuals throughout.

## One-sided percentile candidates

All predictive bounds are one-sided lower bounds at \(\alpha=.05\). No two-sided predictive interval is substituted.

| Output | Interpretation |
|---|---|
| `pLegacy`, `lowerLegacy` | Zero-deletion OOB benchmark, on the comparable outcome scale |
| `pRaw`, `lowerRaw` | Projection replacement; empirical percentile of individual OOB replicate loss differences |
| `pShifted`, `lowerShifted` | Same replicate distribution shifted by \(\widehat C_p\); a diagnostic, not a sampling-distribution repair |
| `pOuterRaw`, `pOuterCorrected` | Outer pairs bootstrap of the entire raw or corrected aggregate; all training, projection and correction repeated |

For draws \(T_b\), \(p=B^{-1}\sum_b1(T_b\le0)\). Bounds use the inverse empirical CDF (order \(\lceil\alpha B\rceil\)). Reject when the lower bound is positive, equivalent to \(p<\alpha\), with a strict inequality at finite-B ties. An observed \(p=0\) means zero nonpositive Monte Carlo draws, not an exactly zero population tail probability.

Ordinary nested bootstrap uses outer-sample **positions**, including duplicate values, as its inner sampling units. It does not claim that repeated original source IDs are independent held-out population observations. Excluding all copies by source ID would define a different estimator. Finite inner-bootstrap noise remains and should be checked by increasing `InnerBoot`.

Nested percentile does not automatically solve first-order degeneracy. See Chen and Fang (2019), [Inference on Functionals under First Order Degeneracy](https://arxiv.org/abs/1901.04861). A future boundary-constrained or second-order bootstrap would be a separate method, not silently presented as an ordinary percentile interval.

## Running

```matlab
addpath('research');
test_gsca_fit;
test_gsca_risk_maps;
test_gsca_predictive;

opts = struct('Boot',999,'Outer',0,'Lambda',0,'Seed',20261106);
result = gsca_predictive(Data,W0,C0,B0,opts);
% Check result.valid before using results.
% result.deltaRaw, result.correction, result.deltaCorrected
% result.pRaw, result.pShifted

% Nested whole-estimator percentile (much more expensive):
opts.Outer=499; opts.InnerBoot=999;
result = gsca_predictive(Data,W0,C0,B0,opts);
% Also check result.outerValid before using outer inference.
```

A failed fit, empty OOB draw, invalid projection, or insufficient OOB coverage invalidates the point estimate. Fits are not silently redrawn or omitted. If any outer replicate fails, outer inference is unavailable. Outputs retain reasons and counts. These strict rules reveal where the specified algorithm is not reliably defined; successful-only results must not be described as unconditional operating characteristics.

## Simulation update

`Run_Simulation_Correction_v3` replaces the old simulation's obsolete `RegGSCA` multi-output calls with the explicit research API. It uses a self-contained controlled composite DGP before returning to the more complex ACSI design. It does **not** reproduce or replace the old ACSI results. It needs neither the private `DGP_SEM` path nor `gramm`.

The controlled design has two correlated predictor components and one endogenous component, with singleton indicators for the OLS bridge and three indicators per component for learned-weight GSCA. Fixed zero and positive regularization are separated. The coefficient-zero case, an estimated finite-N predictive boundary, and a positive coefficient are distinct conditions.

Each fitted residual is affine in raw X. Known population moments give its exact conditional risk; independent training sets then estimate \(\Delta\mathrm{EPE}_p(N)\). Boundary calibration uses separate random streams from truth validation and performance simulation. Calibration MCSE is reported: an approximate boundary is not an exact known null. Type I error/power labels require the sign of validated predictive truth, not just a coefficient-zero/nonzero label.

The planned comparisons include raw/corrected bias and RMSE, one-sided lower-bound coverage, rejection rates with Monte Carlo uncertainty, solver failures, and oracle versus sample-X references. Small nested runs establish executability only. Broader nested simulations, increased B, additional sample sizes, heteroscedasticity and the original ACSI design are needed before choosing a final procedure.

## Scope of solver changes

The research solver preserves the GSCA mean-loss penalty convention and uses covariance algebra equivalent to the original weight-update direction. It refreshes dependent component targets during coordinate updates, refits coefficients before reporting the objective, handles measurement-only models, records convergence, and rejects ill-conditioned required subproblems. These changes mean numerical differences from legacy output can reflect solver fixes as well as inference changes. Production files are retained for explicit comparisons.

No empirical correction coefficient is tuned to make a reported rejection rate equal .05.
