# Fixed-indicator predictive contribution for GSCA

This research extension evaluates structural paths on **fixed observed target indicators**. The old component-target functions, simulations and results remain available for comparison.

## Research progression and estimand

The research program progresses from finite-sample MLR under multivariate normality, through resampling-based MLR and single-outcome ERA (predictor blocks), to GSCA (predictor and outcome blocks).

For an edge $p\to q$, with target indicators $Y_{q1},\ldots,Y_{qJ_q}$, define

$
\mathrm{EPE}^{(q)}_m(N)
=\frac{1}{J_q}E_{\mathcal D_N,Z}
\sum_{j=1}^{J_q}(Y_{qj}-\widehat f_{m,qj}(X;\mathcal D_N))^2,
\quad
\Delta\mathrm{EPE}^{(q)}_p(N)
=\mathrm{EPE}^{(q)}_{-p}(N)-\mathrm{EPE}^{(q)}_1(N).
$

The expectation includes independent training samples of size $N$ and new observations. Outcomes have a fixed identity and unit across fits; weights, coefficients and loadings remain learned. Raw-unit squared errors are averaged equally across target indicators. Use substantively comparable units or prespecified fixed evaluation scales. Training standardization is undone before evaluating losses.

## Prediction and replacement

For each target $q$, observed parent-block scores are passed through the structural coefficients to form a **predicted** target score. Target loadings then map that prediction to the observed indicators:

$
\widehat Y_{qj}
=\widehat\mu_j+
\widehat s_j\,\widehat c_{qj}
\sum_{r\in\mathrm{Pa}(q)}
\widehat b_{rq}\widehat\gamma_r(X).
$

Observed target indicators are never used to compute the target score for its own prediction. This is a local conditional prediction: all specified parent blocks are assumed observed at prediction time. It is not a simultaneous forecast with every endogenous block unavailable.

The replacement predictor is the affine linear projection of the source score onto the other specified parent scores of the same equation. The full fitted $W,B,C$ are retained; no reduced outcome refit occurs. The estimand uses the population projection conditional on fitted weights. Sample-reference calculations plug in full-sample predictor moments and are evaluated as approximations to that same target, so nuisance-estimation error is included in their bias.

Supported initial scope: simple disjoint reflective indicator blocks, weights supported within those blocks, and no structural self paths. Unsupported crossloadings or target-input overlap produce explicit errors.

## Estimators and one-sided inference

The source observation is the sampling unit. In OOB evaluation, each observation first averages its own loss differences over its available OOB fits and then receives equal total weight. For $K_i$ OOB fits:

$
\widehat\Delta_{\rm raw}
=\frac{1}{NJ_q}\sum_{i,j}\overline{e_{-p,bij}^2-e_{1,bij}^2},
\qquad
\widehat C
=\frac{1}{NJ_q}\sum_{i,j}\{s_b^2(e_{1,bij})-s_b^2(e_{-p,bij})\}.
$

The corrected candidate is $\widehat\Delta_{\rm VC}=\widehat\Delta_{\rm raw}+\widehat C$.
Sample variances use denominator $K_i-1$. Since the observed outcome is fixed across OOB fits, residual variance equals prediction variance. No moving-target covariance term remains. Squaring occurs **before** averaging over indicators; cross-indicator independence is not required.

The identity
$
\widehat\Delta_{\rm VC}
=\frac{1}{NJ_q}\sum_{i,j}(\bar e_{-p,ij}^{\,2}-\bar e_{1,ij}^{\,2})
+\frac{1}{NJ_q}\sum_{i,j}
\frac{s_b^2(e_{1,bij})-s_b^2(e_{-p,bij})}{K_i}
$
does **not** establish unbiasedness for the original finite-$N$, individually fitted GSCA target. It identifies a variance/bagging correction candidate that must be checked against independent truth. Corrections may have either sign and may worsen some conditions.

All reported predictive bounds are **one-sided 95% lower bounds** at the default alpha. Individual OOB percentile and shifted-OOB percentile are diagnostic benchmarks. Nested percentile instead resamples and recomputes the entire raw/corrected statistic, including fitting and empirical projection moments. Ordinary empirical bootstrap positions, including duplicate rows, are the sampling units within each outer sample. This procedure is not a universal calibration theorem.

Empirical p-values are fractions of draws at or below zero. The inverse-empirical-CDF lower bound rejects if positive, equivalently $p<\alpha$ with strict handling of ties. Zero reported p-values are Monte Carlo grid values, not exact zero probabilities.

Failed draws are retained in diagnostics and invalidate the estimator; no silent redraw or omission occurs. Every observation needs at least two OOB fits. Any invalid outer replicate suppresses outer inference.

## API and coordinates

    addpath('research');
    out = gsca_indicator_predictive(Data,W0,C0,B0,struct( ...
        'Boot',199,'Outer',199,'InnerBoot',99,'Lambda',0,'Seed',20261107));

- B-only edge order is MATLAB find(B0), available in out.edges.
- Main fields: deltaRaw, correction, deltaCorrected; pRaw/pShifted and their lower bounds; pOuterRaw/pOuterCorrected and their lower bounds.
- pLegacy is zero deletion under the **new indicator outcome**. It is not the previous component-target p-value and it targets a different replacement contrast.
- out.indicator contains coordinate-level estimates. out.aggregation maps coordinates to mean target-block losses.
- gsca_indicator_maps constructs affine residual maps per edge-indicator coordinate.
- gsca_exact_risk integrates their losses exactly from population moments; multiply its coordinate risks by maps.aggregation.
- ReferenceMean/ReferenceCov may specify population moments. Otherwise they are estimated afresh inside each outer resample.
- Lambda is a fixed prespecified scalar on the mean-loss penalty scale. Data-adaptive tuning is outside the present implementation.

## Validation and simulations

    test_gsca_fit; test_gsca_risk_maps; test_gsca_predictive;
    test_gsca_indicator_maps; test_gsca_indicator_predictive;
    run_indicator_experiment('smoke');
    run_indicator_experiment('primary_n80_jp3');
    run_indicator_experiment('nested_n80_sample_boundary');

Independent tests reconstruct singleton OLS bootstrap estimates directly, verify coordinate aggregation, projection, absence of target-input leakage, raw-unit transformations, fixed-outcome variance equality, equal-observation weighting, nested reference recomputation, and common random streams across inner budgets.

V4 uses the previous controlled three-component DGP, with two correlated parents and one endogenous component. Indicators form exact components in the DGP. Normal and heteroscedastic designs share population covariance. Every indicator has population variance one; $J_p=1$ and $3$ distinguish the analytic bridge from learned block weights.

For $J_p=1,\lambda=0$, normal data, the exact target remains
$
\Delta\mathrm{EPE}_p(N)
=.91b^2-\frac{.8775-b^2-.21b}{N-4}.
$
Other boundaries are independently calibrated for the new indicator target. Population affine risk integration eliminates test-sample MC noise, but training-sample MC error remains and is reported. Boundary labels without an analytic formula are near-boundary candidates, not known exact nulls.

Primary jobs cover $N=80,200$, $J_p=1,3$, $\lambda=0,.1$, normal/heteroscedastic errors, coefficient zero/near-boundary/.30 and oracle/sample references. Nested jobs focus on representative and robustness conditions. Inner99/199 comparison uses identical original and outer datasets and the same prefix of inner samples. Each job saves CSV, JSON, MAT, seeds and failure diagnostics.

The old six-component manuscript simulation has not been replaced by this controlled experiment. The optional run_indicator_example function evaluates the existing ACSI replica as a practical integration check; it does not provide simulation truth or validate nested inference.

See GSCA_Indicator_Correction_Report_KO.md for completed results and limitations.

