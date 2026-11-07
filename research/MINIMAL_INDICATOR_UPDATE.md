# Minimal GSCA predictive-evaluation update

This branch starts from the original GSCA.Reg_Prime main commit and changes only predictive evaluation:

1. Fit an auxiliary regression of each focal component on an intercept and the **other component predictors of the same equation**, using in-bag component scores. Apply it to OOB scores in place of mean/zero substitution.
2. For a structural path, evaluate the equally averaged squared errors of the target component's **observed indicators**, using parent scores �� B �� C to predict them. Undo in-bag standardization to retain fixed raw indicator outcomes.

The original ALS function, lambda selection and CV criterion, bootstrap draws, parameter estimates and parameter summaries, per-bootstrap weighting, and percentile calculation are retained. The new helper performs no GSCA refit or resampling. Computational cost remains one GSCA fit per bootstrap draw, plus small auxiliary regressions per edge. No variance correction, extra bootstrap layer, or new inference rule is introduced.

## API

Call GSCA_Reg with its existing arguments. GSCA_Indicator_OOB.m must be on the same path.

- TABLE.B columns 1?4 retain the original coefficient estimate, SE and two-sided coefficient CI.
- TABLE.B columns 5?6 now report the original-style one-sided percentile lower bound and fraction of nonpositive draws for the indicator-block loss difference.
- TABLE.C retains the original standardized-indicator reconstruction loss; its replacement uses the same auxiliary-regression rule.
- Original INI.GoF remains the original fit/OPE diagnostic; it has not been redefined.
- New INI.indicatorOPE_byComponent reports mean OOB raw-indicator loss by structural target.
- ETC.indicatorOPE_Boot stores those target-block losses for each bootstrap.
- ETC.delOPE_b_Boot and ETC.delOPE_c_Boot store the updated predictive differences.
- ETC.delOPE_b_legacy_Boot and ETC.delOPE_c_legacy_Boot retain the original evaluation on exactly the same fitted bootstrap models, for an auditable comparison.

The B legacy metric concerns a component, whereas the new B metric concerns raw indicators. Their numerical magnitudes are different estimands/units and should not be compared as estimator bias. The new indicator losses are averaged **after squaring**, not obtained by squaring a mean residual.

## Scope and units

Current indicator prediction supports simple disjoint measurement blocks. The observed target block must not enter the parent scores used to predict itself. Parent blocks are assumed observed at prediction time. Raw outcome units stay fixed across bootstrap fits; use substantively comparable indicator units if interpreting their average.

With no other predictors in an equation, the auxiliary regression contains only an intercept and reduces to the in-bag component mean. The GSCA coefficient is not set to zero, and the outcome equation is not refitted.

The original quantile rank round(alpha*(B-1))+1 is deliberately preserved. At a finite bootstrap count it need not have exact decision equivalence to a different empirical-CDF p-value cutoff; this update does not silently substitute a new convention.

## Verification

test_gsca_original_flow compares the original and updated functions under identical serial RNG seeds, including CV with multiple lambda values. It checks exact equality of fitted W/C/B, CV selection, original GoF, every bootstrap parameter matrix, parameter-summary columns, and final RNG state.

test_gsca_minimal_update checks auxiliary-regression arithmetic, OOB target-input exclusion, raw indicator prediction, measurement-loss preservation, sign invariance and the singleton OLS/FWL identity.

The small simulation uses N=80/200, J_p=3, lambda=0/.1, focal coefficients 0/.1/.3, 100 datasets per cell and 199 bootstrap draws. It reports empirical distributions and rejection rates under those coefficient conditions, without calling coefficient zero the predictive-null boundary. The existing ACSI replica is also tested using the original model masks.

These checks evaluate the two requested changes. They do not claim that the retained bootstrap percentile procedure is calibrated for the finite-N predictive target.

