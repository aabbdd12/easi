# Replication files for *Estimating the Exact Affine Stone Index demand system: the easi Stata module*

Everything reported in the technical note (`../docs/easi_technical_note.pdf`)
is produced by the scripts in this folder, from the two data sets shipped in
`../examples/` and the module in `../src/`. No file outside the repository is
needed, except R for the optional scripts in `R_reference/`.

## How to run

```stata
cd replication
do master.do                 // everything, about an hour (four bootstraps of 400 replications)
global BOOT 0
do master.do                 // everything except the bootstraps, a few minutes
```

Each script can also be run on its own, from this directory, with `do`.
Logs, CSV files and figures go to `out/`. The LaTeX tables of Section 7.7
are then built with

```
python make_tables.py out
```

Requirements: Stata 14.2 or later (17 or later for the `collect` tables;
`matlist` is used below), the module in `../src` (no installation needed,
the scripts put it on the adopath), Python 3 for `make_tables.py`.

## What produces what

| note | script | output in `out/` |
|---|---|---|
| Appendix A.1, the reproduction lock: `compat` = R package to 1.8e-11 (coefficients) and 1.3e-12 (standard errors) | `01a_compat_lock_coefficients.do`, `01b_compat_lock_elasticities.do` | `step1.log`, `step2.log` |
| Appendix A, effect of each departure from the package (E1–E3, L1–L7), the `legacy()` switches one at a time | `02_effect_of_each_correction.do` | `step3.log` |
| Table 1, invariance to the dropped good; one-step vs iterated Σ | `03_invariance_to_dropped_good.do` | `step10.log` |
| Appendix A.5, the finite-difference benchmark of the elasticity formulas, the Slutsky check; A.6, the range of Φ | `04a_formulas_finite_differences.do` | `audit.log` |
| Appendix A.5, degree-zero homogeneity | `04b_formulas_homogeneity.do` | `homogeneity.log` |
| Table 2, the four coefficient covariances against the bootstrap (400 replications, SRS) | `05a_coefficient_se_bootstrap.do` | `step8.log` |
| Table 6, the elasticity standard errors against the bootstrap; Appendix A, item L7 | `05b_elasticity_se_bootstrap.do` | `step9.log`, `elast_se_compare.dta` |
| Section 5.2, the generated-regressor Jacobian against a numerical Jacobian of an independently written moment vector (1e-9) | `06_generated_regressor.do` | `step11.log` |
| Table 4, rows *pweight*: standard errors under sampling weights against a weighted bootstrap | `07a_pweight_bootstrap.do` | `step14.log` |
| Table 3 and Table 4, row *survey design*: `vce(robust)`, `vce(cluster)`, `vce(svy)` against the design bootstrap (PSUs within strata) | `07b_design_bootstrap.do` | `step15.log` |
| Section 4.4, the analytic Jacobian against forward differences on six specifications and four families; the demographic standard errors cell by cell | `08_analytic_jacobian.do` | `step16.log` |
| Section 4.3, orientation of the price tables, `compensated_q` = Γ, Engel and Cournot aggregation | `09_tables_orientation_aggregation.do` | `step12.log` |
| Section 6 and Table 5: no estimation, near-dependency detected and named, the 1e-9 perturbation, iterations predicted | `10_easidiag_checks.do` | `step13.log` |
| Section 7.7: the header and Table 02 as printed, the `e()` matrices behind Tables 7–10, Figure 1 | `11_hixdata_results.do` + `make_tables.py` | `hixdata_output.log`, `elast_*.csv`, `tab_*.tex`, `engel_*.pdf` |
| Section 6.1, the four `easidiag` cases | `12_easidiag_cases.do` | `easidiag_cases.log` |
| Section 7.6, running times | `13_timing.do` | (screen) |

The scripts `01`–`10` are the module's test suite under descriptive names:
each one asserts what it verifies and stops on any departure, so a clean run
of `master.do` is itself the statement that the note's numbers reproduce.

## The R reference (`R_reference/`)

`out/R_*.csv` is the frozen output of the R package `easi` 0.21 on the
vignette specification, at full double precision, against which scripts
`01a`/`01b` lock the `compat` mode. It was produced by
`make_reference.R` under R 4.3.0, sourcing the package's own functions
(`easi_R_source/`, extracted from the archived package: the compiled package
does not load under R ≥ 4). `check_binary_r34.R` re-runs the *compiled*
package, unmodified, under R 3.4.4 and compares every table with the frozen
reference (2.6e-9 on the coefficients, 1e-9 to 1e-15 on the elasticity
tables); it also shows two of the package's defects directly in its output.
Neither R script is needed to reproduce the note: the frozen CSV files are
what the Stata scripts read.

## Data

* `../examples/hixdata.dta` — the Canadian data of Lewbel and Pendakur
  (2009), as distributed with the R package; 4,847 households, 9 goods,
  prices and expenditure already in logarithms.
* `../examples/mex_bench.dta` — a reduced extract of Mexico's ENIGH 2014
  (cereals): 2,477 households, 703 PSUs, 28 strata, `svyset`; three goods
  (corn, wheat, everything else); centred logs `lx lp1 lp2 lp3` and their
  uncentred copies `lx_raw lp1_raw lp2_raw`; demographics and an exhaustive
  set of composition dummies `nocup0-nocup4` for the `easidiag` example. The
  construction (which strata and PSUs were kept, and why) is documented in
  the script that built it, kept with the test suite.
