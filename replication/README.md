# Replication files for *Estimating the Exact Affine Stone Index demand system: the easi Stata module*

Every table, figure and number reported in the technical note
(`../paper/easi_technical_note.pdf`) is produced by the scripts in this
folder, from the two data sets shipped in `../examples/` and the module in
`../src/`. Scripts are named after what they reproduce: `Table1.do` prints
Table 1 in the note's layout, `Table8_9_10_Figure1.do` produces the results
section, `Section4-4_analytic_jacobian.do` the numbers quoted in Section 4.4.
No file outside the repository is needed, except R for the optional scripts
in `R_reference/`.

## How to run

```stata
cd <path-to-repository>/replication
do master.do                 // everything, about two hours and a quarter
global BOOT 0
global MC 0
do master.do                 // everything without resampling, a few minutes
do Table2.do                 // any single script
```

Every script locates the module, the data and the R reference relative to
the current directory and stops with a message if it is not run from
`replication/`; nothing needs to be edited.  Logs, CSV files and figures go
to `out/`.  `global BOOT 0` skips the bootstraps (Tables 3 to 6, 9 and 10,
about forty minutes); `global MC 0` skips the brute force and the Monte
Carlo experiments (Tables 7, 8 and 9, about an hour and a half).  The LaTeX
tables of Section 8.8 (Tables 14 to 16) are then built with

```
python make_tables.py out
```

Requirements: Stata 14.2 or later, Python 3 for `make_tables.py`.
`dgp_easi.do` holds the data-generating process of Sections 5.9 and 6.8
(`easi_dgp, n() seed() [power() censor]`); the scripts that need it run it.

## What produces what

| result in the note | script | output in `out/` |
|---|---|---|
| **Table 1** — invariance to the dropped good, one-step vs iterated Σ; the dropped good recovered vs estimated | `Table1.do` | `step10.log` |
| **Table 2** — the three types of elasticities on the Mexican survey; ln of the mean against the mean of ln x; the effect of the mean of products on the compensated elasticities (Section 4.3) | `Table2.do` | `types.log` |
| Section 4.3 — the market and individuals elasticities against oracles, the total Jacobian with the y term | `Section4-3_types_checks.do` | (screen) |
| Section 4.4 — orientation of the price tables, `compensated_q` = Γ, Engel and Cournot aggregation | `Section4-4_orientation_aggregation.do` | `step12.log` |
| Section 4.5 — the analytic Jacobian against forward differences, six specifications, four families; the demographic standard errors cell by cell | `Section4-5_analytic_jacobian.do` | `step16.log` |
| **Table 3** — four coefficient covariances against the bootstrap, 500 replications, SRS | `Table3.do` | `step8.log` |
| Section 5.3 — the generated-regressor Jacobian against a numerical Jacobian of an independently written moment vector, and its effect | `Section5-3_generated_regressor.do` | `step11.log` |
| **Table 4** — `vce(robust)`, `vce(cluster)`, `vce(svy)` against the design bootstrap; also the *survey design* row of Table 5 | `Table4.do` | `step15.log` |
| **Table 5** — elasticity standard errors against the bootstrap under `pweight` | `Table5.do` | `step14.log` |
| **Table 6** — elasticity standard errors against the bootstrap, hixdata | `Table6.do` | `step9.log`, `elast_se_compare.dta` |
| **Table 7** — the influence function by brute force, two specifications of the perfect design | `Table7.do` | `bruteforce.log` |
| **Table 8** — Monte Carlo on the perfect EASI design, 700 samples | `Table8.do` | `mc_perfect.log` |
| Section 6.1 — `pimpute()` against an oracle | `Section6-1_pimpute.do` | (screen) |
| Section 6.4 — the probits against Stata's `probit`; the general engine against the 3SLS engine; the estimator against `ivregress 2sls`, `predict` and the completion of y against closed forms | `Section6-4_probits.do`, `Section6-4_general_engine.do`, `Section6-4_selection_oracle.do` | (screen) |
| Section 6.5 — `vce(bootstrap)` | `Section6-5_bootstrap_option.do` | (screen) |
| Section 6.6 — the elasticities of the expected demand against finite differences | `Section6-6_selection_elasticities.do` | (screen) |
| **Table 9** — the variance under selection: brute force on the Mexican survey; bootstrap and Monte Carlo on the censored perfect design | `Table9_bruteforce.do`, `Table9_bootstrap.do`, `Table9_montecarlo.do` | `selection_bruteforce.log`, `selection_bootstrap.log`, `selection_montecarlo.log` |
| **Table 10** — the households that do not buy on the Mexican survey: `easidiag`, the correction, the survey variance against `vce(bootstrap, svy)` | `Table10.do` | `selection_mex.log` |
| **Table 11** — the dummy set perturbed away from exact collinearity; and the four `easidiag` reports of Section 7.1 | `Table11.do` | `easidiag_cases.log` |
| **Table 12** — the reading grid of Section 7.1 | (text, no script) | — |
| Section 7 — `easidiag` estimates nothing, detects and names a near-dependency, predicts the iterations; its selection section against 1/(1 − R²) | `Section7_easidiag_checks.do`, `Section7_easidiag_selection.do` | `step13.log` |
| **Table 13** (Section 8.7) — running times | `Table13.do` | (screen) |
| **Tables 14, 15, 16 and Figure 1** — the reference specification on hixdata: the header and Table 02 as printed, the `e()` matrices behind the three tables, the Engel curves | `Table14_15_16_Figure1.do` then `make_tables.py` | `hixdata_output.log`, `elast_*.csv`, `tab_*.tex`, `engel_atmeans.pdf`, `engel_asobserved.pdf` |
| **Table 17** — the audit summary of Appendix A | (text, no script) | — |
| **Table 18** — the finite-difference benchmark of the elasticity formulas; also the Slutsky check and the range of Φ (Appendix A.5–A.6) | `Table18.do` | `audit.log` |
| Appendix A.1 — the reproduction lock: `compat` = R package to 1.8e-11 (coefficients) and 1.3e-12 (standard errors) | `SectionA1_compat_coefficients.do`, `SectionA1_compat_elasticities.do` | `step1.log`, `step2.log` |
| Appendix A — the effect of each departure from the package, the `legacy()` switches one at a time | `SectionA_corrections_one_at_a_time.do` | `step3.log` |
| Appendix A.5 — degree-zero homogeneity (the `FAIL` printed there is the *package's* formulas failing the test, as the appendix states) | `SectionA5_homogeneity.do` | `homogeneity.log` |

Most scripts are the module's test suite under these names: each asserts what
it verifies and stops on any departure, so a clean run of `master.do` is
itself the statement that the note's numbers reproduce.  The scripts of
Tables 2 and 7 to 10 are measurements rather than tests: they print the
numbers the note reports.

## The R reference (`R_reference/`)

`easi_hixdata.R` is the reference model on hixdata with the R package and
nothing else — the vignette specification, the three elasticity families,
the running time — for a reader who wants to see the package itself run (it
uses the compiled package under R 3.4.4 when `check_binary_r34.R` has
installed it, and sources the package's functions otherwise).

`out/R_*.csv` is the frozen output of the R package `easi` 0.21 on the
vignette specification, at full double precision, against which the
`SectionA1_*` scripts lock the `compat` mode. It was produced by
`make_reference.R` under R 4.3.0, sourcing the package's own functions
(`easi_R_source/`, extracted from the archived package: the compiled package
does not load under R ≥ 4). `check_binary_r34.R` re-runs the *compiled*
package, unmodified, under R 3.4.4 and compares every table with the frozen
reference (2.6e-9 on the coefficients, 1e-9 to 1e-15 on the elasticity
tables); it also shows two of the package's defects directly in its output.
Neither R script is needed to reproduce the note: the frozen CSV files are
what the Stata scripts read. The logs of both runs are kept beside them
(`make_reference.log`, `check_binary_r34.log`); they carry the R running
times quoted in Section 8.7 of the note — 45.8 s under R 4.3.0, 43 to 46 s
for the compiled package under R 3.4.4, on the same machine as the Stata
timings.

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
