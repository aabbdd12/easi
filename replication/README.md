# Replication files for *Estimating the Exact Affine Stone Index demand system: the easi Stata module*

Every table, figure and number reported in the technical note
(`../docs/easi_technical_note.pdf`) is produced by the scripts in this
folder, from the two data sets shipped in `../examples/` and the module in
`../src/`. Scripts are named after what they reproduce: `Table1.do` prints
Table 1 in the note's layout, `Table8_9_10_Figure1.do` produces the results
section, `Section4-4_analytic_jacobian.do` the numbers quoted in Section 4.4.
No file outside the repository is needed, except R for the optional scripts
in `R_reference/`.

## How to run

```stata
cd <path-to-repository>/replication
do master.do                 // everything, about an hour (four bootstraps of 400 replications)
global BOOT 0
do master.do                 // everything except the bootstraps, a few minutes
do Table1.do                 // any single script
```

Every script locates the module, the data and the R reference relative to
the current directory and stops with a message if it is not run from
`replication/`; nothing needs to be edited. Logs, CSV files and figures go
to `out/`. The LaTeX tables of Section 7.7 are then built with

```
python make_tables.py out
```

Requirements: Stata 14.2 or later (17 or later for the `collect` tables;
`matlist` is used below), Python 3 for `make_tables.py`.

## What produces what

| result in the note | script | output in `out/` |
|---|---|---|
| **Table 1** — invariance to the dropped good, one-step vs iterated Σ; the dropped good recovered vs estimated | `Table1.do` | `step10.log` (the table is printed in the note's layout at the end) |
| **Table 2** — four coefficient covariances against the bootstrap, 400 replications, SRS | `Table2.do` | `step8.log` |
| **Table 3** — `vce(robust)`, `vce(cluster)`, `vce(svy)` against the design bootstrap; also the *survey design* row of Table 4 | `Table3.do` | `step15.log` |
| **Table 4** — the ratio of two estimators: rows *SRS* and *pweight* | `Table4.do` | `step14.log` |
| **Table 5** — elasticity standard errors against the bootstrap | `Table5.do` | `step9.log`, `elast_se_compare.dta` |
| **Table 6** — the dummy set perturbed away from exact collinearity; and the four `easidiag` reports of Section 6.1 | `Table6.do` | `easidiag_cases.log` |
| **Table 7** — the reading grid of Section 6.1 | (text, no script) | — |
| **Tables 8, 9, 10 and Figure 1** — the reference specification on hixdata: the header and Table 02 as printed, the `e()` matrices behind the three tables, the Engel curves | `Table8_9_10_Figure1.do` then `make_tables.py` | `hixdata_output.log`, `elast_*.csv`, `tab_*.tex`, `engel_atmeans.pdf`, `engel_asobserved.pdf` |
| **Table 11** — the audit summary of Appendix A | (text, no script) | — |
| **Table 12** — the finite-difference benchmark of the elasticity formulas; also the Slutsky check and the range of Φ (Appendix A.5–A.6) | `Table12.do` | `audit.log` |
| Section 4.3 — orientation of the price tables, `compensated_q` = Γ, Engel and Cournot aggregation | `Section4-3_orientation_aggregation.do` | `step12.log` |
| Section 4.4 — the analytic Jacobian against forward differences, six specifications, four families; the demographic standard errors cell by cell | `Section4-4_analytic_jacobian.do` | `step16.log` |
| Section 5.3 — the generated-regressor Jacobian against a numerical Jacobian of an independently written moment vector (1e-9), and its effect | `Section5-3_generated_regressor.do` | `step11.log` |
| Section 6 — `easidiag` estimates nothing, detects and names a near-dependency, predicts the iterations | `Section6_easidiag_checks.do` | `step13.log` |
| Section 7.6 — running times | `Section7-6_timing.do` | (screen) |
| Appendix A.1 — the reproduction lock: `compat` = R package to 1.8e-11 (coefficients) and 1.3e-12 (standard errors) | `SectionA1_compat_coefficients.do`, `SectionA1_compat_elasticities.do` | `step1.log`, `step2.log` |
| Appendix A — the effect of each departure from the package, the `legacy()` switches one at a time | `SectionA_corrections_one_at_a_time.do` | `step3.log` |
| Appendix A.5 — degree-zero homogeneity (the `FAIL` printed there is the *package's* formulas failing the test, as the appendix states) | `SectionA5_homogeneity.do` | `homogeneity.log` |

Most scripts are the module's test suite under these names: each asserts what
it verifies and stops on any departure, so a clean run of `master.do` is
itself the statement that the note's numbers reproduce.

## The R reference (`R_reference/`)

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
