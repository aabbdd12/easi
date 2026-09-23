# Estimating the Exact Affine Stone Index demand system: the easi Stata module

`easi` estimates the EASI demand system of Lewbel and Pendakur (2009) entirely
within Stata and Mata. No R, no external dependency.

* iterated restricted 3SLS with Slutsky symmetry, residual covariance iterated
  to convergence (the estimator is invariant to the good dropped);
* expenditure, price (compensated and uncompensated) and demographic
  elasticities, with delta-method standard errors and an analytic Jacobian;
* sampling weights (`aw`, `fw`, `pw`, `iw`); `vce(robust)` by default,
  `vce(cluster)`, and `vce(svy)`, which reads the design declared by `svyset`
  (weight, PSU, strata, finite-population correction);
* the generated-regressor term of the implicit utility index in the variance,
  and the mean budget share in the denominator of an elasticity treated as the
  estimator it is;
* `predict` (shares, residuals, `y`), `estat engel`, standard `e()` results;
* `easidiag`, a pre-estimation diagnostic that finds the causes of slow or
  failed convergence — near-collinearity, near-zero variance, sparse
  categories — without running the estimator;
* a dialog box (`db easi`) for both estimation and diagnostic;
* `compat`, which reproduces the R package `easi` 0.21 to the eleventh
  decimal so that earlier results remain reproducible; `sr_easi` is kept as
  a wrapper.

Every variance is validated against a bootstrap of the whole procedure under
simple random sampling, sampling weights and a stratified cluster design. The
technical note in `paper/` gives the model, the estimator, the elasticity
formulas and the validation.

## Installation

```stata
net install easi, from("https://raw.githubusercontent.com/aabbdd12/easi/main") replace
net get easi          // optional: the examples
```

Stata 14.2 or later. The wide elasticity tables use `collect` from Stata 17
and fall back to `matlist` below.

## Quick start

```stata
use hixdata, clear
easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,           ///
     lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers)  ///
     lnexpenditure(log_y) demographics(age hsex carown time tran)          ///
     power(5) py pz zy
matrix list e(elast_price_nc)
predict double shat*, shares
estat engel

* with a survey design
svyset psu [pw=weight], strata(strata) fpc(fpc)
easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) vce(svy)

* before a large model: what will make it slow?
easidiag w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) py pz zy
```

`help easi` and `help easidiag` document every option and stored result.

## Layout

```
src/        the command: easi.ado, easi_p.ado, easi_estat.ado, easi.sthlp,
            easi.dlg, easidiag.ado, easidiag.sthlp, sr_easi.ado, sr_easi_estat.ado
examples/     hixdata.dta (Lewbel and Pendakur's reference data),
              mex_bench.dta (a reduced Mexican survey, stratified two-stage
              design, already svyset), easi_tour.do (guided tour on hixdata),
              mex_tour.do (the survey-design features on the Mexican data),
              hixdata_for_dialog.do
paper/        the technical note (PDF)
replication/  every table and figure of the note, one script each, with a
              master.do and a README that maps them; the frozen output of the
              R package and the R scripts that produced and re-checked it
stata.toc, easi.pkg    net-install manifest
```

The module's test suite and the LaTeX sources of the note are kept outside the
repository; `replication/` carries the tests under descriptive names, so a
clean run of `replication/master.do` is the statement that the note reproduces.

## Citing

If you use `easi`, please cite the technical note,

Araar, A. (2026). *Estimating the Exact Affine Stone Index demand system:
the easi Stata module*. Zenodo. https://doi.org/10.5281/zenodo.22914059

(this version; the DOI 10.5281/zenodo.22914058 covers all versions and
resolves to the latest one — see `CITATION.cff`), and the model's paper,

Lewbel, A., and K. Pendakur. 2009. Tricks with Hicks: the EASI demand system.
*American Economic Review* 99: 827–863.

## License

GNU General Public License v3.0 or later (`LICENSE`). The frozen R reference
in `replication/R_reference/` includes the sources of the archived R package
`easi` 0.21, which is itself under the GPL.

## Author

Abdelkrim Araar, Université Laval / PEP — aabd@ecn.ulaval.ca
