*! Table9_bootstrap.do -- Table 9 of the technical note, second column:
*! the analytic standard errors under selection against a bootstrap of the
*! whole procedure, on one sample of 10,000 households of the censored
*! perfect DGP (dgp_easi.do, censor: an exclusion restriction per good,
*! delta = 0.04).  SRS, no weight, vce(robust); 500 draws of households, the
*! probits and the system re-estimated from scratch on each.  For each block,
*! the ratio SE_analytic / SE_bootstrap (median, min, max over the cells) and
*! the cells outside the band 1.96 sqrt((k-1)/(4B)) of the bootstrap noise
*! (k = kurtosis of the draws).  About fifteen minutes.

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table9_bootstrap.do
* Nothing needs to be edited.  The check below stops with a clear message
* when the working directory is not replication/.
capture confirm file "master.do"
if _rc {
	di as error "run this script from the replication/ directory:"
	di as error "    cd <path-to-repository>/replication"
	exit 601
}
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"
local OUT "`ROOT'/replication/out"
run "`ROOT'/replication/dgp_easi.do"
capture log close t9b
log using "`OUT'/selection_bootstrap.log", replace text name(t9b)

local B 500
easi_dgp, n(10000) seed(1) power(2) censor
di as txt "buyers with a nonpositive latent share: " r(nneg)
local SPEC "w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(2) selvars(w1: q1 ; w2: q2) vce(robust) nolog"

mata:
// the elasticities of one estimation, stacked: households then market
real rowvector _el_vec(string scalar se)
{
	real rowvector v
	real scalar i
	string rowvector sf
	sf = ("", "_mkt")
	v = J(1, 0, .)
	for (i = 1; i <= 2; i++) {
		v = v, st_matrix("e(elast_exp" + sf[i] + se + ")"),
			vec(st_matrix("e(elast_price_nc" + sf[i] + se + ")"))',
			vec(st_matrix("e(elast_price_c" + sf[i] + se + ")"))',
			vec(st_matrix("e(elast_demo" + sf[i] + se + ")"))'
	}
	return(v)
}
real scalar _med(real colvector x)
{
	real colvector y
	real scalar m
	y = sort(x, 1); m = rows(y)
	return(mod(m, 2) ? y[(m + 1) / 2] : (y[m / 2] + y[m / 2 + 1]) / 2)
}
// one line of a report: median, min, max of the ratios r
void _line(string scalar nm, real rowvector r)
{
	printf("{txt}  %-30s {res}%8.4f %8.4f %8.4f %5.0f\n", nm, _med(r'), min(r), max(r), cols(r))
}
// the elasticity blocks of a vector of ratios laid out as _el_vec
void _el_lines(real rowvector rat)
{
	real scalar ty, bl, off, len
	string rowvector bn, tn
	bn = ("expenditure", "uncompensated", "compensated", "demographic")
	tn = ("households", "market")
	off = 0
	for (ty = 1; ty <= 2; ty++) {
		for (bl = 1; bl <= 4; bl++) {
			len = (bl == 1 ? 3 : (bl <= 3 ? 9 : 6))
			_line(tn[ty] + ", " + bn[bl], rat[| off + 1 \ off + len |])
			off = off + len
		}
	}
}
// e(b), the probit coefficients, the elasticities
real rowvector _all(real scalar se)
{
	real rowvector a
	if (se) return((sqrt(diagonal(st_matrix("e(V_sel)")))', _el_vec("_se")))
	a = vec(st_matrix("e(sel_alpha)"))'
	return((st_matrix("e(b)"), select(a, a :< .), _el_vec("")))
}
end

quietly easi `SPEC' compensated notable
local cn : colfullnames e(b)
mata: SEA = _all(1); BT = J(0, cols(SEA), .)
mata: st_numscalar("__nal", cols(st_matrix("e(V_sel)")) - cols(st_matrix("e(b)")))
set seed 7920
local nfail 0
forvalues b = 1/`B' {
	preserve
	bsample
	capture quietly easi `SPEC' noelastse notable
	if _rc | e(converged) != 1 local ++nfail
	else mata: BT = BT \ _all(0)
	restore
}
di as txt "`B' draws, `nfail' failed"

mata:
nm   = tokens(st_local("cn")); K = cols(nm); nal = st_numscalar("__nal")
sd   = sqrt(diagonal(quadvariance(BT)))'
ku   = J(1, cols(BT), .)
for (i = 1; i <= cols(BT); i++) {
	ku[i] = mean((BT[, i] :- mean(BT[, i])):^4) / (mean((BT[, i] :- mean(BT[, i])):^2)^2)
}
rat  = SEA :/ sd
band = 1.96 :* sqrt((ku :- 1) :/ (4 * rows(BT)))
dl   = selectindex(strpos(nm, ":delta") :> 0)
sb   = J(1, K, 1); sb[dl] = J(1, cols(dl), 0)
printf("\n{txt}Table 9, bootstrap: SE analytic / SE bootstrap, B = %f\n", rows(BT))
printf("{txt}  %-30s %8s %8s %8s %5s\n", "", "median", "min", "max", "cells")
_line("latent coefficients", rat[selectindex(sb)])
_line("delta", rat[dl])
_line("probit coefficients", rat[| K + 1 \ K + nal |])
_el_lines(rat[| K + nal + 1 \ cols(rat) |])
printf("{txt}  all cells: median %6.3f, mean |log ratio| %6.3f, outside the band %f of %f (median band %5.3f)\n",
	_med(rat'), mean(abs(ln(rat))'), sum(abs(rat :- 1) :> band), cols(rat), _med(band'))
end
log close t9b
di as res _n "Table9_bootstrap: done"
