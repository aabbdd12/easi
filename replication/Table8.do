*! Table8.do -- Table 8 of the technical note: Monte Carlo of the analytic
*! standard errors of the elasticities (households and market) on the
*! perfect EASI DGP (dgp_easi.do): 700 samples of 3,000 households (seeds
*! 5001 ... 5700), SRS, no weight, vce(robust).  For each block: the mean
*! analytic SE divided by the SD of the estimates across the samples (median,
*! min, max over the cells), the cells outside the band 1.96 sqrt((k-1)/(4R))
*! (k = kurtosis of the estimates), the bias against the population values
*! (the same estimator on 1,000,000 households of the DGP) in units of the
*! SD, and the rejection rate of the Wald test at 5% of H0: theta = population
*! value.  About ten minutes.

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table8.do
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
capture log close t8
log using "`OUT'/mc_perfect.log", replace text name(t8)

local R 700
local seed0 5000
local n 3000
local power 2
local SPEC "w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(`power') nolog"

mata:
real rowvector _bt_vec(string scalar sfx, string scalar se)
{
	return((st_matrix("e(elast_exp" + sfx + se + ")"),
		vec(st_matrix("e(elast_price_nc" + sfx + se + ")"))',
		vec(st_matrix("e(elast_price_c" + sfx + se + ")"))',
		vec(st_matrix("e(elast_demo" + sfx + se + ")"))'))
}
E = J(0, 54, .); S = J(0, 54, .)
end

* the samples
local nfail 0
forvalues r = 1/`R' {
	quietly easi_dgp, n(`n') seed(`=`seed0' + `r'') power(`power')
	capture quietly easi `SPEC' compensated notable
	if _rc local ++nfail
	else {
		mata: E = E \ (_bt_vec("", ""), _bt_vec("_mkt", ""))
		mata: S = S \ (_bt_vec("", "_se"), _bt_vec("_mkt", "_se"))
	}
}
di as txt "`R' samples of `n' households, `nfail' failure(s)"

* the population values
quietly easi_dgp, n(1000000) seed(987654321) power(`power')
quietly easi `SPEC' noelastse notable

mata:
real scalar _med(real colvector x)
{
	real colvector y
	real scalar m
	y = sort(x, 1); m = rows(y)
	return(mod(m, 2) ? y[(m + 1) / 2] : (y[m / 2] + y[m / 2 + 1]) / 2)
}
real scalar _q(real colvector x, real scalar p)
{
	real colvector y
	y = sort(x, 1)
	return(y[max((1, ceil(p * rows(y))))])
}
POP = (_bt_vec("", ""), _bt_vec("_mkt", ""))
R  = rows(E)
sd = sqrt(diagonal(quadvariance(E)))'
ms = mean(S)
ku = J(1, 54, .)
for (i = 1; i <= 54; i++) {
	ku[i] = mean((E[, i] :- mean(E[, i])):^4) / (mean((E[, i] :- mean(E[, i])):^2)^2)
}
rat  = ms :/ sd
band = 1.96 :* sqrt((ku :- 1) :/ (4 * R))
bias = (mean(E) - POP) :/ sd
rej  = mean((abs(E :- POP) :/ S) :> invnormal(0.975))
bb   = 1.96 * sqrt(0.05 * 0.95 / R)
bn = ("expenditure", "uncompensated", "compensated", "demographic")
tn = ("households", "market")
J = 3; T = 2; nb = J + 2 * J * J + T * J
printf("\n{txt}Table 8. Monte Carlo, R = %f samples: mean analytic SE / SD of the estimates\n", R)
printf("{txt}%-11s %-14s %7s %7s %7s %5s %8s %9s %8s\n", "type", "block", "median", "min", "max",
	"cells", "outside", "max|b|/SD", "Wald %")
for (ty = 1; ty <= 2; ty++) {
	off = (ty - 1) * nb
	for (bl = 1; bl <= 4; bl++) {
		len = (bl == 1 ? J : (bl <= 3 ? J * J : T * J))
		lo = off + 1 + (bl >= 2) * J + (bl >= 3) * J * J + (bl >= 4) * J * J
		hi = lo + len - 1
		r = rat[| lo \ hi |]
		printf("{txt}%-11s %-14s {res}%7.3f %7.3f %7.3f %5.0f %8.0f %9.3f %8.2f\n", tn[ty], bn[bl],
			_med(r'), min(r), max(r), len, sum(abs(r :- 1) :> band[| lo \ hi |]),
			max(abs(bias[| lo \ hi |])), 100 * mean(rej[| lo \ hi |]'))
	}
}
printf("{txt}(band 1.96 sqrt((k-1)/(4R)): median %5.3f; kurtosis median %5.2f, max %5.2f)\n",
	_med(band'), _med(ku'), max(ku))
printf("{txt}(Wald: all cells %5.2f%%, one cell 5%% +/- %4.2f, cells outside %f;\n",
	100 * mean(rej'), 100 * bb, sum(abs(rej :- 0.05) :> bb))
printf("{txt} quantiles of t pooled: 2.5%% %6.3f, 97.5%% %6.3f)\n",
	_q(vec((E :- POP) :/ S), 0.025), _q(vec((E :- POP) :/ S), 0.975))
end
log close t8
di as res _n "Table8: done"
