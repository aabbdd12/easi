*! Table7.do -- Table 7 of the technical note: the influence function by
*! brute force.  For each household h of a sample of the perfect DGP
*! (dgp_easi.do) the WHOLE procedure is re-run with its weight at 1 + eps and
*! 1 - eps (iweights, not normalised): the mean shares of the instrument, the
*! first pass and the rebuilt instrument, the first stage, Sigma, the
*! iteration on y, the final refit, the elasticities of the households and of
*! the market.  U_h = d theta / d w_h by central differences, and
*! sqrt(n/(n-1) sum_h U_h^2) is the robust standard error in which no nested
*! estimator is neglected.  The reported robust standard errors are divided
*! by it.  Deterministic: no Monte Carlo.  Two specifications, 1,000
*! households each: power 2 without interaction, power 3 with py zy pz.
*! About two minutes (4,000 estimations).

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table7.do
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
capture log close t7
log using "`OUT'/bruteforce.log", replace text name(t7)

mata:
real rowvector _ni_vec(string scalar sfx, string scalar se)
{
	return((st_matrix("e(elast_exp" + sfx + se + ")"),
		vec(st_matrix("e(elast_price_nc" + sfx + se + ")"))',
		vec(st_matrix("e(elast_price_c" + sfx + se + ")"))',
		vec(st_matrix("e(elast_demo" + sfx + se + ")"))'))
}
real scalar _ni_med(real colvector x)
{
	real colvector v
	real scalar m
	v = sort(x, 1); m = rows(v)
	return(mod(m, 2) ? v[(m + 1) / 2] : (v[m / 2] + v[m / 2 + 1]) / 2)
}
// reported SE / brute-force SE, by block
void _ni_report(real rowvector SEE, real matrix U, real rowvector SEB,
	real matrix UB, real scalar n)
{
	real rowvector rat, ratb, r
	real scalar ty, bl, J, T, nb, off, len, lo, hi
	string rowvector bn, tn
	rat  = SEE :/ sqrt((n / (n - 1)) :* colsum(U:^2))
	ratb = SEB :/ sqrt((n / (n - 1)) :* colsum(UB:^2))
	printf("{txt}  %-26s %9s %9s %9s\n", "", "median", "min", "max")
	printf("{txt}  %-26s {res}%9.4f %9.4f %9.4f\n", "coefficients", _ni_med(ratb'), min(ratb), max(ratb))
	bn = ("expenditure", "uncompensated", "compensated", "demographic")
	tn = ("households", "market")
	J = 3; T = 2; nb = J + 2 * J * J + T * J
	for (ty = 1; ty <= 2; ty++) {
		off = (ty - 1) * nb
		for (bl = 1; bl <= 4; bl++) {
			len = (bl == 1 ? J : (bl <= 3 ? J * J : T * J))
			lo = off + 1 + (bl >= 2) * J + (bl >= 3) * J * J + (bl >= 4) * J * J
			hi = lo + len - 1
			r = rat[| lo \ hi |]
			printf("{txt}  %-26s {res}%9.4f %9.4f %9.4f\n", tn[ty] + ", " + bn[bl],
				_ni_med(r'), min(r), max(r))
		}
	}
	printf("{txt}  %-26s {res}%9.4f %9.4f %9.4f\n", "all elasticities", _ni_med(rat'), min(rat), max(rat))
}
end

capture program drop _bruteforce
program define _bruteforce
	args n power inter seed
	easi_dgp, n(`n') seed(`seed') power(`power')
	quietly gen double w = 1
	local SPEC "w1 w2 w3 [iw = w], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(`power') `inter' nolog notable tolerance(1e-14) iterate(1000)"
	quietly easi `SPEC' compensated
	di as txt _n "n = " e(N) ", power(`power') `inter', K = " colsof(e(b))
	mata: SEE = (_ni_vec("", "_se"), _ni_vec("_mkt", "_se")); SEB = sqrt(diagonal(st_matrix("e(V)")))'
	mata: U = J(`n', cols(SEE), .); UB = J(`n', cols(SEB), .)
	local eps 1e-3
	forvalues h = 1/`n' {
		quietly replace w = 1 + `eps' in `h'
		quietly easi `SPEC' noelastse
		mata: TP = (_ni_vec("", ""), _ni_vec("_mkt", "")); BP = st_matrix("e(b)")
		quietly replace w = 1 - `eps' in `h'
		quietly easi `SPEC' noelastse
		mata: U[`h', ] = ((_ni_vec("", ""), _ni_vec("_mkt", "")) - TP) :/ (-2 * `eps')
		mata: UB[`h', ] = (st_matrix("e(b)") - BP) :/ (-2 * `eps')
		quietly replace w = 1 in `h'
	}
	di as txt "reported SE / brute-force SE, sqrt(n/(n-1) sum_h U_h^2)"
	mata: _ni_report(SEE, U, SEB, UB, `n')
end

_bruteforce 1000 2 "" 21
_bruteforce 1000 3 "py zy pz" 21
log close t7
di as res _n "Table7: done"
