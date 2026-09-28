*! Table9_bruteforce.do -- Table 9 of the technical note, first column:
*! the influence function under selection by brute force.  A random
*! subsample of 800 households of the reduced Mexican survey; the correction
*! for the non-buyers of corn and wheat (probits on ln p - ln x and the
*! demographics, no variable of the probit only: the weakly identified case).
*! For each household the WHOLE procedure is re-run with its weight at
*! 1 + eps and 1 - eps: the probits, the completed implicit utility, the
*! instruments, Sigma, the iteration, the elasticities.  The reported robust
*! standard errors (factor n/(n-1)) are divided by sqrt(n/(n-1) sum_h U_h^2),
*! U_h = d theta / d w_h.  Deterministic; about five minutes (1,600 estimations).

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table9_bruteforce.do
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
capture log close t9a
log using "`OUT'/selection_bruteforce.log", replace text name(t9a)

local n 800
use "`ROOT'/examples/mex_bench.dta", clear
set seed 7
quietly gen double __u = runiform()
sort __u
quietly keep in 1/`n'
drop __u
quietly gen double w = 1
local SPEC "w1 w2 w3 [iw = w], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) selection nolog notable tolerance(1e-14) iterate(1000)"

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
end

quietly easi `SPEC' compensated
di as txt "N = " e(N) ", " e(iter) " iterations, K = " colsof(e(b))
local cn : colfullnames e(b)
mata: SEB = sqrt(diagonal(st_matrix("e(V)")))'; UB = J(`n', cols(SEB), .)
mata: SEE = _el_vec("_se"); UE = J(`n', cols(SEE), .)
local eps 1e-3
forvalues h = 1/`n' {
	quietly replace w = 1 + `eps' in `h'
	quietly easi `SPEC' noelastse compensated
	mata: BP = st_matrix("e(b)"); EP = _el_vec("")
	quietly replace w = 1 - `eps' in `h'
	quietly easi `SPEC' noelastse compensated
	mata: UB[`h', ] = (st_matrix("e(b)") - BP) :/ (-2 * `eps')
	mata: UE[`h', ] = (_el_vec("") - EP) :/ (-2 * `eps')
	quietly replace w = 1 in `h'
}
mata:
f    = `n' / (`n' - 1)
ratb = SEB :/ sqrt(f :* colsum(UB:^2))
nm   = tokens(st_local("cn"))
dl   = selectindex(strpos(nm, ":delta") :> 0)
sb   = J(1, cols(nm), 1); sb[dl] = J(1, cols(dl), 0)
printf("\n{txt}Table 9, brute force: reported SE / sqrt(n/(n-1) sum U_h^2)\n")
printf("{txt}  %-30s %8s %8s %8s %5s\n", "", "median", "min", "max", "cells")
_line("latent coefficients", ratb[selectindex(sb)])
_line("delta", ratb[dl])
_el_lines(SEE :/ sqrt(f :* colsum(UE:^2)))
end
log close t9a
di as res _n "Table9_bruteforce: done"
