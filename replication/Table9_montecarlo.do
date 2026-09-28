*! Table9_montecarlo.do -- Table 9 of the technical note, third and fourth
*! columns: Monte Carlo of the analytic standard errors under selection on
*! the censored perfect DGP (dgp_easi.do, censor): 700 samples of 10,000
*! households (seeds 910001 ... 910700), SRS, no weight, vce(robust), an
*! exclusion restriction per good.  For each block -- the latent
*! coefficients, delta, the probit coefficients, the elasticities of the
*! households and of the market -- the mean analytic SE divided by the SD of
*! the estimates (median, min, max over the cells), the cells outside the
*! band 1.96 sqrt((k-1)/(4R)), the bias against the population values (the
*! same estimator on 500,000 households of the DGP) in units of the SD, and
*! the rejection rate of the Wald test at 5% of H0: theta = population
*! value.  About seventy minutes.  The samples and the population values are
*! saved in out/selection_montecarlo.mmat before the summary; with
*! global MCSUM 1, the script redoes the summary from that file only.

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table9_montecarlo.do
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
capture log close t9c
log using "`OUT'/selection_montecarlo.log", replace text name(t9c)

local R 700
local seed0 910000
local n 10000
local power 2
local SPEC "w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(`power') selvars(w1: q1 ; w2: q2) vce(robust) nolog"

mata:
// one estimation, stacked: e(b), the probit coefficients, the elasticities
// (households then market); se = 1 for the standard errors in the same layout
real rowvector _mc_vec(real scalar se)
{
	real rowvector a, v
	real scalar i
	string rowvector sf
	string scalar s
	s = (se ? "_se" : "")
	if (se) v = sqrt(diagonal(st_matrix("e(V_sel)")))'
	else {
		a = vec(st_matrix("e(sel_alpha)"))'
		v = st_matrix("e(b)"), select(a, a :< .)
	}
	sf = ("", "_mkt")
	for (i = 1; i <= 2; i++) {
		v = v, st_matrix("e(elast_exp" + sf[i] + s + ")"),
			vec(st_matrix("e(elast_price_nc" + sf[i] + s + ")"))',
			vec(st_matrix("e(elast_price_c" + sf[i] + s + ")"))',
			vec(st_matrix("e(elast_demo" + sf[i] + s + ")"))'
	}
	return(v)
}
E = J(0, 0, .); S = J(0, 0, .); NN = J(0, 1, .)
end

if "$MCSUM" != "1" {
* the samples
local nfail 0
forvalues r = 1/`R' {
	quietly easi_dgp, n(`n') seed(`=`seed0' + `r'') power(`power') censor
	local nneg = r(nneg)
	capture quietly easi `SPEC' compensated notable
	if _rc | e(converged) != 1 local ++nfail
	else {
		mata: E = (rows(E) ? E \ _mc_vec(0) : _mc_vec(0))
		mata: S = (rows(S) ? S \ _mc_vec(1) : _mc_vec(1))
		mata: NN = NN \ `nneg'
	}
}
di as txt "`R' samples of `n' households, `nfail' failure(s)"

* the population values
quietly easi_dgp, n(500000) seed(987654321) power(`power') censor
di as txt "population: 500,000 households, buyers with a nonpositive latent share: " r(nneg)
quietly easi `SPEC' noelastse notable
di as txt "population delta: " %8.5f _b[w1:delta] " " %8.5f _b[w2:delta] " (true 0.04)"
local cn : colfullnames e(b)
mata: POP = _mc_vec(0); CN = tokens(st_local("cn"))
mata: mata matsave "`OUT'/selection_montecarlo" E S NN POP CN, replace
}
else mata: mata matuse "`OUT'/selection_montecarlo", replace

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
void _mc_line(string scalar nm, real rowvector r, real rowvector band,
	real rowvector bias, real rowvector rej, real scalar bb)
{
	printf("{txt}  %-28s {res}%7.3f %7.3f %7.3f %5.0f %7.0f %8.3f %7.2f %5.0f\n", nm,
		_med(r'), min(r), max(r), cols(r), sum(abs(r :- 1) :> band),
		max(abs(bias)), 100 * mean(rej'), sum(abs(rej :- 0.05) :> bb))
}
R   = rows(E)
sd  = sqrt(diagonal(quadvariance(E)))'
ms  = mean(S)
ku  = J(1, cols(E), .)
for (i = 1; i <= cols(E); i++) {
	ku[i] = mean((E[, i] :- mean(E[, i])):^4) / (mean((E[, i] :- mean(E[, i])):^2)^2)
}
rat  = ms :/ sd
band = 1.96 :* sqrt((ku :- 1) :/ (4 * R))
bias = (mean(E) - POP) :/ sd
rej  = mean((abs(E :- POP) :/ S) :> invnormal(0.975))
bb   = 1.96 * sqrt(0.05 * 0.95 / R)
K    = cols(CN)
NAL  = cols(E) - K - 54
didx = selectindex(strpos(CN, ":delta") :> 0)
sel  = J(1, K, 1); sel[didx] = J(1, cols(didx), 0); bidx = selectindex(sel)
printf("\n{txt}Table 9, Monte Carlo under selection, R = %f samples (buyers with a nonpositive latent share: %5.2f per sample)\n", R, mean(NN))
printf("{txt}  %-28s %7s %7s %7s %5s %7s %8s %7s %5s\n", "block (SE/SD)", "median", "min", "max",
	"cells", "outside", "max|b|/SD", "Wald %", "out")
_mc_line("latent coefficients", rat[bidx], band[bidx], bias[bidx], rej[bidx], bb)
_mc_line("delta", rat[didx], band[didx], bias[didx], rej[didx], bb)
_mc_line("probit coefficients", rat[| K + 1 \ K + NAL |], band[| K + 1 \ K + NAL |],
	bias[| K + 1 \ K + NAL |], rej[| K + 1 \ K + NAL |], bb)
bn = ("expenditure", "uncompensated", "compensated", "demographic")
tn = ("households", "market")
off = K + NAL
for (ty = 1; ty <= 2; ty++) {
	for (bl = 1; bl <= 4; bl++) {
		len = (bl == 1 ? 3 : (bl <= 3 ? 9 : 6))
		lo  = off + 1
		hi  = lo + len - 1
		_mc_line(tn[ty] + ", " + bn[bl], rat[| lo \ hi |], band[| lo \ hi |],
			bias[| lo \ hi |], rej[| lo \ hi |], bb)
		off = hi
	}
}
printf("{txt}  all cells: median SE/SD %6.3f, mean |log SE/SD| %6.3f, Wald rejection %5.2f%%\n",
	_med(rat'), mean(abs(ln(rat))'), 100 * mean(rej'))
printf("{txt}  (band 1.96 sqrt((k-1)/(4R)): median %5.3f; Wald band 5%% +/- %4.2f; t pooled 2.5%% %6.3f, 97.5%% %6.3f)\n",
	_med(band'), 100 * bb, _q(vec((E :- POP) :/ S), 0.025), _q(vec((E :- POP) :/ S), 0.975))
end
log close t9c
di as res _n "Table9_montecarlo: done"
