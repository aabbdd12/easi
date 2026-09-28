*! test_step27.do -- easi 2.0.0, phase C: the point elasticities under
*! selection against finite differences of the expected demand.
*! Oracle, independent of the elasticity code of easi: for each household,
*!   E_i(l, ln x) = Phi(s_i'a_i) f_i(y, l, z) + delta_i phi(s_i'a_i),
*! the last good by adding up, with y solving the cost identity on the
*! latent shares completed for the non-buyers, their latent residual held
*! fixed:  y = ln x - sum_j (f_j(y) + r_j) l_j + p'Ap/2  (r_j = w~_j - f_j at
*! the data; w~ the observed share of a buyer, f - delta phi/(1 - Phi) for a
*! non-buyer).  Central differences in ln x and in each ln p_k, household by
*! household, averaged with the weights of the type:
*!   expenditure    1 + sum w dE_i/dln x / sum w E_i
*!   uncompensated  sum w dE_i/dln p_k / sum w E_i - 1{i = k}
*!   compensated    sum w (dE_i/dln p_k + E_k (E_i + dE_i/dln x)) / sum w E_i
*!                  - 1{i = k}   (Slutsky, household by household)
*!   demographic    sum w dE_i/dz_t / sum w  (at fixed y)
*! w = the weight (households) or the weight times x (market).
*! Specification: mex_bench [pw], power 2, no interaction, both cereals
*! corrected, age and isMale in the probits only.
clear all
set more off
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Section6-6_selection_elasticities.do
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

global NFAIL 0
capture program drop _chk
program define _chk
	args label value tol
	if (`value' <= `tol') di as txt "  ok    " as res %10.3g `value' as txt "  `label'"
	else {
		di as err "  FAIL  " %10.3g `value' "  `label' (tolerance `tol')"
		global NFAIL = ${NFAIL} + 1
	}
end

use "`ROOT'/examples/mex_bench.dta", clear
quietly easi w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) ///
	demographics(z1 z2) power(2) selvars(age isMale) vce(robust) nolog notable ///
	compensated tolerance(1e-13) iterate(1000)
quietly predict double yhat, y
quietly gen double om = sweight

mata:
// the latent shares of goods 1 and 2 at (y, l, z): columns 1-2
real matrix _f(real colvector y, real matrix L, real matrix Z, real matrix B)
{
	real matrix np
	np = L[, 1..2] :- L[, 3]
	return(J(rows(y), 1, B[1, .]) + y * B[2, .] + (y:^2) * B[3, .] + Z * B[4..5, .] + np * B[6..7, .])
}
// p'Ap with the symmetric price block of the two equations
real colvector _pap(real matrix L, real matrix B)
{
	real matrix np
	np = L[, 1..2] :- L[, 3]
	return(rowsum((np * B[6..7, .]) :* np))
}
// y solving the completed cost identity, the latent residuals R fixed
real colvector _y(real colvector lx, real matrix L, real matrix Z, real matrix B,
	real matrix R, real colvector y0)
{
	real colvector y, yo
	real matrix np
	real scalar it
	np = L[, 1..2] :- L[, 3]
	y  = y0
	for (it = 1; it <= 500; it++) {
		yo = y
		y  = lx - L[, 3] - rowsum((_f(y, L, Z, B) + R) :* np) + 0.5 :* _pap(L, B)
		if (max(abs(y - yo)) < 1e-15) break
	}
	return(y)
}
// the expected shares of the three goods at (l, ln x), y re-solved
real matrix _E(real colvector lx, real matrix L, real matrix Z, real matrix Q,
	real matrix B, real rowvector dl, real matrix AL, real matrix R,
	real colvector y0, real colvector y)
{
	real matrix S, E, F
	real colvector xb
	real scalar i
	y = _y(lx, L, Z, B, R, y0)
	F = _f(y, L, Z, B)
	E = J(rows(lx), 2, .)
	for (i = 1; i <= 2; i++) {
		S  = (J(rows(lx), 1, 1), L :- lx, Z, Q)
		xb = S * AL[, i]
		E[, i] = normal(xb) :* F[, i] + dl[i] :* normalden(xb)
	}
	return((E, 1 :- rowsum(E)))
}
end

* (the estimation sample)
quietly gen byte __es = e(sample)
mata:
b   = st_matrix("e(b)")
// per equation: _cons y1 y2 z1 z2 lp1 lp2 delta
B   = (b[1..7]', b[9..15]')
dl  = (b[8], b[16])
AL  = st_matrix("e(sel_alpha)")
W   = st_data(., ("w1", "w2", "w3"), "__es")
L   = st_data(., ("lp1", "lp2", "lp3"), "__es")
lx  = st_data(., "lx", "__es")
Z   = st_data(., ("z1", "z2"), "__es")
Q   = st_data(., ("age", "isMale"), "__es")
om  = st_data(., "om", "__es")
yh  = st_data(., "yhat", "__es")
n   = rows(W)
// the completed shares and the latent residuals at the data
F0  = _f(yh, L, Z, B)
R   = J(n, 2, .)
for (i = 1; i <= 2; i++) {
	S  = (J(n, 1, 1), L :- lx, Z, Q)
	xb = S * AL[, i]
	wt = (W[, i] :> 0) :* W[, i] + (W[, i] :<= 0) :* (F0[, i] - dl[i] :* normalden(xb) :/ normal(-xb))
	R[, i] = wt - F0[, i]
}
y1 = .
E0 = _E(lx, L, Z, Q, B, dl, AL, R, yh, y1)
st_numscalar("__dy", max(abs(y1 - yh)))
end
_chk "y re-solved at the data = predict, y" __dy 1e-9
quietly predict double __e*, shares
mata: st_numscalar("__de", max(abs(E0 - st_data(., ("__e1", "__e2", "__e3"), "__es"))))
_chk "E re-computed at the data = predict, shares" __de 1e-9

mata:
h  = 1e-6
yt = .
// d E / d ln x and d E / d ln p_k, household by household
DX = (_E(lx :+ h, L, Z, Q, B, dl, AL, R, yh, yt) - _E(lx :- h, L, Z, Q, B, dl, AL, R, yh, yt)) :/ (2 * h)
DP = J(3, 1, NULL)
for (k = 1; k <= 3; k++) {
	Lp = L; Lp[, k] = Lp[, k] :+ h
	Lm = L; Lm[, k] = Lm[, k] :- h
	DP[k] = &((_E(lx, Lp, Z, Q, B, dl, AL, R, yh, yt) - _E(lx, Lm, Z, Q, B, dl, AL, R, yh, yt)) :/ (2 * h))
}
// d E / d z_t at fixed y: the latent derivative and the probit channel
DZ = J(2, 1, NULL)
for (t = 1; t <= 2; t++) {
	Zp = Z; Zp[, t] = Zp[, t] :+ h
	Zm = Z; Zm[, t] = Zm[, t] :- h
	Ep = J(n, 2, .); Em = J(n, 2, .)
	Fp = _f(yh, L, Zp, B); Fm = _f(yh, L, Zm, B)
	for (i = 1; i <= 2; i++) {
		xp = (J(n, 1, 1), L :- lx, Zp, Q) * AL[, i]
		xm = (J(n, 1, 1), L :- lx, Zm, Q) * AL[, i]
		Ep[, i] = normal(xp) :* Fp[, i] + dl[i] :* normalden(xp)
		Em[, i] = normal(xm) :* Fm[, i] + dl[i] :* normalden(xm)
	}
	Ep = Ep, 1 :- rowsum(Ep); Em = Em, 1 :- rowsum(Em)
	DZ[t] = &((Ep - Em) :/ (2 * h))
}
// the oracle elasticities with the weights w
void _oracle(real colvector w, real matrix E0, real matrix DX, pointer colvector DP,
	pointer colvector DZ, real rowvector EX, real matrix EU, real matrix EC, real matrix EZ)
{
	real rowvector sE
	real scalar i, k, t
	sE = quadcolsum(w :* E0)
	EX = 1 :+ quadcolsum(w :* DX) :/ sE
	EU = J(3, 3, .); EC = J(3, 3, .)
	for (i = 1; i <= 3; i++) {
		for (k = 1; k <= 3; k++) {
			EU[i, k] = quadsum(w :* (*DP[k])[, i]) / sE[i] - (i == k)
			EC[i, k] = quadsum(w :* ((*DP[k])[, i] + E0[, k] :* (E0[, i] + DX[, i]))) / sE[i] - (i == k)
		}
	}
	EZ = J(2, 3, .)
	for (t = 1; t <= 2; t++) EZ[t, .] = quadcolsum(w :* (*DZ[t])) :/ quadsum(w)
}
EX = .; EU = .; EC = .; EZ = .
_oracle(om, E0, DX, DP, DZ, EX, EU, EC, EZ)
st_numscalar("__x",  mreldif(EX, st_matrix("e(elast_exp)")))
st_numscalar("__u",  mreldif(EU, st_matrix("e(elast_price_nc)")))
st_numscalar("__c",  mreldif(EC, st_matrix("e(elast_price_c)")))
st_numscalar("__z",  mreldif(EZ, st_matrix("e(elast_demo)")))
_oracle(om :* exp(lx), E0, DX, DP, DZ, EX, EU, EC, EZ)
st_numscalar("__xm", mreldif(EX, st_matrix("e(elast_exp_mkt)")))
st_numscalar("__um", mreldif(EU, st_matrix("e(elast_price_nc_mkt)")))
st_numscalar("__cm", mreldif(EC, st_matrix("e(elast_price_c_mkt)")))
st_numscalar("__zm", mreldif(EZ, st_matrix("e(elast_demo_mkt)")))
end
di as txt _n "the elasticities of easi against finite differences of the expected demand"
_chk "households, expenditure"   __x  1e-7
_chk "households, uncompensated" __u  1e-7
_chk "households, compensated"   __c  1e-7
_chk "households, demographic"   __z  1e-7
_chk "market, expenditure"       __xm 1e-7
_chk "market, uncompensated"     __um 1e-7
_chk "market, compensated"       __cm 1e-7
_chk "market, demographic"       __zm 1e-7
_chk "Engel and Cournot residuals, both types" `=e(chk_engel) + e(chk_cournot) + e(chk_engel_mkt) + e(chk_cournot_mkt)' 1e-12

di as txt _n "{hline 60}"
if ${NFAIL} == 0 di as res "test_step27: all checks passed"
else di as err "test_step27: ${NFAIL} check(s) failed"
