*! test_step11.do -- is the generated-regressor Jacobian right, and how big?
*!
*! The design holds y, which is a function of beta, so the 3SLS covariance is
*! the covariance of the last linear step CONDITIONAL on y.  Pendakur's code
*! flags this twice:
*!     *note that reported standard errors are wrong for iterated estimates
*!
*! -easi- now adds the missing term.  Writing the estimator as GMM in beta with
*! fixed instruments,
*!     m_j(beta) = sum_i wt_i xhat_i sum_m sigma^{jm} e_im(beta),
*!     e_im(beta) = w_im - x_i(y(beta))'b_m,
*! the Jacobian is  -dm/dbeta' = G = (Sigma^-1 # A) + C, and V = P_G B P_G'.
*! The conditional covariance is C = 0.
*!
*! This file checks the analytic G by REIMPLEMENTING the moment vector from
*! scratch -- design matrix, instruments, y(beta), the lot -- and
*! differentiating it numerically.  Nothing from easi.ado is reused except
*! e(b), e(Sigma) and the raw data, so agreement is an independent check and
*! not a tautology.  -legacy(instr)- is used so that the instruments are the
*! single-pass ones and can be rebuilt in one line; the question here is the
*! Jacobian, not the instrument.

clear all
set more off

*------------------------------------------------------------------ Mata
* Defined before -clear all- would wipe them, and with NO top-level
* -version- command: setting the version at do-file scope makes easi fail
* with r(509) further down.
mata:
mata clear

// independent rebuild of the design and of the instruments
real matrix _st11_X(real colvector y, real matrix z, real matrix np,
	real scalar R, real scalar py, real scalar zy, real scalar pz)
{
	real matrix X
	real scalar i, t, neq
	neq = cols(np)
	X   = J(rows(y), 1, 1)
	for (i = 1; i <= R; i++) X = X, y:^i
	X = X, z
	if (zy) X = X, (y :* z)
	X = X, np
	if (py) X = X, (y :* np)
	if (pz) for (t = 1; t <= cols(z); t++)
			for (i = 1; i <= neq; i++) X = X, (np[, i] :* z[, t])
	return(X)
}
real matrix _st11_Z(real colvector yp, real matrix z, real matrix np,
	real scalar R, real scalar py, real scalar zy, real scalar pz)
{
	real matrix Z
	real scalar i, t, neq
	neq = cols(np)
	Z   = J(rows(yp), 1, 1), z, np
	if (pz) for (t = 1; t <= cols(z); t++)
			for (i = 1; i <= neq; i++) Z = Z, (np[, i] :* z[, t])
	for (i = 1; i <= R; i++) Z = Z, yp:^i
	if (zy) Z = Z, (yp :* z)
	if (py) Z = Z, (yp :* np)
	return(Z)
}

// y as a function of beta, straight from the definition
real colvector _st11_y(real colvector b, real colvector ystone, real matrix np,
	real matrix z, real scalar R, real scalar py, real scalar pz,
	real scalar zy)
{
	real matrix B, Aco
	real colvector pAp, pBp
	real scalar neq, n, k, T, onp, oynp, onpz, i, j, t

	neq = cols(np); n = rows(np); T = cols(z)
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = rows(b) / neq
	B    = rowshape(b', neq)'
	pAp  = J(n, 1, 0); pBp = J(n, 1, 0)
	Aco  = B[| onp + 1, 1 \ onp + neq, neq |]'
	for (i = 1; i <= neq; i++)
		for (j = 1; j <= neq; j++)
			pAp = pAp + Aco[i, j] :* np[, i] :* np[, j]
	if (pz) {
		for (t = 1; t <= T; t++) {
			Aco = B[| onpz + (t-1)*neq + 1, 1 \ onpz + t*neq, neq |]'
			for (i = 1; i <= neq; i++)
				for (j = 1; j <= neq; j++)
					pAp = pAp + Aco[i, j] :*
						np[, i] :* np[, j] :* z[, t]
		}
	}
	if (py) {
		Aco = B[| oynp + 1, 1 \ oynp + neq, neq |]'
		for (i = 1; i <= neq; i++)
			for (j = 1; j <= neq; j++)
				pBp = pBp + Aco[i, j] :* np[, i] :* np[, j]
	}
	return((ystone :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp))
}

// the moment vector, equation-major
real colvector _st11_m(real colvector b, real matrix Xhat, real matrix Y,
	real colvector ystone, real matrix np, real matrix z, real matrix Si,
	real colvector wt, real scalar R, real scalar py, real scalar pz,
	real scalar zy)
{
	real matrix X, E
	real colvector y
	y = _st11_y(b, ystone, np, z, R, py, pz, zy)
	X = _st11_X(y, z, np, R, py, zy, pz)
	E = Y - X * rowshape(b', cols(Y))'
	return(vec(quadcross(Xhat, wt, E) * Si))
}

// the Slutsky symmetry restrictions, rebuilt from the block offsets
real matrix _st11_R(real scalar neq, real scalar k, real scalar onp,
	real scalar py, real scalar oynp, real scalar pz, real scalar onpz,
	real scalar T)
{
	real matrix Rr
	real rowvector row
	real scalar i, j, t, o

	Rr = J(0, neq * k, 0)
	for (i = 1; i <= neq; i++) {
		for (j = i + 1; j <= neq; j++) {
			o = onp
			row = J(1, neq * k, 0)
			row[(i - 1) * k + o + j] =  1
			row[(j - 1) * k + o + i] = -1
			Rr = Rr \ row
			if (py) {
				o = oynp
				row = J(1, neq * k, 0)
				row[(i - 1) * k + o + j] =  1
				row[(j - 1) * k + o + i] = -1
				Rr = Rr \ row
			}
			if (pz) {
				for (t = 1; t <= T; t++) {
					o = onpz + (t - 1) * neq
					row = J(1, neq * k, 0)
					row[(i - 1) * k + o + j] =  1
					row[(j - 1) * k + o + i] = -1
					Rr = Rr \ row
				}
			}
		}
	}
	return(Rr)
}

void _st11(string scalar svars, string scalar lpvars, string scalar lxvar,
	string scalar zvars, real scalar R, real scalar py, real scalar pz,
	real scalar zy)
{
	real matrix s, p, z, np, Y, X, Z, Xhat, A, Si, Rr, W, P, G, PG, Vc, Vg
	real matrix ZZiZX, Gnum
	real colvector lnx, ystone, ytil, b, wt, m0, bp, bm, mp, mm, y
	real rowvector ms
	real scalar n, J, neq, T, k, K, r, i, h, onp, oynp, onpz
	real colvector sc, sg

	s   = st_data(., svars);   p = st_data(., lpvars)
	lnx = st_data(., lxvar);   z = st_data(., zvars)
	n = rows(s); J = cols(s); neq = J - 1; T = cols(z)
	wt = J(n, 1, 1)

	np     = p[, 1::neq] :- p[, J]
	ystone = lnx - rowsum(s :* p)
	ms     = mean(s)
	ytil   = lnx - p * ms'
	Y      = s[, 1::neq]

	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? T * neq : 0)
	K    = neq * k

	b  = st_matrix("bb")'
	Si = invsym(st_matrix("Sg"))
	Rr = _st11_R(neq, k, onp, py, oynp, pz, onpz, T)
	r  = rows(Rr)

	// instruments: single pass, built on y_tilda  (-legacy(instr)-)
	Z     = _st11_Z(ytil, z, np, R, py, zy, pz)
	y     = _st11_y(b, ystone, np, z, R, py, pz, zy)
	X     = _st11_X(y, z, np, R, py, zy, pz)
	ZZiZX = cholsolve(quadcross(Z, wt, Z), quadcross(Z, wt, X))
	Xhat  = Z * ZZiZX
	A     = quadcross(Xhat, wt, Xhat)

	// the estimate must zero the moment vector, up to the restrictions
	m0 = _st11_m(b, Xhat, Y, ystone, np, z, Si, wt, R, py, pz, zy)
	st_numscalar("mres", max(abs(m0 - Rr' * qrsolve(Rr', m0))) / max(abs(m0)))

	// conditional covariance, to validate the reimplementation
	W  = (Si # A), Rr' \ Rr, J(r, r, 0)
	P  = luinv(W)[| 1, 1 \ K, K |]
	P  = (P + P') / 2
	Vc = st_matrix("Vcd")
	st_numscalar("dcond", max(abs(P - Vc)) / max(abs(Vc)))

	// numerical Jacobian of the moment vector
	Gnum = J(K, K, 0)
	for (i = 1; i <= K; i++) {
		h = 1e-6 * max((1, abs(b[i])))
		bp = b; bp[i] = bp[i] + h
		bm = b; bm[i] = bm[i] - h
		mp = _st11_m(bp, Xhat, Y, ystone, np, z, Si, wt, R, py, pz, zy)
		mm = _st11_m(bm, Xhat, Y, ystone, np, z, Si, wt, R, py, pz, zy)
		Gnum[, i] = -(mp - mm) / (2 * h)
	}
	W  = (Gnum, Rr' \ Rr, J(r, r, 0))
	PG = luinv(W)[| 1, 1 \ K, K |]
	G  = PG * (Si # A) * PG'
	G  = (G + G') / 2
	Vg = st_matrix("Vgn")
	st_numscalar("dgen", max(abs(G - Vg)) / max(abs(Vg)))

	sc = sqrt(diagonal(Vc)); sg = sqrt(diagonal(Vg))
	st_numscalar("rmin", min(sg :/ sc))
	st_numscalar("rmax", max(sg :/ sc))
}
end

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do 06_generated_regressor.do
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

log using "`ROOT'/replication/out/step11.log", replace text name(ta)

use "`ROOT'/examples/hixdata.dta", clear

local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers

di ""
di as txt "{hline 76}"
di as txt "  Regresseur genere : le jacobien analytique contre les differences finies"
di as txt "{hline 76}"

local fail 0

foreach spec in "light" "inter" {
	if "`spec'" == "light" {
		local XTRA demographics(age hsex carown) power(3)
		local PY 0
		local PZ 0
		local ZY 0
	}
	else {
		local XTRA demographics(age hsex carown) power(3) py pz zy
		local PY 1
		local PZ 1
		local ZY 1
	}
	local OPT lnprices(`PR') lnexpenditure(log_y) `XTRA' nolog noelastse	///
		  vce(conventional)

	qui easi `SH', `OPT' legacy(instr condse)
	matrix bb  = e(b)
	matrix Sg  = e(Sigma)
	matrix Vcd = e(V)
	local R    = e(power)

	qui easi `SH', `OPT' legacy(instr)
	matrix Vgn = e(V)

	mata: _st11("`SH'", "`PR'", "log_y", "age hsex carown", `R',		///
		`PY', `PZ', `ZY')

	di ""
	di as txt "  specification : " as res "`spec'" as txt			///
	   "   (py=`PY' pz=`PZ' zy=`ZY')"
	di as txt "    m(beta_hat) reconstruit, norme max" _col(52)		///
	   as res %11.3e mres
	di as txt "    V conditionnelle, ecart relatif max" _col(52)		///
	   as res %11.3e dcond
	di as txt "    V corrigee contre differences finies" _col(52)		///
	   as res %11.3e dgen
	di as txt "    se(corrigee)/se(conditionnelle), min .. max" _col(46)	///
	   as res %8.4f rmin as txt " .." as res %8.4f rmax

	if dcond > 1e-6 | dgen > 1e-5 | mres > 1e-4 local fail = `fail' + 1
}

di ""
if `fail' == 0 {
	di as txt "  -> le jacobien analytique EST celui du systeme d'estimation."
	di as txt "     Le terme manquant est reel mais petit : y ne depend de beta"
	di as txt "     qu'a travers p'Ap/2 et p'Bp/2, et les prix normalises sont"
	di as txt "     d'ordre 0.1, donc dy/dbeta est d'ordre 0.005."
}
else {
	di as error "  -> `fail' specification(s) en echec"
}
di as txt "{hline 76}"

log close ta
