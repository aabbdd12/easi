*! audit_derivations.do -- are the elasticity formulas of the R package right?
*!
*! Independent check by finite differences.  The EASI model is a fixed point in
*! (w, y) for given (p, x, z, eps):
*!
*!     w_j = c_j + sum_r b_rj y^r + sum_t g_tj z_t + sum_k a_jk(z) l_k
*!           + y sum_k b_jk l_k + y sum_t h_tj z_t + eps_j
*!     y   = [ln x - sum_j w_j l_j + 1/2 sum_jk a_jk(z) l_j l_k]
*!           / [1 - 1/2 sum_jk b_jk l_j l_k]
*!
*! Perturb l_i (or ln x), re-solve the fixed point, and read off dw_j/dl_i.
*! That is the elasticity by construction, with no derivation to get wrong.
*! Compare with what the package reports, and with the analytic formula that
*! this file also derives.
*!
*! The coefficient blocks are rebuilt here from e(b) rather than reused from
*! easi.ado, so this is an independent reimplementation, not a self-check.


clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table12.do
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

log using "`ROOT'/replication/out/audit.log", replace text name(au)

*==========================================================================
mata:
mata set matastrict on

struct blocks {
	real matrix    cc, bjr, gjt, hjt, bjk
	pointer(real matrix) rowvector A
	real scalar    R, T, neq, Jg
}

// rebuild every coefficient block of the full J-good system from e(b),
// completing the dropped good by adding up
struct blocks scalar mkblocks(real colvector b, real scalar R, real scalar T,
	real scalar neq, real scalar py, real scalar zy, real scalar pz,
	real rowvector ipz)
{
	struct blocks scalar B
	real matrix Bm
	real scalar i, j, r, t, s, q, k, onp, oynp, onpz, acc

	B.R = R; B.T = T; B.neq = neq; B.Jg = neq + 1
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? cols(ipz) * neq : 0)
	Bm   = rowshape(b', neq)'

	B.cc = J(1, neq + 1, 0)
	for (i = 1; i <= neq; i++) B.cc[i] = Bm[1, i]
	B.cc[neq + 1] = 1 - sum(B.cc[1::neq])

	B.bjr = J(R, neq + 1, 0)
	for (i = 1; i <= neq; i++) for (r = 1; r <= R; r++) B.bjr[r, i] = Bm[1 + r, i]
	for (r = 1; r <= R; r++) B.bjr[r, neq + 1] = -sum(B.bjr[r, 1::neq])

	B.gjt = J(T, neq + 1, 0)
	for (i = 1; i <= neq; i++) for (t = 1; t <= T; t++) B.gjt[t, i] = Bm[1 + R + t, i]
	for (t = 1; t <= T; t++) B.gjt[t, neq + 1] = -sum(B.gjt[t, 1::neq])

	B.hjt = J(T, neq + 1, 0)
	if (zy) {
		for (i = 1; i <= neq; i++) for (t = 1; t <= T; t++)
			B.hjt[t, i] = Bm[1 + R + T + t, i]
		for (t = 1; t <= T; t++) B.hjt[t, neq + 1] = -sum(B.hjt[t, 1::neq])
	}

	B.bjk = J(neq + 1, neq + 1, 0)
	if (py) {
		for (i = 1; i <= neq; i++) for (j = 1; j <= neq; j++)
			B.bjk[j, i] = Bm[oynp + j, i]
		for (j = 1; j <= neq; j++) B.bjk[j, neq + 1] = -sum(B.bjk[j, 1::neq])
		for (j = 1; j <= neq; j++) B.bjk[neq + 1, j] = B.bjk[j, neq + 1]
		B.bjk[neq + 1, neq + 1] = -sum(B.bjk[neq + 1, 1::neq])
	}

	B.A = J(1, neq + 1, NULL)
	for (i = 1; i <= neq + 1; i++) B.A[i] = &(J(T + 1, neq + 1, 0))
	for (i = 1; i <= neq; i++) for (j = 1; j <= neq; j++)
		(*B.A[i])[1, j] = Bm[onp + j, i]
	for (j = 1; j <= neq; j++) {
		acc = 0
		for (i = 1; i <= neq; i++) acc = acc + (*B.A[i])[1, j]
		(*B.A[neq + 1])[1, j] = -acc
	}
	for (i = 1; i <= neq; i++) (*B.A[i])[1, neq + 1] = (*B.A[neq + 1])[1, i]
	(*B.A[neq + 1])[1, neq + 1] = -sum((*B.A[neq + 1])[1, 1::neq])
	if (pz) {
		for (i = 1; i <= neq; i++) for (s = 1; s <= cols(ipz); s++)
			for (q = 1; q <= neq; q++)
				(*B.A[i])[ipz[s] + 1, q] = Bm[onpz + (s - 1) * neq + q, i]
		for (s = 1; s <= cols(ipz); s++) {
			t = ipz[s] + 1
			for (i = 1; i <= neq; i++) {
				acc = 0
				for (j = 1; j <= neq; j++) acc = acc + (*B.A[j])[t, i]
				(*B.A[neq + 1])[t, i] = -acc
			}
			for (i = 1; i <= neq; i++) (*B.A[i])[t, neq + 1] = (*B.A[neq + 1])[t, i]
			acc = 0
			for (i = 1; i <= neq; i++) acc = acc + (*B.A[i])[t, neq + 1]
			(*B.A[neq + 1])[t, neq + 1] = -acc
		}
	}
	return(B)
}

// a_jk(z), the J x J price matrix of one observation-set, for good i
real matrix azk(struct blocks scalar B, real matrix Zc, real scalar i)
{
	real matrix M
	real scalar k, t
	M = J(rows(Zc), B.Jg, 0)
	for (k = 1; k <= B.Jg; k++) {
		for (t = 1; t <= B.T + 1; t++) M[, k] = M[, k] + (*B.A[i])[t, k] :* Zc[, t]
	}
	return(M)
}

// the systematic part of the shares at a given y
real matrix modelw(struct blocks scalar B, real matrix L, real colvector y,
	real matrix z, real matrix Zc)
{
	real matrix W, az
	real scalar i, r, t, k, n
	n = rows(L)
	W = J(n, B.Jg, 0)
	for (i = 1; i <= B.Jg; i++) {
		az = azk(B, Zc, i)
		W[, i] = J(n, 1, B.cc[i])
		for (r = 1; r <= B.R; r++) W[, i] = W[, i] + B.bjr[r, i] :* y:^r
		for (t = 1; t <= B.T; t++) W[, i] = W[, i] + B.gjt[t, i] :* z[, t]
		for (t = 1; t <= B.T; t++) W[, i] = W[, i] + B.hjt[t, i] :* z[, t] :* y
		for (k = 1; k <= B.Jg; k++) W[, i] = W[, i] + az[, k] :* L[, k]
		for (k = 1; k <= B.Jg; k++) W[, i] = W[, i] + B.bjk[k, i] :* L[, k] :* y
	}
	return(W)
}

// the implicit utility index implied by a given set of shares
real colvector yfrom(struct blocks scalar B, real matrix L, real colvector lnx,
	real matrix w, real matrix Zc)
{
	real matrix az
	real colvector aq, bq
	real scalar i, k, n
	n = rows(L)
	aq = J(n, 1, 0); bq = J(n, 1, 0)
	for (i = 1; i <= B.Jg; i++) {
		az = azk(B, Zc, i)
		for (k = 1; k <= B.Jg; k++) aq = aq + az[, k] :* L[, i] :* L[, k]
		for (k = 1; k <= B.Jg; k++) bq = bq + B.bjk[k, i] :* L[, i] :* L[, k]
	}
	return((lnx - rowsum(w :* L) + 0.5 :* aq) :/ (1 :- 0.5 :* bq))
}

// solve the (w, y) fixed point for given log prices / log expenditure
void solvefp(struct blocks scalar B, real matrix L, real colvector lnx,
	real matrix z, real matrix Zc, real matrix eps, real matrix w, real colvector y)
{
	real matrix az
	real colvector num, den, ysum, aq, bq
	real scalar it, i, r, t, k, n

	n = rows(L)
	if (rows(w) != n) {
		w = J(n, B.Jg, 1 / B.Jg)
		y = lnx
	}
	for (it = 1; it <= 500; it++) {
		// y from the current shares
		aq = J(n, 1, 0); bq = J(n, 1, 0)
		for (i = 1; i <= B.Jg; i++) {
			az = azk(B, Zc, i)
			for (k = 1; k <= B.Jg; k++) aq = aq + az[, k] :* L[, i] :* L[, k]
			for (k = 1; k <= B.Jg; k++) bq = bq + B.bjk[k, i] :* L[, i] :* L[, k]
		}
		num = lnx - rowsum(w :* L) + 0.5 :* aq
		den = 1 :- 0.5 :* bq
		ysum = num :/ den
		if (it > 1 & max(abs(ysum - y)) < 1e-13) {
			y = ysum
			break
		}
		y = ysum
		// shares from the current y
		for (i = 1; i <= B.Jg; i++) {
			az = azk(B, Zc, i)
			w[, i] = J(n, 1, B.cc[i]) + eps[, i]
			for (r = 1; r <= B.R; r++) w[, i] = w[, i] + B.bjr[r, i] :* y:^r
			for (t = 1; t <= B.T; t++) w[, i] = w[, i] + B.gjt[t, i] :* z[, t]
			for (t = 1; t <= B.T; t++) w[, i] = w[, i] + B.hjt[t, i] :* z[, t] :* y
			for (k = 1; k <= B.Jg; k++) w[, i] = w[, i] + az[, k] :* L[, k]
			for (k = 1; k <= B.Jg; k++) w[, i] = w[, i] + B.bjk[k, i] :* L[, k] :* y
		}
	}
}

void audit(string scalar svars, string scalar pvars, string scalar xvar,
	string scalar zvars, string scalar sipz, real scalar R,
	real scalar py, real scalar zy, real scalar pz)
{
	struct blocks scalar B
	real matrix s, L, z, Zc, eps, w0, wp, wm, Lp, az, RP, RI, CP, CI
	real matrix NUMP, ANAP, PKGP
	real colvector lnx, y0, yp, ym, ws, Ci, Di, Gi, Dn2, tot2
	real rowvector ipz, NUMI, ANAI, PKGI
	real scalar n, Jg, neq, T, i, j, k, r, t, d, del, mx1, mx2

	s   = st_data(., svars, 0)
	L   = st_data(., pvars, 0)
	lnx = st_data(., xvar,  0)
	z   = st_data(., zvars, 0)
	ipz = strtoreal(tokens(sipz))
	n = rows(s); Jg = cols(s); neq = Jg - 1; T = cols(z)
	Zc = J(n, 1, 1), z

	B = mkblocks(st_matrix("e(b)")', R, T, neq, py, zy, pz, ipz)

	// y is directly observable from the shares, so the residuals follow in one
	// step: eps = observed shares - model part evaluated at that y.  (s, y0) is
	// then a fixed point of the system by construction; the call below checks it.
	y0  = yfrom(B, L, lnx, s, Zc)
	eps = s - modelw(B, L, y0, z, Zc)
	w0  = s
	solvefp(B, L, lnx, z, Zc, eps, w0, y0)

	printf("\n{hline 78}\n")
	printf("{txt}  AUDIT DES DERIVATIONS -- verification par differences finies\n")
	printf("{hline 78}\n")
	printf("{txt}  Point fixe resolu : max|w - w observe| = {res}%9.2e{txt}\n",
		max(abs(w0 - s)))
	printf("{txt}  (doit etre ~0 : confirme que la representation a J biens est\n")
	printf("{txt}   equivalente au systeme estime a J-1 equations)\n\n")

	del = 1e-5

	// ---- numerical price elasticities -------------------------------------
	NUMP = J(Jg, Jg, 0)
	for (k = 1; k <= Jg; k++) {
		Lp = L; Lp[, k] = L[, k] :+ del
		wp = w0; yp = y0
		solvefp(B, Lp, lnx, z, Zc, eps, wp, yp)
		Lp = L; Lp[, k] = L[, k] :- del
		wm = w0; ym = y0
		solvefp(B, Lp, lnx, z, Zc, eps, wm, ym)
		for (j = 1; j <= Jg; j++) {
			// aggregated the same way the package does: mean(dw)/mean(w)
			NUMP[k, j] = mean((wp[, j] - wm[, j]) :/ (2 * del)) / mean(s[, j]) -
				(j == k)
		}
	}

	// ---- numerical expenditure elasticities --------------------------------
	wp = w0; yp = y0
	solvefp(B, L, lnx :+ del, z, Zc, eps, wp, yp)
	wm = w0; ym = y0
	solvefp(B, L, lnx :- del, z, Zc, eps, wm, ym)
	NUMI = J(1, Jg, 0)
	for (j = 1; j <= Jg; j++) {
		NUMI[j] = 1 + mean((wp[, j] - wm[, j]) :/ (2 * del)) / mean(s[, j])
	}

	// ---- the analytic formula this audit derives ---------------------------
	// dy/dl_i = -w_i / Dn2,  Dn2 = 1 + sum_j (C_j + D_j) l_j + 1/2 sum_jk b_jk l_j l_k
	// dy/dlnx = 1 / Dn2
	Dn2 = J(n, 1, 1)
	for (j = 1; j <= Jg; j++) {
		Ci = J(n, 1, 0)
		for (r = 1; r <= R; r++) Ci = Ci + r * B.bjr[r, j] :* y0:^(r - 1)
		Di = J(n, 1, 0)
		if (zy) for (t = 1; t <= T; t++) Di = Di + B.hjt[t, j] :* z[, t]
		Dn2 = Dn2 + (Ci + Di) :* L[, j]
	}
	tot2 = J(n, 1, 0)
	for (j = 1; j <= Jg; j++) for (k = 1; k <= Jg; k++)
		tot2 = tot2 + B.bjk[j, k] :* L[, j] :* L[, k]
	Dn2 = Dn2 + 0.5 :* tot2

	ANAP = J(Jg, Jg, 0); ANAI = J(1, Jg, 0)
	for (j = 1; j <= Jg; j++) {
		Ci = J(n, 1, 0)
		for (r = 1; r <= R; r++) Ci = Ci + r * B.bjr[r, j] :* y0:^(r - 1)
		Di = J(n, 1, 0)
		if (zy) for (t = 1; t <= T; t++) Di = Di + B.hjt[t, j] :* z[, t]
		Gi = J(n, 1, 0)
		if (py) for (k = 1; k <= Jg; k++) Gi = Gi + B.bjk[k, j] :* L[, k]
		ANAI[j] = 1 + mean((Ci + Di + Gi) :/ Dn2) / mean(s[, j])
		az = azk(B, Zc, j)
		for (k = 1; k <= Jg; k++) {
			ANAP[k, j] = mean((Ci + Di + Gi) :* (-s[, k] :/ Dn2)
				+ az[, k] + B.bjk[k, j] :* y0) / mean(s[, j]) - (j == k)
		}
	}

	printf("{txt}  Facteur Dn2 = 1 + sum_j (C_j+D_j) l_j + 0.5 sum_jk b_jk l_j l_k\n")
	printf("{txt}    (= 1 exactement si tous les log-prix sont nuls ; le paquet R\n")
	printf("{txt}     suppose implicitement Dn2 = 1 partout)\n")
	printf("{txt}    moyenne = {res}%8.5f{txt}   min = {res}%8.5f{txt}   max = {res}%8.5f\n\n",
		mean(Dn2), min(Dn2), max(Dn2))

	RP = st_matrix("RP"); RI = st_matrix("RI")
	CP = st_matrix("CP"); CI = st_matrix("CI")

	printf("{txt}  ELASTICITES-DEPENSE\n")
	printf("{txt}  %-10s %12s %12s %12s\n", "bien", "num. (ref)",
		"easi corrige", "paquet R")
	for (j = 1; j <= Jg; j++) {
		printf("{txt}  %-12.0f {res}%12.5f %12.5f %12.5f\n",
			j, NUMI[j], CI[1, j], RI[1, j])
	}
	printf("\n{txt}  max |analytique - numerique| = {res}%9.2e\n",
		max(abs(ANAI - NUMI)))
	printf("{txt}  max |easi corrige - numerique| = {res}%9.2e\n",
		max(abs(CI - NUMI)))
	printf("{txt}  max |paquet R  - numerique| = {res}%9.2e\n",
		max(abs(RI - NUMI)))

	printf("\n{txt}  ELASTICITES-PRIX (9x9)\n")
	mx1 = max(abs(ANAP - NUMP))
	mx2 = max(abs(RP - NUMP))
	printf("{txt}  max |analytique - numerique| = {res}%9.2e\n", mx1)
	printf("{txt}  max |easi corrige - numerique| = {res}%9.2e\n",
		max(abs(CP - NUMP)))
	printf("{txt}  max |paquet R  - numerique| = {res}%9.2e\n", mx2)
	printf("\n{txt}  diagonale (elasticites-prix propres) :\n")
	printf("{txt}  %-12s %12s %12s %12s\n", "bien", "num. (ref)", "easi corrige",
		"paquet R")
	for (j = 1; j <= Jg; j++) {
		printf("{txt}  %-12.0f {res}%12.5f %12.5f %12.5f\n",
			j, NUMP[j, j], CP[j, j], RP[j, j])
	}
	printf("\n{hline 78}\n")
}

end

use "`ROOT'/examples/hixdata.dta", clear
local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local DE age hsex carown time tran

* the reference specification, in compat mode = the R package bit for bit
* legacy(all) = the R package, formulas and all
qui easi `SH', lnprices(`PR') lnexpenditure(log_y) demographics(`DE')	///
	power(5) py pz zy legacy(all) nolog

matrix RP = e(elast_price)
matrix RI = e(elast_income)

* the same coefficients, but the corrected elasticity formulas: every
* legacy behaviour is kept EXCEPT the derivation and the 1e-6 quantisation
* (which the finite-difference benchmark does not apply either)
qui easi `SH', lnprices(`PR') lnexpenditure(log_y) demographics(`DE')	///
	power(5) py pz zy legacy(interpz instr ez lastse sigma1) nolog
matrix CP = e(elast_price)
matrix CI = e(elast_income)
local ipz "1"

mata: audit("`SH'", "`PR'", "log_y", "`DE'", "`ipz'", `e(power)', 1, 1, 1)

log close au

