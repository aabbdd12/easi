*! dgp_easi.do -- a data-generating process for the EASI demand system, used
*! to validate the standard errors of easi (a perfect case: the model holds,
*! every parameter is strongly identified).
*!
*!   run "replication/dgp_easi.do"
*!   easi_dgp, n(3000) seed(1) [power(2) sde(0.05) censor rho(0.8) lift(0.03)]
*!
*! It leaves in memory three budget shares w1-w3, three log prices lp1-lp3,
*! log total expenditure lx and two demographics z1 (continuous) and z2
*! (binary). The model is the one easi estimates:
*!     w_j = b0_j + sum_r b_rj y^r + sum_t g_tj z_t + sum_k a_jk ln p_k + e_j
*! with adding-up (the constants sum to one, every other coefficient to zero
*! across goods), homogeneity (the rows of A sum to zero) and symmetry
*! (a_jk = a_kj); e is normal, homoskedastic, and sums to zero across goods.
*! As in the model, ln x is exogenous and the implicit utility y is not: it
*! solves, household by household, the cost identity
*!     y = ln x - sum_j w_j ln p_j + (1/2) sum_jk a_jk ln p_j ln p_k,
*! whose shares contain the errors (fixed point, a contraction here).
*! Log prices and log expenditure vary independently across households (sd
*! 0.3 and 0.6), so that every coefficient is identified with a t of 10 or
*! more at n = 3000; the mean shares are about 0.30, 0.25 and 0.45.
*!
*! censor: goods 1 and 2 are not bought by every household (Shonkwiler and
*! Yen 1999, the model of easi's option selection).  Household h buys good
*! i (i = 1, 2) when s_i'alpha_i + v_i > 0, with s_i = (1, ln p_1 - ln x,
*! ln p_2 - ln x, ln p_3 - ln x, z1, z2, q_i): q_i ~ N(0,1) enters the probit
*! of good i only (an exclusion restriction, selvars(w1: q1 ; w2: q2)).  Its
*! latent share is w*_i = f_i(y) + u_i with u_i = rho sd v_i + sqrt(1 -
*! rho^2) e_i, so that E[w_i] = Phi_i f_i + delta_i phi_i, delta_i = rho sd
*! (returned in r(delta1), r(delta2)); the observed share is d_i w*_i and
*! the last good closes the budget.  The implicit utility y is the one easi
*! defines: the cost identity on the latent shares, those of the
*! non-buyers replaced by their expectation f_i - delta_i phi_i/(1 - Phi_i)
*! (a fixed point, as without censoring).  The SY moment then holds up to a
*! term Cov(d_i, f_i(y)) ~ -b_y np_i delta_i phi_i of the order of 1e-4,
*! the price of any two-step censored system with an endogenous index; it
*! is immaterial for the comparison of standard errors.  The draws of the
*! uncensored process are unchanged (the censoring draws come after them).
*! Under normal errors a buyer's latent share can be negative, and easi
*! would read the household as a non-buyer: lift() (default 0.03) raises
*! the latent intercepts of goods 1 and 2 by that much, the last good's
*! falling by twice as much, which makes it rare -- 12 buyers in 200,000
*! at the defaults, counted in r(nneg) -- while the last good's share stays
*! positive (its minimum in r(w3min)).
capture program drop easi_dgp
program define easi_dgp, rclass
	version 14.2
	syntax , n(integer) seed(integer) [ POWer(integer 2) SDe(real 0.05)	///
		CENsor RHO(real 0.8) LIFT(real 0.03) ]
	if !inrange(`power', 1, 3) {
		di as error "easi_dgp: power() 1, 2 or 3"
		exit 198
	}
	if !inrange(`rho', -0.99, 0.99) {
		di as error "easi_dgp: rho() between -0.99 and 0.99"
		exit 198
	}
	clear
	quietly set obs `n'
	set seed `seed'
	quietly {
		gen double lx = rnormal(0, 0.6)
		forvalues j = 1/3 {
			gen double lp`j' = rnormal(0, 0.3)
		}
		gen double z1 = rnormal()
		gen byte   z2 = runiform() < 0.5
		gen double e1 = rnormal(0, `sde')
		gen double e2 = rnormal(0, `sde')
		gen double e3 = -e1 - e2
	}
	if "`censor'" == "" mata: _easi_dgp_solve(`power')
	else {
		quietly {
			gen double q1 = rnormal()
			gen double q2 = rnormal()
			gen double v1 = rnormal()
			gen double v2 = rnormal()
		}
		mata: _easi_dgp_censor(`power', `rho', `sde', `lift')
		drop v1 v2
		label variable q1 "probit of good 1 only"
		label variable q2 "probit of good 2 only"
		return scalar delta1 = `rho' * `sde'
		return scalar delta2 = `rho' * `sde'
		return scalar nneg = _dgp_nneg
		return scalar w3min = _dgp_w3min
		scalar drop _dgp_nneg _dgp_w3min
	}
	drop e1 e2 e3
	label variable w1 "budget share, good 1"
	label variable w2 "budget share, good 2"
	label variable w3 "budget share, good 3"
	label variable lx "log total expenditure"
	label variable z1 "demographic, continuous"
	label variable z2 "demographic, binary"
end

mata:
// the true coefficients, completed by adding-up (columns = goods)
void _easi_dgp_true(real rowvector b0, real matrix Br, real matrix G,
	real matrix A)
{
	b0 = (0.30, 0.25, 0.45)
	Br = ( 0.06, -0.02, -0.04 \			// y
	      -0.02,  0.015, 0.005 \			// y^2
	       0.004, -0.002, -0.002)			// y^3
	G  = ( 0.02, -0.01, -0.01 \			// z1
	      -0.03,  0.02,  0.01)			// z2
	A  = ( 0.08, -0.03, -0.05 \
	      -0.03,  0.06, -0.03 \
	      -0.05, -0.03,  0.08)
}

void _easi_dgp_solve(real scalar R)
{
	real matrix P, Z, E, Br, G, A, W, Yr
	real colvector lx, y, yold, pAp
	real rowvector b0
	real scalar it, r

	lx = st_data(., "lx")
	P  = st_data(., ("lp1", "lp2", "lp3"))
	Z  = st_data(., ("z1", "z2"))
	E  = st_data(., ("e1", "e2", "e3"))
	_easi_dgp_true(b0, Br, G, A)
	Br  = Br[| 1, 1 \ R, 3 |]
	pAp = rowsum((P * A) :* P)
	y   = lx
	for (it = 1; it <= 200; it++) {
		yold = y
		Yr = J(rows(y), R, .)
		for (r = 1; r <= R; r++) Yr[, r] = y:^r
		W = J(rows(y), 1, b0) + Yr * Br + Z * G + P * A + E
		y = lx - rowsum(W :* P) + 0.5 :* pAp
		if (max(abs(y - yold)) < 1e-14) break
	}
	if (max(abs(y - yold)) >= 1e-14) _error("easi_dgp: the fixed point did not converge")
	(void) st_addvar("double", ("w1", "w2", "w3"))
	st_store(., ("w1", "w2", "w3"), W)
}

// the true probit coefficients of goods 1 and 2 (columns): _cons, ln p_1 -
// ln x, ln p_2 - ln x, ln p_3 - ln x, z1, z2, q_i.  About 75% and 68% of
// buyers; the price terms sum to -0.2, so the richer buy more.
real matrix _easi_dgp_alpha()
{
	return(( 1.0,  0.7 \
	        -0.6,  0.2 \
	         0.2, -0.6 \
	         0.2,  0.2 \
	         0.3, -0.2 \
	        -0.3,  0.4 \
	         1.0,  1.0))
}

void _easi_dgp_censor(real scalar R, real scalar rho, real scalar sd,
	real scalar lift)
{
	real matrix P, Z, E, V, Qv, Br, G, A, AL, F, Yr, U, D, M0, Wt, S, W
	real colvector lx, y, yold, pAp, xb
	real rowvector b0
	real scalar it, r, i, dl

	lx = st_data(., "lx")
	P  = st_data(., ("lp1", "lp2", "lp3"))
	Z  = st_data(., ("z1", "z2"))
	E  = st_data(., ("e1", "e2"))
	Qv = st_data(., ("q1", "q2"))
	V  = st_data(., ("v1", "v2"))
	_easi_dgp_true(b0, Br, G, A)
	b0  = b0 + lift :* (1, 1, -2)
	Br  = Br[| 1, 1 \ R, 3 |]
	pAp = rowsum((P * A) :* P)
	AL  = _easi_dgp_alpha()
	dl  = rho * sd
	// buying, and the inverse Mills ratio of the non-buyers
	D  = J(rows(lx), 2, .)
	M0 = J(rows(lx), 2, .)
	for (i = 1; i <= 2; i++) {
		S  = (J(rows(lx), 1, 1), P :- lx, Z, Qv[., i])
		xb = S * AL[., i]
		D[., i]  = (xb + V[., i]) :> 0
		M0[., i] = normalden(xb) :/ normal(-xb)
	}
	// the latent errors, correlated with those of the probits
	U = rho * sd :* V + sqrt(1 - rho^2) :* E
	U = U, -rowsum(U)
	y = lx
	for (it = 1; it <= 200; it++) {
		yold = y
		Yr = J(rows(y), R, .)
		for (r = 1; r <= R; r++) Yr[, r] = y:^r
		F  = J(rows(y), 1, b0) + Yr * Br + Z * G + P * A
		// the completed latent shares of goods 1 and 2
		Wt = F[., 1..2] + D :* U[., 1..2] - (1 :- D) :* (dl :* M0)
		y  = lx - P[., 3] - rowsum(Wt :* (P[., 1..2] :- P[., 3])) + 0.5 :* pAp
		if (max(abs(y - yold)) < 1e-14) break
	}
	if (max(abs(y - yold)) >= 1e-14) _error("easi_dgp: the fixed point did not converge")
	W = D :* (F[., 1..2] + U[., 1..2])
	st_numscalar("_dgp_nneg", sum(D :& (W :<= 0)))
	W = W, 1 :- rowsum(W)
	st_numscalar("_dgp_w3min", min(W[., 3]))
	(void) st_addvar("double", ("w1", "w2", "w3"))
	st_store(., ("w1", "w2", "w3"), W)
}
end
