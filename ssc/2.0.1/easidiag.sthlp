{smcl}
{* *! version 2.0.1  30sep2026}{...}
{vieweralsosee "easi" "help easi"}{...}
{viewerjumpto "Syntax" "easidiag##syntax"}{...}
{viewerjumpto "Description" "easidiag##description"}{...}
{viewerjumpto "What it checks" "easidiag##checks"}{...}
{viewerjumpto "Examples" "easidiag##examples"}{...}
{viewerjumpto "Stored results" "easidiag##results"}{...}
{hline}
{hi:easidiag} {hline 2} Diagnose an EASI specification before estimating it
{hline}
{p 4 4 2}{txt}Package {cmd:easi}, version {res}2.0.1{txt} (30/09/2026) {c |} Stata {res}14.2{txt} or later {c |} first release {res}1.0.0{txt} (23/09/2026){p_end}

{marker syntax}{title:Syntax}

{p 8 15 2}
{cmd:easidiag} {it:sharevars} {ifin} {weight}{cmd:,}
{cmdab:pr:ices(}{it:varlist}{cmd:)}
{cmdab:exp:enditure(}{it:varname}{cmd:)}
{cmdab:demo:graphics(}{it:varlist}{cmd:)}
[{it:options}]

{p 4 4 2}
The syntax is that of {helpb easi}: the same share variables, the same
{opt prices()} or {opt lnprices()}, the same {opt expenditure()} or
{opt lnexpenditure()}, {opt demographics()}, {opt power()}, {opt py}, {opt pz},
{opt zy}, {opt interpz()}, {opt snames()} and weights.  Options that belong to
the estimator alone -- {opt vce()}, {opt dec()}, {opt compensated} and the rest
-- are accepted and ignored, and {cmd:easidiag} says which ones it dropped.

{p 4 4 2}
{opt hhsize()} and {opt elasticities()} are those of {helpb easi}, for the
weights of the estimation the diagnostic prepares: with {opt hhsize()} and
without {opt elasticities()}, or with {cmd:elasticities(individuals)}, each
household counts for its weight times {opt hhsize()}, as {cmd:easi} would
estimate it.

{p 4 4 2}
{opt pimpute()}, {opt selection}, {opt selgoods()} and {opt selvars()} are
those of {helpb easi} too, with the same rules: the missing prices are filled
as the estimation would fill them, and with {opt selection} the probits of the
estimation are run and section 9 reports them.


{marker description}{title:Description}

{p 4 4 2}
{cmd:easidiag} reports whether a specification is fit to be estimated, and
{bf:it does not estimate it}.  Everything in sections 1 to 6 is computed from
the data and from the design matrix evaluated at the starting point
{it:y} = the Stone index, which is a function of the data alone.  That is the
point: a diagnostic that needed the model to fit would be silent exactly when
the model does not fit.

{p 4 4 2}
The fixed-point section is the one exception, since the contraction factor can
only be measured by taking two steps.  It is reported last and is skipped when
the screening has already found something that stops the design being solved.

{p 4 4 2}
Two things make an EASI system fail to converge, or converge to something
meaningless: a regressor with almost no variance, and a near-linear dependency
among regressors.  A rank test finds neither reliably.  An exact dependency is
in fact the benign case -- the generalised inverse pins one column at zero and
the identified quantities come out right -- while a dependency that is exact to
one part in a billion leaves the matrix full rank and multiplies the standard
errors by 3e+08.  {cmd:easidiag} therefore measures condition indexes, which
are continuous, rather than testing rank.


{marker checks}{title:What it checks}

{dlgtab:1. Data admissibility}

{p 4 4 2}
Budget shares summing to one, shares outside [0,1], observations lost to
missing values.

{dlgtab:2. Screening of the design columns}

{p 4 4 2}
The standard deviation and the variance inflation factor of every column of the
design, including the constructed ones -- the powers of {it:y}, the normalised
log prices and all the interactions.  Then the condition indexes of Belsley,
Kuh and Welsch: a component whose index exceeds 30 and on which two or more
columns carry more than half their variance is a near-dependency, and the
variance-decomposition proportions {bf:name the columns involved}.  A group
that includes the constant is the classic exhaustive set of category dummies.

{dlgtab:3. Categorical and discrete variables}

{p 4 4 2}
For each demographic: the number of levels, the storage type, and for a binary
variable the size of the smaller cell.  A cell below 2% is flagged, because
under {opt pz} such a variable contributes {it:J}-1 price interactions that are
non-zero for a handful of households; merging categories is usually the answer.
A variable with a few integer levels is flagged as entering {it:linearly},
which imposes that a move from level 1 to 2 has the same effect as from 10 to
11.  A non-integer regressor stored as {cmd:float} is flagged because float
carries a relative error of 6e-08, sixty times what it takes to turn a harmless
exact dependency into an explosive one.

{dlgtab:4. Normalisation and conditioning}

{p 4 4 2}
The mean, the standard deviation and the ratio |mean|/sd of log expenditure and
of every log price, then the condition number of X'X as given and as it would
be if the variables were centred.  EASI is written for prices and expenditure
measured around a base period.  Fed raw levels, the columns 1, {it:y},
{it:y}{sup:2}, {it:y}{sup:3} correlate above 0.999 and the normal matrix goes
singular: on a Mexican survey, cond(X'X) was 4.5e+10 in levels and 4.5e+01
centred.

{dlgtab:5. Structural sources of a flat objective}

{p 4 4 2}
{bf:The polynomial basis.} The conditioning of the Gram matrix of
1, {it:y}, ..., {it:y}{sup:R}, reported for every power up to {opt power()}.
Centred, each extra power costs a factor of about four; uncentred, about a
hundred.  The number at the chosen power tells you how far you can push
{opt power()} on your own data.

{p 4 4 2}
{bf:Distinct price vectors.} When prices are regional or temporal indexes,
every household in a cell shares one price vector and the price parameters are
identified off the number of cells, not the sample size.

{p 4 4 2}
{bf:Rank of the price variation}, the {bf:first-stage R{sup:2}} for {it:y}, and
the conditioning of the covariance of the estimated shares, which stands in for
the residual covariance before there are any residuals.

{dlgtab:6. Scale and variability}

{p 4 4 2}
Mean and standard deviation of each budget share, and the standard deviation of
each normalised log price.  A mean share below 0.01 is flagged: every
semi-elasticity of that good is divided by it.

{dlgtab:7. Convergence of the fixed point}

{p 4 4 2}
Phi = 1 - p'Bp/2 is the denominator of {it:y}; as it approaches zero the map
stops contracting, and a negative Phi means the cost function has left its
regular region.  The contraction factor is measured by running the iteration
twice from the two-stage least squares start and taking
|y2-y1|/|y1-y0|; if it is below one, the number of iterations needed follows
from log(tol/|y1-y0|)/log(ratio).  On the reference data that prediction is
exact.

{dlgtab:8. Identification}

{p 4 4 2}
Ranks of Z'WZ and of Xhat'W Xhat, and the count of instruments against
coefficients per equation.

{dlgtab:9. Selection of the buyers}

{p 4 4 2}
With {opt selection}: for each corrected good, the probit that {helpb easi}
would estimate (same sample, same weights), with the percentage of buyers, its
pseudo-R2 (McFadden), the households predicted with probability 0 or 1
(separation), the number of variables of the probit only ({opt selvars()}), and
{it:VIF(delta)} at the starting point {it:y} = the Stone index.  In the
expected share Phi {it:f} + {it:delta} phi, {it:delta} is identified through
the variables of the probit only or, without them, through the nonlinearity of
the normal distribution alone; {it:VIF(delta)} = 1/(1 - {it:R}^2), {it:R}^2 the
uncentred R-squared of phi on the other regressors of the equation (those of
EASI times Phi), measures how far phi is from a combination of them.  Above 10
the good is flagged: {it:delta} rests on little more than the curvature of the
probit; add a variable of the probit only in {opt selvars()}, or leave the good
uncorrected with {opt selgoods()}.  {cmd:easi} reports the same measure at the
converged {it:y}.  On the Mexican cereals without a variable of the probit only,
it is 33 for corn and about 2,000 for wheat.


{marker examples}{title:Examples}

{p 4 4 2}Before estimating:{p_end}
{phang2}{cmd:. easidiag w1-w9, prices(p1-p9) expenditure(totexp) demographics(age hsex carown) power(5)}{p_end}

{p 4 4 2}The same specification, with a survey weight:{p_end}
{phang2}{cmd:. easidiag w1-w9 [pw=sweight], prices(p1-p9) expenditure(totexp) demographics(age hsex) power(3)}{p_end}

{p 4 4 2}On the Mexican cereals that come with the package (an ancillary file;
see the examples of {helpb easi}), raw and centred logarithms (the example 7
of {helpb easi}):{p_end}
{phang2}{cmd:. use mex_bench, clear}{p_end}
{phang2}{cmd:. easidiag w1 w2 w3, lnprices(lp1_raw lp2_raw lp3) lnexpenditure(lx_raw) demographics(z1 z2) power(3)}{p_end}
{phang2}{cmd:. easidiag w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3)}{p_end}
{p 8 8 2}{txt}({stata "easi_examples 7":click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 7, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 7, do":open as a do-file}){p_end}

{p 4 4 2}
{cmd:easidiag} is also reached from the dialog box: {cmd:db easi}, then
{bf:Action: 2 - Diagnostic}.


{marker results}{title:Stored results}

{pstd}{cmd:easidiag} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}observations used{p_end}
{synopt:{cmd:r(nprob)}}number of points flagged{p_end}
{synopt:{cmd:r(cond_raw)}}cond(X'X) as given{p_end}
{synopt:{cmd:r(cond_ctr)}}cond(X'X) if centred{p_end}
{synopt:{cmd:r(ratio)}}contraction factor of the fixed point{p_end}
{synopt:{cmd:r(iter)}}iterations predicted{p_end}
{synopt:{cmd:r(phimin)}}smallest Phi{p_end}
{synopt:{cmd:r(wmin)}}smallest mean budget share{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(selgoods)}}the goods corrected ({opt selection}){p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(var)}}mean, sd and |mean|/sd of expenditure and the prices{p_end}
{synopt:{cmd:r(good)}}mean share, sd of the share, sd of the normalised price{p_end}
{synopt:{cmd:r(sel_diag)}}by good: percentage of buyers, pseudo-R2, perfectly predicted, VIF(delta) at the start ({opt selection}){p_end}


{title:Author}

{pstd}Abdelkrim Araar, Universite Laval / PEP{break}
{browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}{p_end}

{pstd}Version 2.0.1, 30 September 2026 (first release 1.0.0, 23 September 2026).  License: GPL-3.0-or-later.
{browse "https://github.com/aabbdd12/easi"}{p_end}


{title:Also see}

{psee}
{helpb easi}, {helpb easi##predict:predict after easi},
{helpb easi##engel:estat engel}
{p_end}
