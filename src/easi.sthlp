{smcl}
{* 21sep2026}{...}
{vieweralsosee "[R] nlsur" "help nlsur"}{...}
{viewerjumpto "Syntax" "easi##syntax"}{...}
{viewerjumpto "Description" "easi##description"}{...}
{viewerjumpto "Options" "easi##options"}{...}
{viewerjumpto "Remarks" "easi##remarks"}{...}
{viewerjumpto "Differences from sr_easi" "easi##compat"}{...}
{viewerjumpto "Examples" "easi##examples"}{...}
{viewerjumpto "Stored results" "easi##results"}{...}
{hline}
{hi:easi} {hline 2} Exact Affine Stone Index (EASI) demand system
{hline}

{marker syntax}{title:Syntax}

{p 8 15 2}
{cmd:easi} {it:sharevars} {ifin} {weight}{cmd:,}
{cmdab:pr:ices(}{it:varlist}{cmd:)}
{cmdab:exp:enditure(}{it:varname}{cmd:)}
{cmdab:demo:graphics(}{it:varlist}{cmd:)}
[{it:options}]

{p 4 4 2}
{it:sharevars} is the list of budget-share variables, one per good.  They must
sum to 1 in every observation.  The {bf:last} good is the one dropped from the
system and recovered by adding up.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt pr:ices(varlist)}}prices, in levels; logged internally{p_end}
{synopt:{opt lnpr:ices(varlist)}}prices, already in logarithms{p_end}
{synopt:{opt exp:enditure(varname)}}total expenditure, in levels{p_end}
{synopt:{opt lnexp:enditure(varname)}}total expenditure, already logged{p_end}
{synopt:{opt demo:graphics(varlist)}}demographic variables; at least one{p_end}
{synopt:{opt pow:er(#)}}highest power of the implicit utility index; default
{cmd:power(5)}{p_end}

{syntab:Interactions}
{synopt:{opt py}}interact prices with the implicit utility index{p_end}
{synopt:{opt pz}}interact prices with demographics{p_end}
{synopt:{opt zy}}interact demographics with the implicit utility index{p_end}
{synopt:{opt interpz(varlist)}}which demographics interact with prices; default
is all of them{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}{opt r:obust} (default), {opt cl:uster}
{it:clustvar}, {opt svy} or {opt conv:entional}{p_end}

{syntab:Reporting}
{synopt:{opt sn:ames(namelist)}}short names of the goods, used in the
elasticity tables{p_end}
{synopt:{opt dec(#)}}decimals in the displayed tables; default {cmd:dec(4)}{p_end}
{synopt:{opt dislas(#)}}{cmd:dislas(0)} hides the last good; default
{cmd:dislas(1)}{p_end}
{synopt:{opt dregres(#)}}{cmd:dregres(1)} also shows the coefficient table{p_end}
{synopt:{opt compensat:ed}}add the compensated (Hicksian) price elasticity
table{p_end}
{synopt:{opt demoel:ast}}add the demographic elasticity table{p_end}
{synopt:{opt checks}}show the aggregation identities{p_end}
{synopt:{opt det:ail}}all three of the above{p_end}
{synopt:{opt noelastse}}skip the elasticity standard errors and their tables{p_end}
{synopt:{opt nolog}}suppress the iteration log{p_end}

{syntab:Advanced}
{synopt:{opt compat}}reproduce the R package {bf:easi} 0.21 bit for bit,
including its bugs{p_end}
{synopt:{opt tol:erance(#)}}convergence tolerance; default {cmd:tolerance(1e-6)}{p_end}
{synopt:{opt iter:ate(#)}}maximum iterations; default {cmd:iterate(100)}{p_end}
{synoptline}
{p 4 6 2}{opt aweight}s, {opt fweight}s, {opt pweight}s and {opt iweight}s are
allowed; see {help weight}.  {opt pweight}s imply {cmd:vce(robust)}.{p_end}

{p 4 6 2}
{cmd:easi} is {help estcom:e-class}; {cmd:predict} is available, see
{help easi postestimation##predict:below}.{p_end}


{marker description}{title:Description}

{p 4 4 2}
{cmd:easi} estimates the Exact Affine Stone Index demand system of Lewbel and
Pendakur (2009).  Budget shares are linear in the parameters conditional on a
measure of real expenditure -- the {it:implicit utility index} {it:y} -- which
itself depends on the parameters.  The system is estimated by iterated
three-stage least squares with Slutsky symmetry imposed as cross-equation
restrictions, iterating on {it:y} until it stops moving.

{p 4 4 2}
Everything is computed in Stata and Mata.  Earlier versions of this module
({helpb sr_easi}) wrote an R script and shelled out to R; that is no longer the
case, and R need not be installed.

{p 4 4 2}
The command reports expenditure elasticities and uncompensated price
elasticities with their standard errors, and stores demographic elasticities,
semi-elasticities, the Slutsky matrix and compensated quantity derivatives in
{cmd:e()}.


{marker options}{title:Options}

{dlgtab:Model}

{phang}
{opt prices(varlist)} / {opt lnprices(varlist)} give the price of each good, in
the same order as {it:sharevars}.  Use {cmd:prices()} for prices in levels (they
are logged internally) or {cmd:lnprices()} if they are already logarithms.
Exactly one of the two is required.

{phang}
{opt expenditure(varname)} / {opt lnexpenditure(varname)} give total household
expenditure, in levels or already logged.  Exactly one is required.

{phang}
{opt demographics(varlist)} lists the demographic variables.  At least one is
required.

{phang}
{opt power(#)} sets the highest power of {it:y} in the Engel curves.  The
default is 5.  Higher powers allow more flexible Engel curves at the cost of
parameters.

{dlgtab:Interactions}

{phang}
{opt py}, {opt pz} and {opt zy} add, respectively, price x expenditure, price x
demographic and demographic x expenditure interactions.  The legacy spellings
{cmd:inpy(1)}, {cmd:inpz(1)} and {cmd:inzy(1)} are also accepted.

{phang}
{opt interpz(varlist)} restricts the price x demographic interactions to the
listed demographics.  The default is all of them.  {bf:This option behaves
differently from the R package and from} {helpb sr_easi} -- see
{help easi##compat:Differences} below.

{marker vce}{dlgtab:SE/Robust}

{phang}
{opt vce(vcetype)} selects the variance estimator.  {cmd:robust}, a sandwich
estimator, is the default; {cmd:cluster} {it:clustvar} allows for correlation
within clusters; {cmd:svy} reads the design declared by {helpb svyset};
{cmd:conventional} is the homoskedastic 3SLS covariance matrix that
{helpb reg3}, {cmd:systemfit} and the R package report, and is what
{opt compat} uses.

{phang}
{bf:If your data come from a clustered survey, use {cmd:vce(svy)} or at least}
{bf:{cmd:vce(cluster} {it:psu}{cmd:)}.}  Against a design bootstrap -- PSUs
resampled within strata, weights travelling with the rows -- {cmd:vce(robust)}
understates the standard errors by 7.5% on average and by 36% at worst, and the
coefficients it misses most are those on {it:y}.  Both design-aware estimators
close that gap to within 2%.  This is the largest error we have measured in the
command, and it is silent.

{phang}
{cmd:vce(svy)} takes the sampling weight, the primary sampling unit, the strata
and the finite-population correction from the {cmd:svyset} characteristics.  The
weight is adopted automatically when none is given on the command line.  The
estimator is the stratified ultimate-cluster linearization: scores are summed
within each PSU, the PSU totals are centred on {it:their own stratum mean}, and
each stratum contributes {it:n_h}/({it:n_h}-1) times the cross-product of those
deviations, times 1-{it:f_h} when an {opt fpc()} is declared.  Only the
first-stage correction enters -- that is what makes the ultimate-cluster
approximation legitimate.  Without an {opt fpc()} the PSUs are treated as drawn
with replacement, which overstates the variance: safe, but not exact when the
sampling fraction is appreciable.  A stratum holding a single PSU contributes
nothing and the command says so rather than letting the variance quietly
shrink.

{phang}
The default is {cmd:robust} because budget-share equations are heteroskedastic
and the cost of ignoring it falls on the coefficients that matter most.  Against
a 400-replication bootstrap of the whole iterated procedure, the conventional
standard errors of the coefficients on {it:y} -- the ones that shape the Engel
curves -- are 23% too small on average and up to 72% too small; the robust ones
are within 1% on average and 9% at worst.

{phang}
A second, smaller correction is applied by default.  The design contains
{it:y}, which is itself a function of the coefficients, so the 3SLS covariance
is the covariance of the last linear step {it:conditional} on {it:y}; Pendakur's
code flags this in a comment.  {cmd:easi} adds the missing
{it:dy}/{it:d}{bf:b} term to the Jacobian of the estimating equations.  It is
worth less than half a percent here, because {it:y} depends on the coefficients
only through p'A(z)p/2 and p'Bp/2 and normalised log prices are small, but it is
the right Jacobian.  {opt compat} restores the conditional covariance.

{dlgtab:Reporting}

{phang}
{opt snames(namelist)} supplies short labels for the goods, one per share
variable, used to head the elasticity tables.

{phang}
{opt dislas(0)} drops the last good from the displayed tables.  Its elasticities
are still stored in {cmd:e()}.

{phang}
{opt compensated} adds Table 04, the compensated (Hicksian) price elasticities,
and computes their standard errors.  They are off by default because the output
is already long and not everyone wants them; leaving them off also skips their
cost, since the standard errors are only computed for what is shown.

{phang}
{opt demoelast} adds Table 05, the demographic elasticities.  They are always in
{cmd:e(elast_demo)}.

{phang}
{opt checks} shows the aggregation identities.  They are {it:always computed},
and a failure is reported whether or not you ask for it: Engel aggregation
(sum of w_j times the expenditure elasticities equals one) and Cournot
aggregation (the share-weighted sum of a price column equals minus that share)
hold exactly for these formulas, so a departure means the elasticity code is
wrong on your data, which you should know about either way.  Residuals of 1e-16
to 1e-10 are normal.

{phang}
{opt detail} is {opt compensated}, {opt demoelast} and {opt checks} together.

{phang}
All four can also be given when replaying results, as in
{cmd:. easi, compensated}.

{dlgtab:Advanced}

{phang}
{opt compat} makes {cmd:easi} reproduce the R package {bf:easi} 0.21 exactly,
including the defects listed under {help easi##compat:Differences}.  Use it to
reproduce results published with {helpb sr_easi}.  Do not use it for new work.


{marker remarks}{title:Remarks}

{dlgtab:Reading the price elasticity tables}

{p 4 4 2}
Tables 03 and 04 are indexed {bf:row = good, column = price}: the entry in row
{it:j} and column {it:k} is the elasticity of the quantity of good {it:j} with
respect to the price of good {it:k}.  Own-price elasticities are on the
diagonal.  {cmd:e(elast_price_nc)} and {cmd:e(elast_price_c)} follow the same
convention.

{p 4 4 2}
On Stata 17 and later the elasticity matrices are displayed with {helpb table},
which ignores {help linesize} and prints the whole J x J block in one piece with
the good names intact, however many goods there are; the Results window scrolls
it sideways.  On Stata 14 to 16 the display falls back to {helpb matlist}, which
wraps at {cmd:c(linesize)}: with more than about six goods the matrix is split
into blocks and the names are abbreviated.  Widening the Results window, or
{cmd:set linesize}, reduces the splitting.  Everything else is identical.

{p 4 4 2}
The legacy matrices do not.  {cmd:e(elast_price)}, kept for backward
compatibility, is indexed {bf:[price, good]} -- the transpose -- and
{cmd:e(compensated_q)} is indexed [good, price] but holds the compensated
elasticities {bf:plus the identity matrix}, so its diagonal is every own-price
elasticity plus one.  Both come from the R package.  Because both axes carry
the same good names, nothing in a printed table signals which way to read it,
which is why the current names exist and the display states the orientation.

{dlgtab:Standard errors of the elasticities}

{p 4 4 2}
An elasticity is a nonlinear function of the coefficients averaged over
households.  {cmd:easi} applies the delta method to the quantity it actually
reports: se = sqrt(g'Vg) with g the gradient of that average elasticity.  The
gradient of the whole elasticity vector is taken at once, so the dropped good is
treated like any other -- its elasticity is a function of the coefficients too,
and the covariances between equations enter as they should.

{p 4 4 2}
An expenditure elasticity is a {it:ratio} of two estimators,
{it:eta_j} = 1 + {it:alephbar_j}/{it:wbar_j}, the denominator being a mean
budget share estimated from the same sample.  Differentiating with respect to
the coefficients alone, with {it:wbar} held fixed, keeps one term of three.  The
missing pieces are not a forgotten variance: {it:alephbar} and {it:wbar} move
together, so the cross term is negative and the omission {it:overstates} the
standard error -- on a survey with mean shares near 0.02, by 18%.  {cmd:easi}
uses the full delta method, in closed form, and the result then agrees with a
bootstrap to within 3%.  The same influence function carries the survey design,
so this works unchanged under {cmd:vce(svy)}.

{p 4 4 2}
The R package instead reports, for each elasticity, the median across households
of the pointwise standard error of the corresponding {it:semi}-elasticity divided
by the mean budget share.  That is not the standard error of anything it
reports, and against a bootstrap it is out by a factor of fourteen on the
dropped good.  {cmd:legacy(elastse)} restores it; {opt compat} implies it.

{p 4 4 2}
{opt noelastse} skips the elasticity standard errors.  They used to be the
expensive part of the command: the Jacobian of the elasticities was taken by
finite differences, one full re-evaluation per coefficient, and on the reference
example that was two thirds of the run time.  The Jacobian is now analytic and
the standard errors cost nothing measurable (8.6 seconds with them, 8.8
without, on the reference example), so {opt noelastse} is a choice about the
output rather than about speed.  When they are computed, every elasticity table
is followed by a table of its standard errors; with {opt noelastse} those tables,
and the standard-error column of Table 02, are simply absent.

{p 4 4 2}
{cmd:easi} reports its own execution time and stores it in {cmd:e(time)}.  It
measures this with Stata timer number 100.  Stata's timers are a single global
pool and a running timer cannot be read back, so there is no way to save one and
put it back: if you are using timer 100 yourself, {cmd:easi} will overwrite it.

{dlgtab:Weights}

{p 4 4 2}
{opt aweight}s and {opt pweight}s are normalised to sum to the number of
observations; {opt fweight}s and {opt iweight}s are used as they stand.  All the
population means used to form the elasticities are weighted, as is the median
above -- which makes {opt fweight}s exactly equivalent to duplicating
observations.  The R package has no weights at all.


{marker compat}{title:Differences from sr_easi and the R package}

{p 4 4 2}
{cmd:easi} fixes six defects of the R package that {helpb sr_easi} inherited.
{cmd:compat} switches them all back on.  Only the first two change results
materially.

{phang}
{bf:1. interpz.}  In the R package the line
{cmd:interpz <- ifelse(length(interpz) > 1, interpz, 1:nsoc)} silently collapses
the vector to its {bf:first element}, because {cmd:ifelse()} returns a result
the length of its test.  So {cmd:pz} only ever interacted prices with the
{bf:first} demographic variable, whatever was requested, and
{cmd:interpz(}{it:third variable}{cmd:)} selected the first one.  {cmd:easi}
interacts prices with every variable in {opt interpz()}.  {bf:This changes
results substantially}: in the reference example the system grows from 320 to
576 coefficients and price elasticities move by up to 2.0 in absolute value.

{phang}
{bf:2. Instruments.}  The R estimation loop restarts from the untouched data and
never refreshes its instrument columns, so the instruments are frozen at the
mean-share Stone index.  {cmd:easi} refreshes them each iteration.  Effect on
the reference example: about 4e-04 on the expenditure elasticities.

{phang}
{bf:3. Quantisation.}  The R code rounds the price quadratic forms to 1e-6.
{cmd:easi} does not.  Effect: about 2e-07.

{phang}
{bf:4. Demographic elasticities.}  The R accumulator is initialised outside the
loop over demographics, so effects accumulate from one variable to the next, and
it is indexed one position short, so the constant block is read instead of the
first demographic.  Effect: about 7e-03, on the demographic elasticities only.

{phang}
{bf:5. Standard error of the last good.}  The R code returns it with a minus
sign, and omits the squares in two closing cells.  Effect: on those cells only.

{phang}
{bf:6. Elasticity formulas.}  The R derivation differentiates the implicit
Marshallian shares holding the shares themselves fixed, so it misses that
{it:y} responds to prices and to total expenditure.  {cmd:easi} uses the correct
derivatives, which are also simpler; they agree with finite differences to nine
significant digits and satisfy degree-zero homogeneity to 6e-09, against 6e-03
for the R formulas.  {bf:This changes reported elasticities.}

{p 4 4 2}
Four further changes are improvements rather than fixes, and are {bf:not}
reverted by {opt compat} except where stated.

{phang}
{bf:a. Iterated residual covariance.}  {cmd:systemfit} and {cmd:reg3} form Sigma
once and take one GLS step.  For a singular demand system the estimator is
invariant to the deleted equation only if Sigma is iterated to convergence
(Barten 1969); with one step the price elasticities move by up to 0.20 depending
on which good is dropped.  {cmd:easi} iterates.  {cmd:legacy(sigma1)}, implied by
{opt compat}, takes the single step.

{phang}
{bf:b. Robust standard errors by default}, and {opt vce(cluster)}.  See
{help easi##vce:SE/Robust} above.  {opt compat} uses {cmd:vce(conventional)}.

{phang}
{bf:c. Generated regressor.}  The Jacobian of the estimating equations includes
the dependence of {it:y} on the coefficients.  {opt compat} reverts it.

{phang}
{bf:d. Weights.}  {cmd:easi} supports {opt aweight}s, {opt fweight}s,
{opt pweight}s and {opt iweight}s.  The R package has none.


{marker examples}{title:Examples}

{p 4 4 2}A basic system with five demographics:{p_end}
{phang2}{cmd:. easi w1-w9, prices(p1-p9) expenditure(totexp) demographics(age hsex carown)}{p_end}

{p 4 4 2}Prices and expenditure already in logs, all interactions, power 5:{p_end}
{phang2}{cmd:. easi w1-w9, lnprices(lp1-lp9) lnexpenditure(ly) demographics(age hsex) power(5) py pz zy}{p_end}

{p 4 4 2}Survey data, clustered standard errors, only two demographics crossed
with prices:{p_end}
{phang2}{cmd:. easi w1-w5 [pw=sweight], prices(p1-p5) expenditure(x) demographics(age sex educ) pz interpz(age sex) vce(cluster psu)}{p_end}

{p 4 4 2}Reproduce an old {helpb sr_easi} result exactly:{p_end}
{phang2}{cmd:. easi w1-w5, prices(p1-p5) expenditure(x) demographics(age sex) power(3) py compat}{p_end}

{p 4 4 2}Postestimation:{p_end}
{phang2}{cmd:. predict double shat*, shares}{p_end}
{phang2}{cmd:. predict double y, y}{p_end}
{phang2}{cmd:. matrix list e(elast_price_nc)}{p_end}
{phang2}{cmd:. easi, dec(6)}{p_end}


{marker results}{title:Stored results}

{pstd}{cmd:easi} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_eff)}}effective sample size used as the divisor of Sigma{p_end}
{synopt:{cmd:e(ngoods)}}number of goods{p_end}
{synopt:{cmd:e(neq)}}number of estimated equations, {cmd:e(ngoods)}-1{p_end}
{synopt:{cmd:e(nsoc)}}number of demographics{p_end}
{synopt:{cmd:e(power)}}highest power of {it:y}{p_end}
{synopt:{cmd:e(k_eq)}}coefficients per equation{p_end}
{synopt:{cmd:e(py)}, {cmd:e(pz)}, {cmd:e(zy)}}interaction flags{p_end}
{synopt:{cmd:e(iter)}}iterations used{p_end}
{synopt:{cmd:e(crit)}}final convergence criterion{p_end}
{synopt:{cmd:e(converged)}}1 if converged{p_end}
{synopt:{cmd:e(time)}}execution time in seconds{p_end}
{synopt:{cmd:e(N_psu)}, {cmd:e(N_strata)}}PSUs and strata, with {cmd:vce(svy)}{p_end}
{synopt:{cmd:e(chk_engel)}}|sum_j w_j eta^x_j - 1|, Engel aggregation{p_end}
{synopt:{cmd:e(chk_cournot)}}max_k |sum_j w_j eta^k_j + w_k|, Cournot
aggregation{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:easi}{p_end}
{synopt:{cmd:e(shares)}, {cmd:e(prices)}, {cmd:e(demographics)}}variable lists{p_end}
{synopt:{cmd:e(interpz)}}demographics crossed with prices{p_end}
{synopt:{cmd:e(vce)}, {cmd:e(clustvar)}}variance estimator{p_end}
{synopt:{cmd:e(svyunit)}, {cmd:e(svystrata)}}PSU and strata, with {cmd:vce(svy)}{p_end}
{synopt:{cmd:e(svywvar)}, {cmd:e(svyfpc)}}survey weight and fpc variables{p_end}
{synopt:{cmd:e(vcetype2)}}{cmd:generated regressor} or {cmd:conditional on y}{p_end}
{synopt:{cmd:e(wtype)}, {cmd:e(wexp)}}weight type and expression{p_end}
{synopt:{cmd:e(mode)}}{cmd:corrected}, {cmd:compat} or {cmd:mixed}{p_end}
{synopt:{cmd:e(report)}}reporting options in force{p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}coefficients and their variance matrix{p_end}
{synopt:{cmd:e(elast_exp)}, {cmd:e(elast_exp_se)}}expenditure elasticities{p_end}
{synopt:{cmd:e(elast_price_nc)}, {cmd:e(elast_price_nc_se)}}uncompensated
(Marshallian) price elasticities, [good, price]{p_end}
{synopt:{cmd:e(elast_price_c)}, {cmd:e(elast_price_c_se)}}compensated
(Hicksian) price elasticities, [good, price]; the standard errors are missing
unless {opt compensated} was specified{p_end}
{synopt:{cmd:e(elast_demo)}, {cmd:e(elast_demo_se)}}demographic elasticities,
[demographic, good]{p_end}
{synopt:{cmd:e(semi_exp)}, {cmd:e(semi_price)}}semi-elasticities{p_end}
{synopt:{cmd:e(slutsky)}}Slutsky matrix in share form{p_end}
{synopt:{cmd:e(Sigma)}}cross-equation residual covariance{p_end}

{p2col 5 24 28 2: Matrices, legacy}{p_end}
{p 4 6 2}Kept unchanged so that existing do-files keep working.  See
{help easi##remarks:Reading the price elasticity tables}.{p_end}
{synopt:{cmd:e(elast_income)}, {cmd:e(elast_income_se)}}same as
{cmd:e(elast_exp)}{p_end}
{synopt:{cmd:e(elast_price)}, {cmd:e(elast_price_se)}}uncompensated price
elasticities, [price, good] -- the {bf:transpose} of {cmd:e(elast_price_nc)}{p_end}
{synopt:{cmd:e(semi_income)}}same as {cmd:e(semi_exp)}{p_end}
{synopt:{cmd:e(compensated_q)}}{cmd:e(elast_price_c)} {bf:plus the identity
matrix}{p_end}


{marker predict}{title:Postestimation: predict}

{p 8 15 2}
{cmd:predict} [{it:type}] {it:stub}{cmd:*} {ifin} [{cmd:,} {opt sh:ares} | {opt res:iduals}]{p_end}
{p 8 15 2}
{cmd:predict} [{it:type}] {it:newvar} {ifin}{cmd:,} {opt y}{p_end}

{phang}{opt shares} (the default) stores one fitted budget share per good.  The
last one is obtained by adding up, so the fitted shares sum to 1 exactly.{p_end}

{phang}{opt residuals} stores observed minus fitted shares, one per good.{p_end}

{phang}{opt y} stores the implicit utility index.{p_end}

{pstd}Everything is recomputed from {cmd:e(b)}, so predictions are available
outside the estimation sample.


{marker engel}{title:Postestimation: estat engel}

{p 8 15 2}
{cmd:estat engel} {ifin} [{cmd:,} {it:options}]{p_end}

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt atm:eans}}demographics and prices at their weighted means; the
default{p_end}
{synopt:{opt asob:served}}demographics and prices as observed{p_end}
{synopt:{opt n(#)}}number of points at which the curve is evaluated; default
{cmd:n(100)}{p_end}
{synopt:{opt bw:idth(#)}}smoothing bandwidth; {opt asobserved} only{p_end}
{synopt:{opt trim(#)}}percent trimmed from each tail of the grid; default
{cmd:trim(1)}{p_end}
{synopt:{opt l:evel(#)}}confidence level for the band{p_end}
{synopt:{opt noci}}omit the confidence band{p_end}
{synopt:{opt data(filename)}}save the plotted curves as a dataset{p_end}
{synopt:{opt sav:ing(filename)}}save the graph{p_end}
{synopt:{opt nodraw}}compute but do not draw{p_end}
{synoptline}

{pstd}
{cmd:estat engel} traces the fitted budget share of each good against total
expenditure. Two objects are available, and only one of them involves
smoothing.

{dlgtab:atmeans -- the Engel curve}

{pstd}
The default. The model's fitted share is evaluated over a grid of total
expenditure with the demographics and the prices held at their weighted means.
At each grid point the fixed point in the budget shares and the implicit utility
index is solved exactly, so this is the {bf:exact function}, not an estimate of
it: it is smooth by construction and nothing is fitted to it. {opt n()} is
therefore a {bf:resolution}, not a smoothing parameter, and {opt bwidth()} is
refused -- there is nothing to smooth. This isolates the expenditure effect,
which is what an Engel curve is.

{dlgtab:asobserved -- the sample profile}

{pstd}
The fitted shares as they actually vary with expenditure, leaving the
demographics and prices as observed. This also picks up however those covary
with expenditure, so it is a scatter and does need smoothing. {cmd:easi} uses a
{bf:local linear smoother} rather than binning: bins are a rectangular kernel
with arbitrary breakpoints, and a bandwidth is the better-behaved choice.

{pstd}
{opt bwidth()} sets the bandwidth. Left out, it is the rule of thumb that
{helpb lpoly} computes for local polynomial {bf:regression}; note that
Silverman's rule is a rule for {bf:density} estimation and is not the
appropriate one for a smoother. One bandwidth is used for every panel so that
they are comparable, and the value is reported in the figure note and in
{cmd:r(bwidth)}.

{dlgtab:Both}

{pstd}
The confidence band is the delta-method standard error of the fitted share,
{it:x'Vx}, at each evaluation point. For the good dropped from the system it
uses the full sum of the covariance blocks, so covariances across equations are
accounted for. The band treats the implicit utility index as fixed.

{pstd}
A polynomial of order {opt power()} oscillates near the ends of the range of the
implicit utility index, which can make the extreme percentiles wild. That is the
model, not an error, but it dominates the vertical scale, so {opt trim()} drops
a little from each tail of the grid by default; {cmd:trim(0)} shows everything.

{pstd}
The panels are laid out on a near-square grid whatever the number of goods --
2 x 2 for four goods, 3 x 3 for nine, 4 x 3 for twelve -- and the canvas grows
with the grid so that every panel keeps the same size.  Any {it:twoway_options}
given to {cmd:estat engel} are passed to {helpb graph combine}; {opt cols()},
{opt rows()}, {opt xsize()} and {opt ysize()} override the automatic layout.

{pstd}
{cmd:estat engel} stores {cmd:r(n)} and, with {opt asobserved},
{cmd:r(bwidth)}.

{phang2}{cmd:. estat engel}{p_end}
{phang2}{cmd:. estat engel, n(200) level(90) saving(engel, replace)}{p_end}
{phang2}{cmd:. estat engel, asobserved}{p_end}
{phang2}{cmd:. estat engel, asobserved bwidth(.15) trim(2)}{p_end}


{title:References}

{p 4 8 2}Lewbel, A., and K. Pendakur. 2009.
Tricks with Hicks: The EASI demand system.
{it:American Economic Review} 99: 827-863.{p_end}

{p 4 8 2}Pendakur, K. 2008. EASI made Easier.
{browse "http://www.sfu.ca/~pendakur/"}{p_end}

{p 4 8 2}Hoareau, S., G. Lacroix, M. Hoareau, and L. Tiberti. 2012.
Exact Affine Stone Index Demand System in R: The easi Package.
CIRPEE working paper.{p_end}


{title:Author}

{pstd}Abdelkrim Araar{break}
{browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}{p_end}


{title:Also see}

{psee}{helpb sr_easi} (compatibility wrapper for the old syntax){p_end}
