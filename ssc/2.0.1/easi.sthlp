{smcl}
{* *! version 2.0.1  30sep2026}{...}
{vieweralsosee "[R] nlsur" "help nlsur"}{...}
{viewerjumpto "Syntax" "easi##syntax"}{...}
{viewerjumpto "Description" "easi##description"}{...}
{viewerjumpto "Options" "easi##options"}{...}
{viewerjumpto "Remarks" "easi##remarks"}{...}
{viewerjumpto "Which elasticities" "easi##types"}{...}
{viewerjumpto "The selection of the buyers" "easi##selection"}{...}
{viewerjumpto "Differences from sr_easi" "easi##compat"}{...}
{viewerjumpto "Examples" "easi##examples"}{...}
{viewerjumpto "Stored results" "easi##results"}{...}
{hline}
{hi:easi} {hline 2} Exact Affine Stone Index (EASI) demand system
{hline}
{p 4 4 2}{txt}Package {cmd:easi}, version {res}2.0.1{txt} (30/09/2026) {c |} Stata {res}14.2{txt} or later {c |} first release {res}1.0.0{txt} (23/09/2026){p_end}

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
is all of them; implies {opt pz}{p_end}

{syntab:Missing prices and non-buyers}
{synopt:{opt pim:pute(varlist)}}fill a missing price by the mean log price of the same group, the groups tried in order{p_end}
{synopt:{opt sel:ection}}correct for the households that do not buy (Shonkwiler and Yen 1999){p_end}
{synopt:{opt selg:oods(namelist)}}the shares corrected; default: those with zeros, the last good excepted; implies {opt selection}{p_end}
{synopt:{cmdab:selv:ars(}[{it:good}{cmd::}] {it:varlist} [{cmd:;} ...]{cmd:)}}variables of the probits only,
for every corrected good or for one good; implies {opt selection}{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}{opt r:obust} (default), {opt cl:uster}
{it:clustvar}, {opt svy} or {opt conv:entional}{p_end}
{synopt:{cmd:vce(bootstrap} [{cmd:,} {it:{help easi##bootopts:boot_opts}}]{cmd:)}}bootstrap of the whole procedure{p_end}

{syntab:Elasticities}
{synopt:{opt elas:ticities(type)}}{cmd:households} (the default), {cmd:individuals} or {cmd:market}{p_end}
{synopt:{opt hh:size(varname)}}size of the household; with it, the default is {cmd:individuals}{p_end}

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
{synopt:{opt st:ars}}significance stars on the elasticities{p_end}
{synopt:{opt saveres(filename)}}write the tables to a file: {cmd:.docx}, {cmd:.tex}, {cmd:.xlsx}, {cmd:.csv} or {cmd:.md}{p_end}
{synopt:{opt notab:le}}do not display the tables{p_end}
{synopt:{opt nolog}}suppress the iteration log{p_end}

{syntab:Advanced}
{synopt:{opt compat}}reproduce the R package {bf:easi} 0.21 bit for bit,
including its bugs{p_end}
{synopt:{opt tol:erance(#)}}convergence tolerance; default {cmd:tolerance(1e-6)}{p_end}
{synopt:{opt iter:ate(#)}}maximum iterations; default {cmd:iterate(100)}{p_end}
{synoptline}
{p 4 6 2}{opt aweight}s, {opt pweight}s and {opt iweight}s are allowed; see
{help weight}.  {opt pweight}s imply {cmd:vce(robust)}.  {opt fweight}s are not
allowed (see {help easi##weights:Weights}).{p_end}

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
{cmd:e()}.  The elasticities are means over the households of the elasticities
of each household, of three types: those of the household, of the individual
or of the market, the total demand (see
{help easi##types:Which elasticities}).

{p 4 4 2}
Survey data record zero shares for the households that do not buy a good, and
no price for them.  {opt pimpute()} fills the missing prices from the other
households of the same group, and {opt selection} corrects the system for the
non-buyers by the two-step method of Shonkwiler and Yen (1999); see
{help easi##selection:The selection of the buyers}.  {helpb easidiag} diagnoses
a specification, the probits included, before estimating it.


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
listed demographics, which must be among {opt demographics()}.  The default is
all of them.  {opt interpz()} implies {opt pz}.
{bf:This option behaves differently from the R package and from}
{helpb sr_easi} -- see
{help easi##compat:Differences} below.

{phang}
The lists are checked before anything is estimated: a variable may appear only
once in a list (shares, prices, demographics, {opt interpz()}), and in one role
only -- a share, a price, the expenditure or a demographic; the names of
{opt snames()} must be distinct, one per good.

{dlgtab:Missing prices and non-buyers}

{phang}
{opt pimpute(varlist)} fills the missing prices.  Without it, a household with
a missing price leaves the sample; with unit values as prices, that drops the
households that do not buy some good, and a correction for the non-buyers is
then impossible.  Each missing log price is replaced by the weighted mean of
the log prices of the households of the sample in the same group of the first
variable (for example the PSU); when nobody in the group has a price, the
group of the next variable is used (for example urban/rural), and so on.
Households still without a price leave the sample, and the means are computed
again until the sample no longer changes: the donors are the households of the
final sample.  The weight is that of the estimation (times {opt hhsize()} under
{cmd:elasticities(individuals)}).  A note reports, good by good, the prices
filled at each level.  The analytic standard errors include the imputation:
a household that gives its price moves the mean of its group, hence the
prices it fills, and its influence function carries that change through the
system, the probits and the elasticities.  The term sums to zero within each
group, so it cancels under {cmd:vce(cluster)} or {cmd:vce(svy)} at the level
of the first grouping variable; {cmd:vce(bootstrap)} fills the prices again on
every replication.  {cmd:predict} and {cmd:estat engel} fill them in the same
way on {cmd:e(sample)}.

{phang}
{opt selection} corrects for the households that do not buy (Shonkwiler and
Yen 1999).  For every corrected good {it:i}, a probit of purchase
({it:w_i} > 0) is estimated first, with the weights of the estimation, on

{p 12 12 2}
{it:s_i} = (1, ln {it:p_1} - ln {it:x}, ..., ln {it:p_J} - ln {it:x}, {it:z}, the variables of {opt selvars()} for {it:i}),

{pmore}
and the system is then estimated on all the households, buyers or not, with
the expected shares

{p 12 12 2}
E[{it:w_i}] = Phi({it:s_i}'{it:a_i}) {it:f_i} + {it:delta_i} phi({it:s_i}'{it:a_i}),

{pmore}
{it:f_i} the share of the EASI model and {it:delta_i} a new parameter, the
covariance of the error of the share with that of the probit; {it:delta_i} = 0
is no selection.  The last good closes the system, E[{it:w_J}] = 1 - the sum
of the others: it must be bought by every household (place last a good that
everybody buys, such as the rest of the budget), and it is not corrected.  The
regressors of a corrected equation are those of EASI times Phi, and phi; its
instruments, those of EASI times Phi, and phi.  The implicit utility {it:y} is
computed as without selection, from the latent shares: those of the non-buyers
are replaced by their expectation {it:f_i} - {it:delta_i} phi/(1 - Phi), the
last good closing the budget (a fixed point in {it:y}).  Symmetry is imposed on
the latent shares {it:f}, not on the expected shares.  The probit is written in
prices relative to expenditure, so that it is homogeneous of degree zero, as
are the expected shares.  The elasticities are those of the expected demand of
all the households, buyers and non-buyers, with the expected shares as
denominators: Engel and Cournot aggregation hold exactly for every type.  The
standard errors include the estimation of the probits (their influence
functions are stacked into those of the system).  {opt selection} is refused
with {cmd:vce(conventional)} and with {opt compat}.

{phang}
{opt selgoods(namelist)} restricts the correction to the goods listed (names
of the shares).  By default, every good with zero shares in the sample is
corrected, the last good excepted; a good listed that every household buys is
left uncorrected, with a note.

{phang}
{cmd:selvars(}[{it:good}{cmd::}] {it:varlist} [{cmd:;} ...]{cmd:)} gives the
variables of the probits only: they move the purchase of the good but not the
share of those who buy (exclusion restrictions).  Segments are separated by
{cmd:;}.  A segment without a good name goes to the probit of every corrected
good; {it:good}{cmd::} {it:varlist} to the probit of that good only; the two
add up.  With the goods {cmd:w1} and {cmd:w2} corrected:

{p2colset 12 44 46 2}{...}
{p2col:{cmd:selvars(age)}}{cmd:age} in the two probits{p_end}
{p2col:{cmd:selvars(w1: isMale)}}{cmd:isMale} in the probit of {cmd:w1} only{p_end}
{p2col:{cmd:selvars(age ; w1: isMale)}}{cmd:age} in the two, {cmd:isMale} for {cmd:w1}{p_end}
{p2col:{cmd:selvars(w1: age ; w2: isMale)}}{cmd:age} for {cmd:w1}, {cmd:isMale} for {cmd:w2}{p_end}
{p2colreset}{...}

{pmore}
Rules, checked before estimating: a good named must be one of the shares, not
the last one, and corrected (in {opt selgoods()} when it is given); one good
per segment; a variable of {opt selvars()} must not be in the model (a share,
price, expenditure or demographic).  Households with a missing value of these
variables leave the sample, with a note.  Without {opt selvars()}, the
correction is identified by the nonlinearity of the probit only.  The names
and the rules are those of {cmd:equaids} and {cmd:duvm}.

{marker vce}{dlgtab:SE/Robust}

{phang}
{opt vce(vcetype)} selects the variance estimator.  {cmd:robust}, a sandwich
estimator, is the default; {cmd:cluster} {it:clustvar} allows for correlation
within clusters; {cmd:svy} reads the design declared by {helpb svyset};
{cmd:conventional} is the homoskedastic 3SLS covariance matrix that
{helpb reg3}, {cmd:systemfit} and the R package report, and is what
{opt compat} uses.  {cmd:robust} and {cmd:cluster} carry the small-sample
factors {it:n}/({it:n}-1) and {it:n_c}/({it:n_c}-1), so that {cmd:robust} is the
{cmd:cluster} variance with one household per cluster and the {cmd:svy} variance
of a simple random sample.

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

{marker bootopts}{...}
{phang}
{cmd:vce(bootstrap} [{cmd:,} {it:boot_opts}]{cmd:)} estimates the full sample,
then resamples the households with replacement, and on each replication runs
the whole procedure again: the prices filled by {opt pimpute()}, the probits of
{opt selection}, the iteration on {it:y}, the system.  The households keep
their weights.  The point estimates are those of the full sample; the variance
of the coefficients (and of delta), the standard errors of the elasticities of
every type computed, and {cmd:e(V_sel)} with {opt selection}, are those of the
replications.  A replication that does not converge is dropped and counted
({cmd:e(N_reps_ok)}).  {it:boot_opts} are

{p2colset 12 30 32 2}{...}
{p2col:{opt r:eps(#)}}replications; default {cmd:reps(200)}{p_end}
{p2col:{opt seed(#)}}the random-number seed, for reproducible draws{p_end}
{p2col:{opt str:ata(varname)}}resample within strata{p_end}
{p2col:{opt psu(varname)}}resample the primary sampling units, whole{p_end}
{p2col:{opt svy}}take the strata, the PSUs and the weight from {helpb svyset}{p_end}
{p2colreset}{...}

{pmore}
It costs one estimation per replication, and its standard errors carry the
noise of a finite number of draws: about 1/sqrt(2B) of their value, 7% at 100
replications.  The analytic standard errors are exact derivatives of the
procedure (see {help easi##selection:The selection of the buyers}); the
bootstrap is the reference when identification is weak.

{dlgtab:Elasticities}

{phang}
{opt elasticities(type)} sets which elasticities are reported.  All three are
means over the households of the elasticities of each household, at its own
prices, expenditure and demographics; they differ by the weight of each
household in the mean (see {help easi##types:Which elasticities}).

{phang2}
{cmd:households}, the default, gives those of the household: each household
counts for its weight.

{phang2}
{cmd:individuals} gives those of the individual: each household counts for its
weight times its size, {opt hhsize()}, in the whole estimation (the
coefficients, their variance and the means).  It is an estimation, not a
display: after it, only {cmd:individuals} can be displayed, and it cannot be
displayed after another estimation.  The weight used is stored in
{cmd:e(wexp)}, so that {cmd:estat engel} uses the same.

{phang2}
{cmd:market} gives those of the total demand: each household counts for its
weight times its total expenditure, which makes each elasticity a ratio of
totals, the response of the total demand of the population.  It comes from the
same estimation as {cmd:households} and both are stored, so that
{cmd:easi, elasticities(market)} displays the one after the other.

{pmore}
Any abbreviation is accepted, and {cmd:mkt} for {cmd:market}.  {opt compat}
reproduces the R package, which has the elasticities of the households only:
the other two are refused with it.

{phang}
{opt hhsize(varname)} gives the size of the household, positive.  With it and
without {opt elasticities()}, the elasticities are those of the individuals.
It does not enter the model (add it to {opt demographics()} for that); with
{cmd:elasticities(households)} or {cmd:elasticities(market)} it is not used.

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
{opt stars} marks the elasticities with significance stars (* 10%, ** 5%,
*** 1%), from the ratio of each elasticity to its standard error: normal
under {cmd:vce(robust)} and {cmd:vce(cluster)}, t on the design degrees of
freedom (PSUs minus strata) under {cmd:vce(svy)}.  The standard errors keep
their own table (02-b, 03-b, ...).

{phang}
{opt saveres(filename)} writes all the elasticity tables of the display to one
file, whose extension gives the format: {cmd:.docx} (Word), {cmd:.tex} (LaTeX,
booktabs), {cmd:.xlsx} (Excel, the values as numbers), {cmd:.csv} or
{cmd:.md} (Markdown).  An existing file is replaced.

{phang}
{opt notable} displays the header only, not the tables; with {opt saveres()},
the tables are still written to the file.

{phang}
All of these can also be given when replaying results, as in
{cmd:. easi, compensated} or {cmd:. easi, saveres(elasticities.docx)}, and so
can {opt elasticities()}, within the limits above.

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
The elasticity matrices are displayed in one piece, whatever the number of
goods and whatever {help linesize}, with the good names intact; the Results
window scrolls a wide table sideways.  The display is the same on every
version of Stata from 14.2.

{p 4 4 2}
The legacy matrices do not.  {cmd:e(elast_price)}, kept for backward
compatibility, is indexed {bf:[price, good]} -- the transpose -- and
{cmd:e(compensated_q)} is indexed [good, price] but holds the compensated
elasticities {bf:plus the identity matrix}, so its diagonal is every own-price
elasticity plus one.  Both come from the R package.  Because both axes carry
the same good names, nothing in a printed table signals which way to read it,
which is why the current names exist and the display states the orientation.

{marker types}{...}
{dlgtab:Which elasticities}

{p 4 4 2}
Each household has its own elasticities, at its own prices, expenditure and
demographics.  {cmd:easi} reports their mean over the households, each
weighted by the share of the good in the budget of the household and by the
weight of the household: {it:omega_h} for {cmd:households},
{it:omega_h n_h} for {cmd:individuals} ({it:n_h} = {opt hhsize()}),
{it:omega_h x_h} for {cmd:market} ({it:x_h} = total expenditure).  For the
expenditure elasticity, with {it:a_hj} = d{it:w_hj}/d ln {it:x},

{p 8 8 2}
sum_h omega_h w_hj eta_hj / sum_h omega_h w_hj = 1 + sum_h omega_h a_hj / sum_h omega_h w_hj,

{p 4 4 2}
a ratio of two means, which never divides by the share of a single household
(the plain mean of the household elasticities does, and explodes when some
shares are near zero).  With the weights {it:omega_h x_h}, numerator and
denominator are totals of spending: the ratio is exactly the elasticity of the
total demand for the good, whatever the model.  The price elasticities are
formed the same way; the demographic semi-elasticities are plain weighted
means.  Engel and Cournot aggregation hold exactly for every type, with the
mean shares of that type.

{p 4 4 2}
They are not the elasticities of a {it:reference household}, at the means of
prices, expenditure and demographics.  Such a household depends on which mean
is taken: the mean of ln {it:x} and the log of the mean of {it:x} differ by about
half the variance of ln {it:x}, 17% of expenditure on the Canadian data and 45%
on the Mexican data.  The model is nonlinear in {it:y}, so the elasticity at the
mean is not the mean of the elasticities; and {it:y} being implicit, no
household has the mean one.

{p 4 4 2}
The compensated elasticity of household {it:h} is
(d{it:w_hj}/d ln {it:p_k} at constant {it:y} + {it:w_hj w_hk})/{it:w_hj} minus
one on the diagonal.  Its mean so weighted carries the mean of the products of
the shares, not the product of the mean shares that the R
package use; the two differ by cov({it:w_j},{it:w_k})/{it:wbar_j}, from 0.03 to
0.07 on the own-price elasticities of the Canadian data.  {cmd:e(slutsky)} is
likewise the mean of the Slutsky matrices of the households, negative
semidefinite whenever each of them is.  {opt compat} keeps the product of the
means.

{dlgtab:Standard errors of the elasticities}

{p 4 4 2}
Every reported elasticity is a ratio of two means over the households (a mean,
for the demographic semi-elasticities), and both are estimated twice over: the
coefficients inside them are estimated, and the households they are taken over
are a sample.  The standard error comes from the influence function, which
has the two terms: the sampling of the households, and the estimation of the
coefficients times the gradient of the elasticity.  The gradient of the whole
elasticity vector is taken at once, in closed form, so the dropped good is
treated like any other and the covariances between equations enter as they
should; it includes the dependence through the implicit utility {it:y}, which
is itself a function of the coefficients.  The influence function is aggregated by the estimator of
{opt vce()}: robust, clustered, or the survey design under {cmd:vce(svy)}.  Under
{cmd:vce(conventional)}, the coefficient term is the conventional one.

{p 4 4 2}
The coefficient term alone misstates the standard error.  Numerator and
denominator move together: for the expenditure elasticities the cross term is
negative and the omission {it:overstates} the standard error when the budget
shares are small.  Re-running the whole estimation
for each household with its weight perturbed -- the influence function by
brute force, in which every estimate nested in the procedure moves (the mean
shares of the instrument, the first pass that rebuilds it, the first stage,
Sigma, the iteration on {it:y}) -- gives the reported standard errors within
0.16%.  On a case built so that the model holds (700 samples of 3,000
households), the mean standard error is between 0.96 and 1.03 times the
standard deviation of the estimates, for every elasticity of the households
and of the market, and the Wald test at 5% rejects 5.4% of the time.  The
technical note gives the tables.

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
example that was two thirds of the run time.  The Jacobian is now analytic;
the standard errors -- the full influence functions of every elasticity, for
the households and for the market -- take 2.9 of the 12.4 seconds of the
reference example (9.5 without them).  When they are computed, every
elasticity table is followed by a table of its standard errors; with
{opt noelastse} those tables are simply absent.

{p 4 4 2}
{cmd:easi} reports its own execution time and stores it in {cmd:e(time)}.  It
measures this with Stata timer number 100.  Stata's timers are a single global
pool and a running timer cannot be read back, so there is no way to save one and
put it back: if you are using timer 100 yourself, {cmd:easi} will overwrite it.

{marker selection}{...}
{dlgtab:The selection of the buyers}

{p 4 4 2}
With {opt selection}, a table under the header reports, for each corrected
good, the percentage of buyers, the pseudo-R2 of its probit (McFadden), the
households predicted with probability 0 or 1 ({it:Perfect}: separation), the
variance inflation {it:VIF(delta)} and delta with its z.  {it:VIF(delta)} is
[({it:X}'{it:WX})^-1]_dd times phi'{it:W}phi, {it:X} the regressors of the
equation (those of EASI times Phi, and phi): the factor by which the
collinearity of phi with Phi times the other regressors inflates the variance of
delta, 1/(1 - {it:R}^2) with {it:R}^2 the uncentred R-squared of phi on them.
Above 10, delta rests on little more than the curvature of the probit; the
good can be left uncorrected with {opt selgoods()}, or a variable of the probit
only added in {opt selvars()}.  {helpb easidiag} reports the same table before
estimating, at the starting point (the Stone index for {it:y}).  Under the
table, the size of the completion of {it:y}: the percentage of the households
with a zero share of a corrected good, and the mean absolute change of their
{it:y} against the {it:y} computed with the observed zeros, in units and in
standard deviations of {it:y}.  The probits
are estimated good by good (a probit for each good, not a multivariate probit):
the two-step estimator is consistent, not efficient.

{p 4 4 2}
{bf:Standard errors.}  The influence function of the coefficients stacks the
system and the probits: the moments of the system, plus their derivative with
respect to the probit coefficients times the influence functions of the
probits, Phi, phi and {it:y} moving with them.  The influence function of a
probit uses the observed Hessian, the derivative of the score that the
estimator solves; the expected information that the scoring iterations use is
equal to it only in expectation, and on survey data the two differed by up to
4.6% of the standard errors of the system.  The elasticities carry the same
terms, plus the sampling of the households.  Re-running the whole procedure
for each household with its weight perturbed (probits, iteration on {it:y},
system, elasticities) gives the reported standard errors within 0.5% on the
Mexican data, without a variable of the probit only.  On a censored case built
so that the model holds (10,000 households, an exclusion restriction per good;
{cmd:replication/dgp_easi.do}), the mean analytic standard error is 0.93 to
1.05 times the standard deviation of the estimates over 700 samples (median
0.99), and the Wald test at 5% rejects 5.5% of the time; against a
500-replication bootstrap of the whole procedure on one sample, the median
ratio is 1.006.  delta carries a finite-sample bias of a third of its
standard deviation at that size.  The elasticities of the expected demand
agree with finite differences of E[w] to 2e-11.  The technical note gives the
tables.

{marker weights}{...}
{dlgtab:Weights}

{p 4 4 2}
{opt aweight}s and {opt pweight}s are normalised to sum to the number of
observations; {opt iweight}s are used as they stand.  All the population means
used to form the elasticities are weighted, as is the median above.  With
{cmd:elasticities(individuals)} the weight, or the {helpb svyset} weight under
{cmd:vce(svy)}, is multiplied by {opt hhsize()}, in the whole estimation (an
{opt aweight} of {opt hhsize()} when there is none).  The R package has no
weights at all.

{p 4 4 2}
{opt fweight}s are not allowed.  On household surveys they are almost always
sampling weights in disguise -- such as {cmd:fw = int(pw*10000)}, a frequency
weight built from a sampling weight for commands that accept no other -- and a
frequency weight counts each household that many times, which would make the
standard errors far too small.  A sampling weight goes in {cmd:[pweight=]} (or
{cmd:[aweight=]}), whatever its scale.


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
interacts prices with every variable in {opt interpz()}.
{bf:This changes results substantially}: in the reference example the system grows from 320 to
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
Six further changes are improvements rather than fixes, and are {bf:not}
reverted by {opt compat} except where stated.

{phang}
{bf:a. Iterated residual covariance.}  {cmd:systemfit} and {cmd:reg3} form Sigma
once and take one GLS step.  For a singular demand system the estimator is
invariant to the deleted equation only if Sigma is iterated to convergence
(Barten 1969); with one step the price elasticities move by up to 0.20 depending
on which good is dropped.  {cmd:easi} iterates.  {cmd:legacy(sigma1)}, implied by
{opt compat}, takes the single step.

{phang}
{bf:b. Robust standard errors by default}, with the factor {it:n}/({it:n}-1), and {opt vce(cluster)}.  See
{help easi##vce:SE/Robust} above.  {opt compat} uses {cmd:vce(conventional)}.

{phang}
{bf:c. Generated regressor.}  The Jacobian of the estimating equations includes
the dependence of {it:y} on the coefficients.  {opt compat} reverts it.

{phang}
{bf:d. Weights.}  {cmd:easi} supports {opt aweight}s, {opt pweight}s and
{opt iweight}s.  The R package has none.

{phang}
{bf:e. Compensated elasticities.}  The mean over the households of the
compensated elasticities, which carries the mean of the products of the shares
(see {help easi##types:Which elasticities}).  {opt compat} reverts it.

{phang}
{bf:f. Types of elasticities and standard errors.}  The elasticities of the
individual and of the market, and standard errors from the full influence
function for every elasticity, the dependence of {it:y} on the coefficients
included.  Not available with {opt compat}.

{phang}
{bf:g. Missing prices, non-buyers, bootstrap.}  {opt pimpute()},
{opt selection} and {cmd:vce(bootstrap)}.  The R package has none of them;
{opt selection} is not available with {opt compat}.


{marker examples}{title:Examples}

{pstd}
The examples use the two data sets that come with the package as ancillary
files: the Canadian data of Lewbel and Pendakur (2009) ({cmd:hixdata.dta}: nine
goods, prices and expenditure in logarithms) and the Mexican cereals
({cmd:mex_bench.dta}: three goods, a stratified two-stage design already
{helpb svyset}).  {stata "ssc install easi, all replace"} (or {cmd:net get easi})
copies them into the current folder; the links read them from there, else from
the SSC archive, else from GitHub, and write nothing to disk.  Each one
runs from its blue links: in the command window, as a do-file opened in the
Do-file Editor, and, except the fifth (commands after an estimation), in the
dialog box (filled in, prices and expenditure in levels; click OK).  The data in memory are not lost: the command window and the do-file give
them back at the end, even after an error; the dialog box, which needs the
example data in memory, refuses to replace data that have unsaved changes.
Files written by the examples go to Stata's temporary folder.  The links call
{cmd:easi_examples} {it:#} [{cmd:, db} | {cmd:do}].

{title:Example 1: The Canadian data of Lewbel and Pendakur (2009)}

{phang2}{cmd:. use hixdata, clear}{p_end}
{phang2}{cmd:. easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,}
{cmd:lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers) lnexpenditure(log_y)}
{cmd:demographics(age hsex carown) power(3)}{p_end}
{phang2}{cmd:. easi, compensated checks stars}{p_end}
{p 8 8 2}{txt}({stata "easi_examples 1":example 1: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 1, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 1, do":open as a do-file}){p_end}

{title:Example 2: The specification of Lewbel and Pendakur (2009)}

{pstd}Power 5 and all the interactions; about ten seconds.{p_end}
{phang2}{cmd:. use hixdata, clear}{p_end}
{phang2}{cmd:. easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,}
{cmd:lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers) lnexpenditure(log_y)}
{cmd:demographics(age hsex carown time tran) power(5) py zy pz}{p_end}
{p 8 8 2}{txt}({stata "easi_examples 2":example 2: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 2, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 2, do":open as a do-file}){p_end}

{title:Example 3: Survey design (Mexican cereals)}

{pstd}The design is read from {helpb svyset}: PSUs within strata, the sampling weight.{p_end}
{phang2}{cmd:. use mex_bench, clear}{p_end}
{phang2}{cmd:. svyset}{p_end}
{phang2}{cmd:. easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) py vce(svy) stars}{p_end}
{p 8 8 2}{txt}({stata "easi_examples 3":example 3: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 3, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 3, do":open as a do-file}){p_end}

{title:Example 4: Reproduce the R package easi 0.21}

{pstd}{opt compat} reproduces the R package, and results published with {helpb sr_easi}, bit for bit.{p_end}
{phang2}{cmd:. use hixdata, clear}{p_end}
{phang2}{cmd:. easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,}
{cmd:lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers) lnexpenditure(log_y)}
{cmd:demographics(age hsex carown time tran) power(5) py zy pz compat}{p_end}
{p 8 8 2}{txt}({stata "easi_examples 4":example 4: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 4, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 4, do":open as a do-file}){p_end}

{title:Example 5: After estimation}

{pstd}Fitted shares, the implicit utility index and the Engel curves; run from its link, the example writes the curves to Stata's temporary folder.{p_end}
{phang2}{cmd:. use hixdata, clear}{p_end}
{phang2}{cmd:. easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,}
{cmd:lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers) lnexpenditure(log_y)}
{cmd:demographics(age hsex carown) power(3) notable}{p_end}
{phang2}{cmd:. predict double what*, shares}{p_end}
{phang2}{cmd:. predict double yhat, y}{p_end}
{phang2}{cmd:. estat engel, n(60) data(curves)}{p_end}
{p 8 8 2}{txt}({stata "easi_examples 5":example 5: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 5, do":open as a do-file}){p_end}

{title:Example 6: The tables in a file}

{pstd}The extension of {opt saveres()} gives the format: Word, Excel, LaTeX, CSV or Markdown.{p_end}
{phang2}{cmd:. use hixdata, clear}{p_end}
{phang2}{cmd:. easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,}
{cmd:lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers) lnexpenditure(log_y)}
{cmd:demographics(age hsex carown) power(3) compensated notable saveres(tables.docx)}{p_end}
{phang2}{cmd:. easi, stars saveres(tables.xlsx) notable}{p_end}
{p 8 8 2}{txt}({stata "easi_examples 6":example 6: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 6, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 6, do":open as a do-file}){p_end}

{title:Example 7: Diagnose a specification before estimating it}

{pstd}The same specification with raw and with centred logarithms: {helpb easidiag} predicts the slow convergence of the first.{p_end}
{phang2}{cmd:. use mex_bench, clear}{p_end}
{phang2}{cmd:. easidiag w1 w2 w3, lnprices(lp1_raw lp2_raw lp3) lnexpenditure(lx_raw) demographics(z1 z2) power(3)}{p_end}
{phang2}{cmd:. easidiag w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3)}{p_end}
{p 8 8 2}{txt}({stata "easi_examples 7":example 7: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 7, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 7, do":open as a do-file}){p_end}

{title:Example 8: The elasticities of the households, of the market and of the individuals}

{pstd}The same estimation gives those of the households (the default) and of the market; {opt hhsize()} gives those of the individuals.{p_end}
{phang2}{cmd:. use mex_bench, clear}{p_end}
{phang2}{cmd:. easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) py vce(svy)}{p_end}
{phang2}{cmd:. easi, elasticities(market)}{p_end}
{phang2}{cmd:. easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) py vce(svy) hhsize(hhsize)}{p_end}
{p 8 8 2}{txt}({stata "easi_examples 8":example 8: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 8, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 8, do":open as a do-file}){p_end}

{title:Example 9: The households that do not buy}

{pstd}A household in ten buys no corn, one in six no wheat.  The diagnostic runs the probits first, with the age
and the sex of the head in the probits only; the system is then estimated on all the households with the expected
shares.  On these data delta is weakly identified (VIF(delta) above 10: it rests mostly on the curvature of the
probit), and easi warns so after the estimation with {cmd:vce(svy)}; this is expected here.  A bootstrap of the
whole procedure then gives standard errors that take it into account (about two minutes).{p_end}
{phang2}{cmd:. use mex_bench, clear}{p_end}
{phang2}{cmd:. easidiag w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) selvars(age isMale)}{p_end}
{phang2}{cmd:. easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) vce(svy) selvars(age isMale)}{p_end}
{phang2}{cmd:. predict double Ew*, shares}{p_end}
{phang2}{cmd:. predict double f*, shares latent}{p_end}
{phang2}{cmd:. summarize w1 Ew1 f1 w2 Ew2 f2}{p_end}
{phang2}{cmd:. easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) selvars(age isMale) vce(bootstrap, reps(50) seed(1) svy)}{p_end}
{p 8 8 2}{txt}({stata "easi_examples 9":example 9: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 9, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "easi_examples 9, do":open as a do-file}){p_end}


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
{synopt:{cmd:e(chk_engel_mkt)}, {cmd:e(chk_cournot_mkt)}}the same for the market elasticities{p_end}
{synopt:{cmd:e(N_reps)}, {cmd:e(N_reps_ok)}}bootstrap replications, successful ones ({cmd:vce(bootstrap)}){p_end}
{synopt:{cmd:e(sel_ycomp_pct)}}with {opt selection}: percentage of the households whose {it:y} is completed (a zero share of a corrected good){p_end}
{synopt:{cmd:e(sel_ycomp_mean)}, {cmd:e(sel_ycomp_sd)}}their mean absolute change of {it:y}, and the standard deviation of {it:y}{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:easi}{p_end}
{synopt:{cmd:e(shares)}, {cmd:e(prices)}, {cmd:e(demographics)}}variable lists{p_end}
{synopt:{cmd:e(interpz)}}demographics crossed with prices{p_end}
{synopt:{cmd:e(vce)}, {cmd:e(clustvar)}}variance estimator: {cmd:robust}, {cmd:cluster}, {cmd:svy}, {cmd:conventional} or {cmd:bootstrap}{p_end}
{synopt:{cmd:e(boot_design)}, {cmd:e(boot_seed)}}resampling unit and seed ({cmd:vce(bootstrap)}){p_end}
{synopt:{cmd:e(svyunit)}, {cmd:e(svystrata)}}PSU and strata, with {cmd:vce(svy)}{p_end}
{synopt:{cmd:e(svywvar)}, {cmd:e(svyfpc)}}survey weight and fpc variables{p_end}
{synopt:{cmd:e(vcetype2)}}{cmd:generated regressor} or {cmd:conditional on y}{p_end}
{synopt:{cmd:e(wtype)}, {cmd:e(wexp)}}weight type and expression (times {opt hhsize()} under {cmd:individuals}){p_end}
{synopt:{cmd:e(elasticities)}}{cmd:households}, {cmd:individuals} or {cmd:market}{p_end}
{synopt:{cmd:e(hhsize)}}the variable of {opt hhsize()}{p_end}
{synopt:{cmd:e(pimpute)}}the grouping variables of {opt pimpute()}{p_end}
{synopt:{cmd:e(selection)}}{cmd:shonkwiler-yen} with {opt selection}{p_end}
{synopt:{cmd:e(selgoods)}}the goods corrected; {cmd:e(selvars)} the {opt selvars()} specification{p_end}
{synopt:{cmd:e(sel_z_}{it:good}{cmd:)}}the variables of the probit of {it:good} only{p_end}
{synopt:{cmd:e(sel_zall)}}all the variables of {opt selvars()}; {cmd:e(sel_flags)}, {cmd:e(sel_zmask)} which equation
is corrected and which of them enter it (for {cmd:predict} and {cmd:estat engel}){p_end}
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
{synopt:{cmd:e(slutsky)}}Slutsky matrix in share form, mean over the households{p_end}
{synopt:{cmd:e(elast_exp_mkt)}, ...}the market elasticities and their standard errors, with the names above
followed by {cmd:_mkt} ({cmd:e(elast_exp_mkt_se)}...); not after {cmd:individuals} or under {opt compat}{p_end}
{synopt:{cmd:e(Sigma)}}cross-equation residual covariance{p_end}
{synopt:{cmd:e(sel_alpha)}}probit coefficients by good: constant, ln {it:p}, demographics, {opt selvars()}
(missing where a variable does not enter or the good is not corrected){p_end}
{synopt:{cmd:e(sel_diag)}}by good: percentage of buyers, pseudo-R2, perfectly predicted, VIF(delta){p_end}
{synopt:{cmd:e(V_sel)}}variance of the coefficients and of the probits together ({cmd:estat engel}){p_end}

{p2col 5 24 28 2: Matrices, legacy}{p_end}
{p 4 6 2}Kept unchanged so that existing do-files keep working.  See
{help easi##remarks:Reading the price elasticity tables}.{p_end}
{synopt:{cmd:e(elast_income)}, {cmd:e(elast_income_se)}}same as
{cmd:e(elast_exp)}{p_end}
{synopt:{cmd:e(elast_price)}, {cmd:e(elast_price_se)}}uncompensated price
elasticities, [price, good] -- the {bf:transpose} of {cmd:e(elast_price_nc)}{p_end}
{synopt:{cmd:e(semi_income)}}same as {cmd:e(semi_exp)}{p_end}
{synopt:{cmd:e(compensated_q)}}{cmd:e(elast_price_c)} {bf:plus the identity matrix}{p_end}


{marker predict}{title:Postestimation: predict}

{p 8 15 2}
{cmd:predict} [{it:type}] {it:stub}{cmd:*} {ifin} [{cmd:,} {opt sh:ares} [{opt lat:ent}] | {opt res:iduals}]{p_end}
{p 8 15 2}
{cmd:predict} [{it:type}] {it:newvar} {ifin}{cmd:,} {opt y}{p_end}

{phang}{opt shares} (the default) stores one fitted budget share per good.  The
last one is obtained by adding up, so the fitted shares sum to 1 exactly.{p_end}

{phang}{opt residuals} stores observed minus fitted shares, one per good.{p_end}

{phang}{opt y} stores the implicit utility index.{p_end}

{phang}After {opt selection}, {opt shares} are the expected shares
E[{it:w_i}] = Phi {it:f_i} + {it:delta_i} phi (they sum to 1 exactly),
{opt residuals} are {it:w} - E[{it:w}], {opt y} is the implicit utility of the
completed latent shares, and {opt latent} with {opt shares} stores the latent
shares {it:f} of the model instead.{p_end}

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
{bf:After selection}, the curves are those of the expected shares,
Phi {it:f} + {it:delta} phi, the probits at the means of prices, demographics
and variables of {opt selvars()} too; the band includes the estimation of the
probits ({cmd:e(V_sel)}).

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

{p 4 8 2}Shonkwiler, J. S., and S. T. Yen. 1999.
Two-step estimation of a censored system of equations.
{it:American Journal of Agricultural Economics} 81: 972-982.{p_end}

{p 4 8 2}Zhen, C., E. A. Finkelstein, J. M. Nonnemaker, S. A. Karns, and J. E. Todd. 2014.
Predicting the effects of sugar-sweetened beverage taxes on food and beverage demand in a large demand system.
{it:American Journal of Agricultural Economics} 96: 1-25.{p_end}

{p 4 8 2}Pendakur, K. 2008. EASI made Easier.
{browse "http://www.sfu.ca/~pendakur/"}{p_end}

{p 4 8 2}Hoareau, S., G. Lacroix, M. Hoareau, and L. Tiberti. 2012.
Exact Affine Stone Index Demand System in R: The easi Package.
CIRPEE working paper.{p_end}


{title:Author}

{pstd}Abdelkrim Araar, Universite Laval / PEP{break}
{browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}{p_end}

{pstd}Version 2.0.1, 30 September 2026 (first release 1.0.0, 23 September 2026).  License: GPL-3.0-or-later.
Source, examples, technical note and replication files:
{browse "https://github.com/aabbdd12/easi"}{p_end}


{title:Also see}

{psee}{helpb sr_easi} (compatibility wrapper for the old syntax){p_end}
