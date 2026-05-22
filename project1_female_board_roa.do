********************************************************************************
* PROJECT: Female Board Representation and Firm Performance
* METHOD:  Propensity Score Matching, CEM, IPWRA (Doubly Robust)
* AUTHOR:  Mariia Dorfman
* DATA:    Cross-sectional panel of publicly listed firms, 2010-2019
*          Source: Refinitiv (TRBC classification, board composition, financials)
* OUTPUT:  Tables 1-3 (summary stats, propensity score, treatment effects)
*          Figures: PS overlap, ROA density after matching, dynamic effects
********************************************************************************

********************************************************************************
* SETUP
********************************************************************************

clear all
set more off
set varabbrev off

* ----> Set your working directory here <----
* cd "/your/path/here"

* Data file should be named "data_ps2.dta" in your working directory
use "data_ps2.dta", clear


********************************************************************************
* STEP 1: DATA CLEANING
********************************************************************************

* Drop anticipation period (firms may respond to rumored policy before adoption)
drop if year >= 2011 & year <= 2014

* Drop observations with missing outcome or treatment
drop if missing(female)
drop if missing(roa)

* Winsorize all continuous variables at 1st and 99th percentiles
* to limit influence of outliers
winsor2 roa trev cogs gp gm sga labcosts sgoa rd capex opex opin opm ni ///
    div roe ta ca cl tinv ppe intang stdebt ltdebt tdebt tliab teq        ///
    marcap emp beta esgs ind dirnetwork age nquals boardsize,              ///
    cuts(1 99) replace

summarize
tab year


********************************************************************************
* STEP 2: CONVERT TO CROSS-SECTION
********************************************************************************
* Unit of observation: one firm
* Outcome: average ROA over post-period 2015-2019
* Treatment: firm has at least one female director post-recommendation
* Covariates: firm characteristics as of 2010 (pre-treatment baseline)

* Average ROA in post-period (2015-2019)
gen post = inrange(year, 2015, 2019)
bys id: egen roa_post_mean = mean(roa) if post == 1
bys id: egen roa_post = max(roa_post_mean)
drop roa_post_mean

* Baseline ROA in 2010
gen roa_2010 = roa if year == 2010
bys id: egen roa2010 = max(roa_2010)
drop roa_2010

* Change in ROA (used as robustness outcome)
gen delta_roa = roa_post - roa2010

* Treatment indicator: firm appoints female director in 2015 or later
gen treated_temp = female if year >= 2015
bys id: egen treated = max(treated_temp)
drop treated_temp

* Keep only 2010 baseline row (one observation per firm)
keep if year == 2010

keep id treated roa_post delta_roa roa2010 trev cogs gp gm sga labcosts   ///
    sgoa rd capex opex opin opm ni div roe ta ca cl tinv ppe intang         ///
    stdebt ltdebt tdebt tliab teq marcap emp beta esgs ind dirnetwork age   ///
    nquals boardsize TRBCEconomicSectorName TRBCBusinessSectorName          ///
    TRBCIndustryGroupName

* Remove firms with missing treatment or outcome
duplicates report id
drop if missing(treated)
drop if missing(roa_post)

* Median imputation for covariates with missing values
global controls roa2010 trev gm capex ta marcap emp boardsize ind dirnetwork age nquals

foreach var of global controls {
    summarize `var', detail
    replace `var' = r(p50) if missing(`var')
}


********************************************************************************
* STEP 3: SUMMARY STATISTICS (Table 1)
********************************************************************************

global table1 roa_post delta_roa roa2010 trev gm capex ta marcap emp boardsize ind dirnetwork age nquals

estpost tabstat $table1 if treated == 1, statistics(mean sd min max) columns(statistics)
esttab . using "table1_treated.tex", cells("mean sd min max") label          ///
    mtitles("Treated") title("Table 1: Summary Statistics — Treated Firms")  ///
    replace bfmt(%10.3f) nonum

estpost tabstat $table1 if treated == 0, statistics(mean sd min max) columns(statistics)
esttab . using "table1_control.tex", cells("mean sd min max") label          ///
    mtitles("Control") title("Table 1: Summary Statistics — Control Firms")  ///
    replace bfmt(%10.3f) nonum

* T-tests for pre-treatment covariate balance (check for selection bias)
display "T-test p-values (Treated vs Control):"
foreach var of global controls {
    ttest `var', by(treated)
    local p = r(p)
    display "`var': p = " %6.3f `p'
}


********************************************************************************
* STEP 4: PROPENSITY SCORE ESTIMATION (Table 2)
********************************************************************************

* Encode industry sector dummies
encode TRBCEconomicSectorName, gen(sec_id)
encode TRBCBusinessSectorName,  gen(bus_id)
encode TRBCIndustryGroupName,   gen(ind_id)

* Logit for propensity score: treatment ~ baseline controls + industry FE
logit treated $controls i.sec_id i.bus_id i.ind_id
predict pscore, pr

eststo ps_logit
esttab ps_logit using "table2_pscore.tex",                                   ///
    cells("b(star fmt(3)) se(fmt(3)) z(fmt(3)) p(fmt(3))")                   ///
    label nonum title("Table 2: Propensity Score Logit") replace

summarize pscore if treated == 1
summarize pscore if treated == 0


********************************************************************************
* STEP 5: OVERLAP CHECK
********************************************************************************

* Kernel density of propensity scores by treatment status
kdensity pscore if treated == 1, generate(treated_x treated_d)
kdensity pscore if treated == 0, generate(control_x control_d)

twoway                                                                         ///
    (area treated_d treated_x, color(lavender%60) lcolor(lavender)            ///
        lwidth(medium) lpattern(solid))                                        ///
    (area control_d control_x, color(eltblue%50) lcolor(eltblue)              ///
        lwidth(medium) lpattern(solid)),                                       ///
    legend(order(1 "Treated" 2 "Control") pos(6) ring(0) col(1))              ///
    xlabel(0(0.1)1) ylabel(, angle(horizontal))                                ///
    title("Propensity Score Density: Treated vs Control")                      ///
    ytitle("Density") xtitle("Propensity Score")                               ///
    graphregion(color(white))
graph export "fig_overlap.png", replace width(1600)


********************************************************************************
* STEP 6: MATCHING AND COVARIATE BALANCE
********************************************************************************

* Trim propensity score support to [0.05, 0.95]
drop if pscore < 0.05
drop if pscore > 0.95

* Nearest-neighbor PSM (1:1, no replacement)
* ssc install psmatch2  // uncomment if not installed
psmatch2 treated, pscore(pscore) neighbor(1) common

* Rename controls for cleaner pstest labels
rename (roa2010 trev gm capex ta marcap emp boardsize ind dirnetwork age nquals) ///
       (roa2010_ revenue grossmargin capex_ totalassets marketcap employees      ///
        boardsize_n indep_directors networksize directors_age directors_qual)

global controls_clean roa2010_ revenue grossmargin capex_ totalassets           ///
    marketcap employees boardsize_n indep_directors networksize                  ///
    directors_age directors_qual

* Balance test: standardized bias before and after matching
pstest $controls_clean, both graph


********************************************************************************
* STEP 7: TREATMENT EFFECT ESTIMATION (Table 3)
********************************************************************************

* (1) Naive OLS — no controls, for comparison only
reg roa_post treated
eststo naive_ols

* (2) OLS with full controls and industry FE
reg roa_post treated $controls i.sec_id i.bus_id i.ind_id
eststo ols_controls

* (3) IPWRA — doubly robust (consistent if either outcome or treatment model correct)
logit treated $controls i.sec_id i.bus_id i.ind_id
predict pscore_ipwra, pr
keep if pscore_ipwra >= 0.05 & pscore_ipwra <= 0.95

teffects ipwra (roa_post $controls) (treated $controls i.sec_id i.bus_id i.ind_id), ///
    vce(robust)
eststo ipwra_teffects

* (4) Nearest-neighbor matching with caliper (1:1, caliper = 0.05 SD of pscore)
psmatch2 treated, outcome(roa_post) pscore(pscore) neighbor(1) caliper(0.05) common ties

gen weight_nn  = _weight
gen matched_nn = (_support == 1)

reg roa_post treated [pw=weight_nn] if matched_nn == 1, robust cluster(id)
eststo nn_matching

esttab nn_matching, b(3) se(3) star(* 0.10 ** 0.05 *** 0.01)                 ///
    keep(treated) scalars(N r2)                                                ///
    title("Nearest-Neighbor Matching (1:1, caliper = 0.05)")

* (5) Coarsened Exact Matching (CEM)
cem roa2010_ grossmargin totalassets marketcap boardsize_n directors_age,     ///
    treatment(treated)

reg roa_post treated [iweight = cem_weights], robust cluster(id)
eststo cem_weighted


********************************************************************************
* STEP 8: ROBUSTNESS CHECKS
********************************************************************************

* Panel A: alternative outcome — change in ROA (ΔROA = ROA_post - ROA_2010)
teffects ipwra (delta_roa $controls) (treated $controls i.sec_id i.bus_id i.ind_id), ///
    vce(robust)
eststo panelA

* Panel B: placebo test — "effect" on pre-treatment ROA (should be zero)
global controls_no_roa2010 trev grossmargin capex_ totalassets marketcap        ///
    employees boardsize_n indep_directors networksize directors_age directors_qual

teffects ipwra (roa2010_ $controls_no_roa2010)                                 ///
    (treated $controls_no_roa2010 i.sec_id i.bus_id i.ind_id), vce(robust)
eststo panelB

* Panels C/D: heterogeneity by firm size (split at median market cap)
summarize marketcap, detail
local med = r(p50)
gen large = (marketcap > `med') if !missing(marketcap)

reg roa_post treated $controls i.sec_id i.bus_id if large == 1, vce(robust)
eststo panelC_large

reg roa_post treated $controls i.sec_id i.bus_id if large == 0, vce(robust)
eststo panelD_small

esttab panelA panelB panelC_large panelD_small using "table_robustness.tex",   ///
    keep(treated) b(3) se(3) star(* 0.10 ** 0.05 *** 0.01)                    ///
    mtitles("ΔROA" "Placebo: ROA2010" "Large firms" "Small firms")             ///
    title("Robustness Checks: Female Board → ROA") label replace


********************************************************************************
* STEP 9: FIGURES
********************************************************************************

preserve
set more off

use "data_ps2.dta", clear

* Reconstruct variables for dynamic effects plot
drop if year >= 2011 & year <= 2014
drop if missing(roa) | missing(female)

gen treated_temp = female if year >= 2015
bys id: egen treated = max(treated_temp)
drop treated_temp

gen roa2010_temp = roa if year == 2010
bys id: egen roa2010 = max(roa2010_temp)
drop roa2010_temp

foreach y in 2015 2016 2017 2018 2019 {
    gen roa_temp = roa if year == `y'
    bys id: egen roa`y' = max(roa_temp)
    drop roa_temp
}

gen post = inrange(year, 2015, 2019)
bys id: egen roa_post_mean = mean(roa) if post == 1
bys id: egen roa_post = max(roa_post_mean)
drop roa_post_mean

keep if year == 2010

global controls roa2010 trev gm capex ta marcap emp boardsize ind dirnetwork age nquals

* IPW weights for dynamic effects
logit treated $controls
predict pscore, pr
gen ipw = treated/pscore + (1 - treated)/(1 - pscore)

eststo clear

* Year-by-year IPW-weighted regressions (event-study style)
foreach y in 2015 2016 2017 2018 2019 {
    gen treated`y' = treated
    reg roa`y' treated`y' [pw=ipw], robust
    eststo model`y'
}

* Dynamic treatment effects plot
coefplot model2015 model2016 model2017 model2018 model2019,                    ///
    keep(treated2015 treated2016 treated2017 treated2018 treated2019)          ///
    rename(treated2015="2015" treated2016="2016" treated2017="2017"            ///
           treated2018="2018" treated2019="2019")                              ///
    vertical                                                                   ///
    yline(0, lcolor(black) lpattern(dash) lwidth(medium))                     ///
    ciopts(recast(rcap) lcolor(lavender) lwidth(medium))                      ///
    mcolor(lavender) msymbol(circle) msize(medium)                             ///
    title("Dynamic Treatment Effects on ROA", size(medium) color(black))      ///
    ytitle("Treatment Effect", size(small))                                    ///
    graphregion(color(white)) legend(off)
graph export "fig_dynamic_effects.png", replace width(1600)

* ROA outcome distribution after matching
psmatch2 treated $controls, outcome(roa_post) neighbor(1) common
gen matched = (_weight > 0)

twoway                                                                         ///
    (kdensity roa_post if treated==1 & matched==1 & inrange(roa_post,-100,100), ///
        lwidth(medthick) lcolor(lavender))                                     ///
    (kdensity roa_post if treated==0 & matched==1 & inrange(roa_post,-100,100), ///
        lpattern(dash) lwidth(medthick) lcolor(eltblue)),                      ///
    xscale(range(-100 100)) xlabel(-100(25)100)                                ///
    legend(order(1 "Treated" 2 "Control") pos(5) ring(0))                     ///
    title("ROA Distribution After Matching")                                   ///
    xtitle("Average ROA (2015-2019)") ytitle("Density")                        ///
    graphregion(color(white))
graph export "fig_roa_density.png", replace width(1600)

restore

********************************************************************************
* END OF FILE
********************************************************************************
