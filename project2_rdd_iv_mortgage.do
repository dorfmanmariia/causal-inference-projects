********************************************************************************
* PROJECT: Credit Scores, Interest Rates, and Mortgage Default
* METHOD:  Regression Discontinuity Design (Fuzzy RDD) + IV / 2SLS
* AUTHOR:  Mariia Dorfman
* DATA:    US mortgage-level data with FICO scores, loan characteristics,
*          and delinquency outcomes
* DESIGN:  FICO = 620 threshold creates discontinuity in interest rates;
*          used as instrument for 2SLS estimation of rate → default
* OUTPUT:  Tables 1-5 (descriptives, first stage, validity checks,
*          main results, heterogeneity)
*          Figures 1-4 (binned scatter, McCrary, delinquency by broker,
*          bandwidth sensitivity)
********************************************************************************

********************************************************************************
* SETUP
********************************************************************************

clear all
set more off
capture log close

* ----> Set your working directory here <----
* cd "/your/path/here"

* Data file should be named "data_ps4.dta" in your working directory
use "data_ps4.dta", clear


********************************************************************************
* STEP 0: DATA LOADING AND DIAGNOSTICS
********************************************************************************

describe
codebook, compact
count

* --- Duplicates ---
duplicates report
duplicates drop
* Note: 172 duplicate rows (0.08%) dropped — likely a data generation artifact

* --- Missing values ---
misstable summarize
misstable patterns
tab del if missing(rate)
tab del if missing(fico)

* --- Check whether missingness in covariates is discontinuous at FICO=620 ---
* Discontinuous missingness would threaten RDD validity if covariates are included
gen above = (fico >= 620)

foreach var of varlist age gender income ltv {
    gen miss_`var' = missing(`var')
    reg miss_`var' above if abs(fico - 620) <= 40, cluster(area)
    di "Discontinuity in missingness of `var' at FICO=620:"
}
* Result: age, gender, income show discontinuous missingness at threshold;
* LTV is smooth (p=0.48). Main specs run with and without these covariates.

* --- Summary statistics ---
summarize age area balance broker del fico gender hard income ltv rate refinance white, detail
tabstat age balance fico income ltv rate, by(del) stat(mean sd p25 p50 p75 n) columns(statistics) longstub

foreach var of varlist del broker gender hard refinance white {
    di _n "Frequency: `var'"
    tab `var', miss
}

* --- Logical consistency checks ---
count if fico < 300 | fico > 850   // valid FICO range
count if age < 18
count if rate <= 0 | rate >= 100
count if balance <= 0
count if ltv <= 0
count if income < 0
* All return 0 — data passes integrity checks

* --- Sample size around threshold ---
foreach bw in 10 20 30 40 {
    count if abs(fico - 620) <= `bw'
    di "N within ±`bw' FICO points: " r(N)
}

* --- Cross-tabulations ---
tab del broker,    cell row col
tab del hard,      cell row col
tab del refinance, cell row col

save "mortgage_clean.dta", replace


* --- Distribution checks and variable transformations ---
foreach var of varlist balance income ltv rate {
    histogram `var', normal title("`var'") kdensity
    graph export "hist_`var'.png", replace
    di _n "Skewness: `var'"
    sktest `var'
    quietly summarize `var', detail
    di "  Skewness coefficient = " r(skewness)
}

* Log transformations for skewed variables
gen log_balance = log(balance)
label variable log_balance "Log loan balance"

gen log_income = log(income + 1)
label variable log_income "Log(income + 1)"

* --- Correlation and multicollinearity ---
pwcorr rate fico log_balance log_income ltv age broker gender hard refinance white, star(0.05)

reg rate fico log_balance log_income ltv age broker gender hard refinance white ///
    if abs(fico - 620) <= 40, cluster(area)
vif

save "mortgage_clean.dta", replace

* Descriptive statistics table
estpost summarize del rate fico log_balance log_income ltv age broker gender hard refinance white, detail
esttab using "table1_descriptives.tex", replace                                ///
    cells("mean(fmt(3)) sd(fmt(3)) p25(fmt(3)) p50(fmt(3)) p75(fmt(3))")      ///
    label title("Table 1: Descriptive Statistics")                             ///
    collabels("Mean" "Std. Dev." "p25" "Median" "p75") nomtitle nonumber


********************************************************************************
* STEP 1: VISUAL DIAGNOSTICS
********************************************************************************

use "mortgage_clean.dta", clear

* Center running variable at threshold
gen fico_c = fico - 620
label variable fico_c "FICO score (centered at 620)"

* Bin running variable for scatter plot
gen fico_bin = floor(fico / 2) * 2

bysort fico_bin above: egen mean_rate_bin = mean(rate)

* Figure 1: Binned scatter — interest rate vs FICO score
twoway                                                                         ///
    (scatter mean_rate_bin fico_bin if above==0 & abs(fico_c)<=40,            ///
        msize(small) mcolor(eltblue) msymbol(circle))                         ///
    (scatter mean_rate_bin fico_bin if above==1 & abs(fico_c)<=40,            ///
        msize(small) mcolor(lavender) msymbol(circle))                        ///
    (lfit mean_rate_bin fico_bin if above==0 & abs(fico_c)<=40,               ///
        lcolor(eltblue) lwidth(medthick))                                      ///
    (lfit mean_rate_bin fico_bin if above==1 & abs(fico_c)<=40,               ///
        lcolor(lavender) lwidth(medthick)),                                    ///
    xline(620, lcolor(black) lpattern(dash))                                   ///
    legend(order(1 "Below 620" 2 "Above 620") position(1) ring(0))            ///
    xtitle("FICO Score") ytitle("Mean Interest Rate")
graph export "fig1_binscatter_rate_fico.png", replace width(1600)

* FICO distribution histogram (visual check for bunching)
histogram fico, discrete frequency                                             ///
    xline(620, lcolor(red) lpattern(dash) lwidth(medthick))                   ///
    xtitle("FICO Score") ytitle("Frequency")                                   ///
    title("Distribution of FICO Scores")                                       ///
    color(navy%70)
graph export "fig_fico_histogram.png", replace width(1600)

save "mortgage_clean.dta", replace


* --- First stage: effect of FICO≥620 on interest rate ---
* Polynomial orders p=1,2,3 × bandwidths

estimates clear

foreach bw in 20 30 40 60 80 {
    eststo linear_`bw': reg rate above c.fico_c above#c.fico_c               ///
        if abs(fico_c) <= `bw', cluster(area)
}

foreach bw in 20 30 40 60 80 {
    eststo quad_`bw': reg rate above c.fico_c above#c.fico_c                 ///
        c.fico_c#c.fico_c above#c.fico_c#c.fico_c                            ///
        if abs(fico_c) <= `bw', cluster(area)
}

foreach bw in 40 60 80 {
    eststo cubic_`bw': reg rate above c.fico_c above#c.fico_c                ///
        c.fico_c#c.fico_c above#c.fico_c#c.fico_c                            ///
        c.fico_c#c.fico_c#c.fico_c above#c.fico_c#c.fico_c#c.fico_c         ///
        if abs(fico_c) <= `bw', cluster(area)
}

esttab linear_20 linear_30 linear_40 linear_60 linear_80                      ///
       quad_40 quad_60 quad_80 cubic_40 cubic_60 cubic_80,                    ///
    keep(above) b(4) se(4) star(* 0.1 ** 0.05 *** 0.01)                       ///
    mtitles("±20" "±30" "±40" "±60" "±80" "±40" "±60" "±80" "±40" "±60" "±80") ///
    mgroups("Linear" "Quadratic" "Cubic", pattern(1 0 0 0 0 1 0 0 1 0 0))    ///
    scalars("N Observations" "r2 R-squared")                                   ///
    title("Table 2: First Stage — FICO≥620 → Interest Rate")                  ///
    note("Clustered SE by area.")

* MSE-optimal bandwidth selection (rdrobust)
rdrobust rate fico, c(620) kernel(uniform)    bwselect(mserd)
rdrobust rate fico, c(620) kernel(triangular) bwselect(mserd)
rdrobust rate fico, c(620) kernel(triangular) bwselect(cerrd)


********************************************************************************
* STEP 2: RDD VALIDITY CHECKS
********************************************************************************

* --- McCrary density test (no bunching at threshold) ---
* ssc install rddensity   // uncomment if not installed
* ssc install lpdensity   // uncomment if not installed
rddensity fico, c(620)

preserve
keep if abs(fico_c) <= 23   // bandwidth chosen by rddensity
collapse (count) n=del, by(fico)
quietly sum n
gen density = n / r(sum)

twoway                                                                         ///
    (bar density fico if fico < 620,  color(eltblue%60) barwidth(0.8))        ///
    (bar density fico if fico >= 620, color(lavender%80) barwidth(0.8)),      ///
    xline(620, lcolor(black) lpattern(dash) lwidth(medthick))                 ///
    legend(order(1 "Below 620" 2 "Above 620") position(1) ring(0) cols(1))   ///
    xtitle("FICO Score") ytitle("Density")
graph export "fig2_mccrary.png", replace width(1600)
restore

* --- Covariate balance at threshold (Table 3) ---
* Covariates should not jump at FICO=620 under local randomization assumption

estimates clear
local covariates "log_balance log_income ltv age broker gender hard refinance white"

foreach var of local covariates {
    eststo bal_`var': reg `var' above c.fico_c above#c.fico_c                ///
        if abs(fico_c) <= 40, cluster(area)
}

esttab bal_log_balance bal_log_income bal_ltv bal_age bal_broker               ///
       bal_gender bal_hard bal_refinance bal_white,                            ///
    keep(above) b(4) se(4) star(* 0.1 ** 0.05 *** 0.01)                       ///
    mtitles("log(balance)" "log(income)" "LTV" "Age" "Broker"                 ///
            "Gender" "Hard inq." "Refinance" "White")                         ///
    scalars("N Observations" "r2 R-squared")                                   ///
    title("Table 3: Covariate Balance at FICO=620")                           ///
    note("BW=±40, linear polynomial, clustered SE by area.")

* --- Placebo cutoff tests ---
* No jump in rate should exist at fake thresholds

estimates clear
foreach cutoff in 590 600 610 630 640 650 {
    gen above_pl   = (fico >= `cutoff')
    gen fico_c_pl  = fico - `cutoff'
    eststo placebo_`cutoff': reg rate above_pl c.fico_c_pl above_pl#c.fico_c_pl ///
        if abs(fico_c_pl) <= 5, cluster(area)
    drop above_pl fico_c_pl
}

esttab placebo_590 placebo_600 placebo_610 placebo_630 placebo_640 placebo_650, ///
    keep(above_pl) b(4) se(4) star(* 0.1 ** 0.05 *** 0.01)                    ///
    mtitles("c=590" "c=600" "c=610" "c=630" "c=640" "c=650")                  ///
    title("Placebo Cutoff Tests — BW=±5")

* --- Donut RDD: exclude observations exactly at FICO=620 ---
eststo donut: reg rate above c.fico_c above#c.fico_c                          ///
    if abs(fico_c) <= 40 & fico != 620, cluster(area)
eststo full:  reg rate above c.fico_c above#c.fico_c                          ///
    if abs(fico_c) <= 40, cluster(area)

esttab full donut, keep(above) b(4) se(4) star(* 0.1 ** 0.05 *** 0.01)       ///
    mtitles("Full sample" "Donut (excl. FICO=620)")                            ///
    title("Donut RDD Robustness Check")


********************************************************************************
* STEP 3: MAIN RESULTS — REDUCED FORM AND 2SLS (Table 4)
********************************************************************************

estimates clear

* Reduced form: FICO≥620 → delinquency (intent-to-treat)
eststo rf_1: reg del above c.fico_c above#c.fico_c                            ///
    if abs(fico_c) <= 40, cluster(area)

eststo rf_2: reg del above c.fico_c above#c.fico_c                            ///
    log_balance ltv hard refinance white broker                                ///
    if abs(fico_c) <= 40, cluster(area)

eststo rf_3: reg del above c.fico_c above#c.fico_c                            ///
    log_balance log_income ltv age hard refinance white broker gender          ///
    if abs(fico_c) <= 40, cluster(area)

* First stage F-statistics
foreach spec in 1 2 3 {
    if `spec' == 1 {
        reg rate above c.fico_c above#c.fico_c if abs(fico_c) <= 40, cluster(area)
    }
    else if `spec' == 2 {
        reg rate above c.fico_c above#c.fico_c                                ///
            log_balance ltv hard refinance white broker                        ///
            if abs(fico_c) <= 40, cluster(area)
    }
    else {
        reg rate above c.fico_c above#c.fico_c                                ///
            log_balance log_income ltv age hard refinance white broker gender  ///
            if abs(fico_c) <= 40, cluster(area)
    }
    scalar fs`spec' = e(F)
}
di "First-stage F-stats: " fs1 " / " fs2 " / " fs3

* 2SLS: rate instrumented by FICO≥620 → delinquency
eststo iv_1: ivregress 2sls del (rate = above) c.fico_c above#c.fico_c       ///
    if abs(fico_c) <= 40, cluster(area)

eststo iv_2: ivregress 2sls del (rate = above) c.fico_c above#c.fico_c       ///
    log_balance ltv hard refinance white broker                                ///
    if abs(fico_c) <= 40, cluster(area)

eststo iv_3: ivregress 2sls del (rate = above) c.fico_c above#c.fico_c       ///
    log_balance log_income ltv age hard refinance white broker gender          ///
    if abs(fico_c) <= 40, cluster(area)

* Bandwidth and polynomial robustness
eststo rob_rf_60:   reg del above c.fico_c above#c.fico_c                    ///
    log_balance ltv hard refinance white broker if abs(fico_c) <= 60, cluster(area)

eststo rob_rf_quad: reg del above c.fico_c above#c.fico_c                    ///
    c.fico_c#c.fico_c above#c.fico_c#c.fico_c                                ///
    log_balance ltv hard refinance white broker if abs(fico_c) <= 80, cluster(area)

eststo rob_iv_60:   ivregress 2sls del (rate = above) c.fico_c above#c.fico_c ///
    log_balance ltv hard refinance white broker if abs(fico_c) <= 60, cluster(area)

eststo rob_iv_quad: ivregress 2sls del (rate = above) c.fico_c above#c.fico_c ///
    c.fico_c#c.fico_c above#c.fico_c#c.fico_c                                 ///
    log_balance ltv hard refinance white broker if abs(fico_c) <= 80, cluster(area)

* Probit marginal effect at threshold (robustness to functional form)
gen above_d = above
probit del above_d c.fico_c log_balance ltv hard refinance white broker       ///
    if abs(fico_c) <= 40, cluster(area)
margins, dydx(above_d) atmeans

esttab rf_1 rf_2 rf_3 iv_1 iv_2 iv_3,                                        ///
    keep(above rate) b(4) se(4) star(* 0.1 ** 0.05 *** 0.01)                  ///
    mtitles("RF(1)" "RF(2)" "RF(3)" "IV(1)" "IV(2)" "IV(3)")                  ///
    scalars("N Observations" "r2 R-squared")                                   ///
    title("Table 4: Main Results — Reduced Form and 2SLS")                    ///
    note("BW=±40, linear polynomial, clustered SE by area.")

esttab rob_rf_60 rob_rf_quad rob_iv_60 rob_iv_quad,                           ///
    keep(above rate) b(4) se(4) star(* 0.1 ** 0.05 *** 0.01)                  ///
    mtitles("RF BW=60" "RF Quad BW=80" "IV BW=60" "IV Quad BW=80")            ///
    title("Robustness: Bandwidth and Polynomial Order")


********************************************************************************
* STEP 4: HETEROGENEITY ANALYSIS (Table 5)
********************************************************************************

estimates clear

* Bank vs. broker channel
eststo rf_bank:   reg del above c.fico_c above#c.fico_c                      ///
    log_balance ltv hard refinance white if abs(fico_c)<=40 & broker==0, cluster(area)
eststo rf_broker: reg del above c.fico_c above#c.fico_c                      ///
    log_balance ltv hard refinance white if abs(fico_c)<=40 & broker==1, cluster(area)
eststo iv_bank:   ivregress 2sls del (rate=above) c.fico_c above#c.fico_c   ///
    log_balance ltv hard refinance white if abs(fico_c)<=40 & broker==0, cluster(area)
eststo iv_broker: ivregress 2sls del (rate=above) c.fico_c above#c.fico_c   ///
    log_balance ltv hard refinance white if abs(fico_c)<=40 & broker==1, cluster(area)

* Formal interaction test
eststo interact: reg del above c.fico_c above#c.fico_c broker above#broker   ///
    log_balance ltv hard refinance white if abs(fico_c)<=40, cluster(area)

* Racial heterogeneity
eststo rf_nonwhite: reg del above c.fico_c above#c.fico_c                    ///
    log_balance ltv hard refinance broker if abs(fico_c)<=40 & white==0, cluster(area)
eststo rf_white:    reg del above c.fico_c above#c.fico_c                    ///
    log_balance ltv hard refinance broker if abs(fico_c)<=40 & white==1, cluster(area)
eststo iv_nonwhite: ivregress 2sls del (rate=above) c.fico_c above#c.fico_c ///
    log_balance ltv hard refinance broker if abs(fico_c)<=40 & white==0, cluster(area)
eststo iv_white:    ivregress 2sls del (rate=above) c.fico_c above#c.fico_c ///
    log_balance ltv hard refinance broker if abs(fico_c)<=40 & white==1, cluster(area)

* Refinance vs. purchase
eststo rf_norefi: reg del above c.fico_c above#c.fico_c                      ///
    log_balance ltv hard white broker if abs(fico_c)<=40 & refinance==0, cluster(area)
eststo rf_refi:   reg del above c.fico_c above#c.fico_c                      ///
    log_balance ltv hard white broker if abs(fico_c)<=40 & refinance==1, cluster(area)
eststo iv_norefi: ivregress 2sls del (rate=above) c.fico_c above#c.fico_c   ///
    log_balance ltv hard white broker if abs(fico_c)<=40 & refinance==0, cluster(area)
eststo iv_refi:   ivregress 2sls del (rate=above) c.fico_c above#c.fico_c   ///
    log_balance ltv hard white broker if abs(fico_c)<=40 & refinance==1, cluster(area)

esttab rf_bank rf_broker iv_bank iv_broker,                                   ///
    keep(above rate) b(4) se(4) star(* 0.1 ** 0.05 *** 0.01)                  ///
    mtitles("RF Bank" "RF Broker" "IV Bank" "IV Broker")                       ///
    title("Table 5a: Heterogeneity — Bank vs Broker")

esttab rf_nonwhite rf_white iv_nonwhite iv_white,                             ///
    keep(above rate) b(4) se(4) star(* 0.1 ** 0.05 *** 0.01)                  ///
    mtitles("RF Non-White" "RF White" "IV Non-White" "IV White")               ///
    title("Table 5b: Heterogeneity — Race")

esttab rf_norefi rf_refi iv_norefi iv_refi,                                   ///
    keep(above rate) b(4) se(4) star(* 0.1 ** 0.05 *** 0.01)                  ///
    mtitles("RF Purchase" "RF Refinance" "IV Purchase" "IV Refinance")         ///
    title("Table 5c: Heterogeneity — Refinance vs Purchase")

* Figure 3: Delinquency by FICO and broker channel
preserve
keep if abs(fico_c) <= 40
bysort fico_bin broker: egen mean_del_bin = mean(del)

twoway                                                                         ///
    (scatter mean_del_bin fico_bin if broker==0,                               ///
        mcolor(eltblue) msymbol(circle) msize(small))                         ///
    (scatter mean_del_bin fico_bin if broker==1,                               ///
        mcolor(navy) msymbol(triangle) msize(small))                           ///
    (lfit mean_del_bin fico_bin if broker==0 & above==0, lcolor(eltblue) lwidth(medthick)) ///
    (lfit mean_del_bin fico_bin if broker==0 & above==1, lcolor(eltblue) lwidth(medthick)) ///
    (lfit mean_del_bin fico_bin if broker==1 & above==0, lcolor(navy) lwidth(medthick)) ///
    (lfit mean_del_bin fico_bin if broker==1 & above==1, lcolor(navy) lwidth(medthick)), ///
    xline(620, lcolor(black) lpattern(dash))                                   ///
    legend(order(1 "Bank" 2 "Broker") position(1) ring(0) cols(1))            ///
    xtitle("FICO Score") ytitle("Mean Delinquency Rate")
graph export "fig3_delinquency_broker.png", replace width(1600)
restore


********************************************************************************
* STEP 5: BANDWIDTH SENSITIVITY (Figure 4)
********************************************************************************

matrix bw_results = J(8, 4, .)
local i = 1

foreach bw of numlist 25 30 35 40 45 50 55 60 {
    quietly ivregress 2sls del (rate = above) c.fico_c above#c.fico_c        ///
        log_balance ltv hard refinance white broker                            ///
        if abs(fico_c) <= `bw', cluster(area)
    matrix bw_results[`i', 1] = `bw'
    matrix bw_results[`i', 2] = _b[rate]
    matrix bw_results[`i', 3] = _b[rate] - 1.96 * _se[rate]
    matrix bw_results[`i', 4] = _b[rate] + 1.96 * _se[rate]
    local i = `i' + 1
}

svmat bw_results, names(bw_)
gen bw_sig = (bw_3 > 0) if !missing(bw_3)

twoway                                                                         ///
    (rarea bw_3 bw_4 bw_1, color(navy%15) lwidth(none))                       ///
    (line bw_2 bw_1, lcolor(navy) lwidth(medthick))                           ///
    (scatter bw_2 bw_1 if bw_sig==1, mcolor(navy) msymbol(circle) msize(medlarge)) ///
    (scatter bw_2 bw_1 if bw_sig==0, mcolor(white) msymbol(circle)            ///
        msize(medlarge) mlcolor(navy) mlwidth(medthick)),                      ///
    yline(0, lcolor(black) lpattern(dash) lwidth(medium))                     ///
    xline(40, lcolor(navy) lpattern(shortdash) lwidth(medium))                ///
    legend(order(1 "95% CI" 2 "2SLS estimate" 3 "Significant" 4 "Insignificant") ///
        position(11) ring(0) cols(1))                                          ///
    xtitle("Bandwidth (FICO points)") ytitle("2SLS Coefficient on Interest Rate")
graph export "fig4_bandwidth_sensitivity.png", replace width(1600)

drop bw_1 bw_2 bw_3 bw_4 bw_sig

* Local randomization tests (small window around threshold)
* net install rdlocrand, from(https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata) replace

rdrandinf del  fico, c(620) wl(615) wr(625)   // window ±5
rdrandinf del  fico, c(620) wl(610) wr(630)   // window ±10
rdrandinf rate fico, c(620) wl(615) wr(625)   // first stage, ±5
rdrandinf rate fico, c(620) wl(610) wr(630)   // first stage, ±10

********************************************************************************
* END OF FILE
********************************************************************************
