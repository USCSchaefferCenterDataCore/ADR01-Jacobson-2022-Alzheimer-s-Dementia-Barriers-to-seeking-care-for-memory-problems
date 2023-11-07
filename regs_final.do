* This STATA Do File Reads in Data Collected from the Understanding America Study (SURVEY 284)
* The data can be obtained from UAS: https://uasdata.usc.edu/index.php

global datadir "/Users/mireillejacobson/Dropbox/Aging Schaeffer/AD projects/cog dissonance/"
global logdir "/Users/mireillejacobson/Dropbox/Aging Schaeffer/AD projects/cog dissonance/log/"
global gphdir "/Users/mireillejacobson/Dropbox/Aging Schaeffer/AD projects/cog dissonance/gph/"
global outdir "/Users/mireillejacobson/Dropbox/Aging Schaeffer/AD projects/cog dissonance/out/"

capture log close
use "${datadir}uas284_weights.dta", clear

svyset uashhid  [pw=final_weight]

drop if age<65

gen eval0 = .
gen eval1 = sc001 ==1 | sc004==1
gen eval2 = sc002 ==1 | sc005==1
gen eval3 = sc003 ==1 | sc006==1

replace eval0 = eval1 + eval2 + eval3
su eval0 eval1 - eval3

gen own = sc_ran==2

** RANDOMIZATION STRATA
gen age_order = age*10 + sc_order

drop if age_order==.

** PHQ-2 REVERSED CODED 
replace ba006a = 4-ba006a 
replace ba006b = 4-ba006b
gen phq2 = ba006a + ba006b

** CONSOLIDATE EDUCATION CODES - GROUP SEP BELOW HS, ASSOCS, AND MASTERS/PHD
tab education
replace education = 8 if education<8
replace education = 12 if education==11
replace education = 14 if education ==15 | education==16

gen college = education>=11
replace college = . if education==.

** CONSOLIDATE ASIAN AND PI
replace race = 4 if race==5

gen nonwhite = race~=1
replace nonwhite = . if race>4

** AGE GROUPS **
gen age_grp = 1 if age>=65 & age<75
replace age_grp = 2 if age>=75 & age<85
replace age_grp = 3 if age>=85


** CONSOLIDATE HH INCOME
replace hhincome = 5 if hhincome<5
replace hhincome=6 if hhincome==7
replace hhincome=8 if hhincome==9
replace hhincome=10 if hhincome==11
replace hhincome=12 if hhincome==13
replace hhincome=14 if hhincome==15

grstyle init
grstyle set plain

matrix drop _all

forvalues i=0(1)3 {
matrix beta_tx`i' = J(1, 2,.)
matrix colnames beta_tx`i' = Other Self
matrix ci_tx`i' = J(2, 2, .) 
matrix colnames ci_tx`i' = ll95 ul9
}


forvalues i=0(1)3 {
  reg eval`i' if own==0, cluster(age_order)
  matrix beta_tx`i'[1, 1] =  _b[_cons] 
  matrix ci_tx`i'[1, 1] = _b[_cons]  - invttail(e(df_r),0.025)*_se[_cons] 
  matrix ci_tx`i'[2, 1] = _b[_cons]  + invttail(e(df_r),0.025)*_se[_cons] 
  
  reg eval`i' if own==1, cluster(age_order)
  matrix beta_tx`i'[1, 2] =  _b[_cons] 
  matrix ci_tx`i'[1, 2] = _b[_cons]  - invttail(e(df_r),0.025)*_se[_cons] 
  matrix ci_tx`i'[2, 2] = _b[_cons]  + invttail(e(df_r),0.025)*_se[_cons]
}

coefplot (matrix(beta_tx0), ci(ci_tx0) label(Eval0)) ,  recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs5) xlabel(0(1)3) title("Vignette Assessment Index", size(medium)) saving("${gphdir}eval0", replace)
coefplot (matrix(beta_tx1), ci(ci_tx1) label(Eval1)) ,  recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs5) xlabel(0(0.25)1) title("Vignette 1 Assessments", size(medium)) saving("${gphdir}eval1", replace)
coefplot (matrix(beta_tx2), ci(ci_tx2) label(Eval1)) ,  recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs5) xlabel(0(0.25)1) title("Vignette 2 Assessments", size(medium)) saving("${gphdir}eval2", replace)
coefplot (matrix(beta_tx3), ci(ci_tx3) label(Eval1)) ,  recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs5) xlabel(0(0.25)1) title("Vignette 3 Assessments", size(medium)) saving("${gphdir}eval3", replace)

graph combine "${gphdir}eval1.gph" "${gphdir}eval2.gph" "${gphdir}eval3.gph", saving("${gphdir}eval_panel", replace)

* PLOT FOR FIGURE 1**
graph combine  "${gphdir}eval1.gph" "${gphdir}eval2.gph" "${gphdir}eval3.gph" "${gphdir}eval0.gph", saving("${gphdir}eval_panelv2", replace)


** PLOTS FOR HETEROGENEITY **


** RUN REGS BY CHARACTERISTIC **

foreach Y in eval0 eval1 eval2 eval3 { 
matrix beta_`Y' = J(1, 7,.)

matrix ci_`Y' = J(2, 7, .)
 
matrix colnames beta_`Y' = All Female Male "No College" College "Ages 65-74" "Ages 75 plus" 


matrix colnames ci_`Y' = All Female Male "No College" College "Ages 65-74" "Ages 75 plus" 

matrix rownames ci_`Y' = ll95 ul9

}

reg phq2 own , cluster(age_order) 
outreg2 using "${outdir}heterog", se ci noaster nolabel bdec(3) rdec(3) replace
  
   foreach Y in eval0 eval1 eval2 eval3 {
local i 0

foreach Z in "" "if gender==0" "if gender" "if college==0" "if college"  "if age_grp==1" "if age_grp==2 | age_grp==3"   {
  local ++ i
  reg `Y' own `Z', cluster(age_order) 
  outreg2 using "${outdir}heterog", se ci noaster nolabel bdec(3) rdec(3) append
  matrix beta_`Y'[1, `i'] = _b[own] 
  matrix ci_`Y'[1, `i'] = _b[own]  - invttail(e(df_r),0.025)*_se[own] \ _b[own]  + invttail(e(df_r),0.025)*_se[own]
   }
}


foreach Y in eval0 eval1 eval2 eval3 {
 local all_`Y' = beta_`Y'[1,1]
  }

 
* PLOT FOR FIGURE 2 
coefplot (matrix(beta_eval0), ci(ci_eval0) label(Index)) ,  xline(0, lcolor(black) lwidth(thin)) xline(`all_eval0', lpattern(dash) lcolor(navy) lwidth(vthin)) recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs10) title("Heterogeneity in Index-measure of Cognitive Dissonance", size(medium)) saving("${gphdir}eval0_heterog", replace)
coefplot (matrix(beta_eval1), ci(ci_eval1) label(Vignette 1)) , recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs10) title("Heterogeneity in Cognitive Dissonance and Vignette 1", size(medium)) xline(0, lcolor(black) lwidth(thin)) xline(`all_eval1', lpattern(dash) lcolor(navy) lwidth(vthin)) xlabel(-0.08(0.02)0)  saving("${gphdir}eval1_heterog", replace)
coefplot (matrix(beta_eval2), ci(ci_eval2) label(Vignette 2)),  xline(0, lcolor(black) lwidth(thin)) xline(`all_eval2', lpattern(dash) lcolor(navy) lwidth(vthin)) xlabel(0(0.1)0.4) recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs10) title("Heterogeneity in Cognitive Dissonance and Vignette 2", size(medium)) saving("${gphdir}eval2_heterog", replace)
coefplot (matrix(beta_eval3), ci(ci_eval3) label(Vignette 3)),  xline(0, lcolor(black) lwidth(thin)) xline(`all_eval3', lpattern(dash) lcolor(navy) lwidth(vthin)) xlabel(-.2(0.05)0) recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs10) title("Heterogeneity in Cognitive Dissonance and Vignette 3", size(medium)) saving("${gphdir}eval3_heterog", replace)



** PLOTS FOR APPENDIX

forvalues i=0(1)3 {
  reg eval`i' [pw=final_] if own==0, cluster(age_order)
  matrix beta_tx`i'[1, 1] =  _b[_cons] 
  matrix ci_tx`i'[1, 1] = _b[_cons]  - invttail(e(df_r),0.025)*_se[_cons] 
  matrix ci_tx`i'[2, 1] = _b[_cons]  + invttail(e(df_r),0.025)*_se[_cons] 
  
  reg eval`i' [pw=final_] if own==1, cluster(age_order)
  matrix beta_tx`i'[1, 2] =  _b[_cons] 
  matrix ci_tx`i'[1, 2] = _b[_cons]  - invttail(e(df_r),0.025)*_se[_cons] 
  matrix ci_tx`i'[2, 2] = _b[_cons]  + invttail(e(df_r),0.025)*_se[_cons]
}

coefplot (matrix(beta_tx0), ci(ci_tx0) label(Eval0)) ,  recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs5) xlabel(0(1)3) title("Vignette Assessment Index", size(medium)) saving("${gphdir}app_eval0", replace)
coefplot (matrix(beta_tx1), ci(ci_tx1) label(Eval1)) ,  recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs5) xlabel(0(0.25)1) title("Vignette 1 Assessments", size(medium)) saving("${gphdir}app_eval1", replace)
coefplot (matrix(beta_tx2), ci(ci_tx2) label(Eval1)) ,  recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs5) xlabel(0(0.25)1) title("Vignette 2 Assessments", size(medium)) saving("${gphdir}app_eval2", replace)
coefplot (matrix(beta_tx3), ci(ci_tx3) label(Eval1)) ,  recast(bar) ciopts(recast(rcap)) citop barwidt(0.3) bcolor(gs5) xlabel(0(0.25)1) title("Vignette 3 Assessments", size(medium)) saving("${gphdir}app_eval3", replace)

graph combine "${gphdir}app_eval1.gph" "${gphdir}app_eval2.gph" "${gphdir}app_eval3.gph", saving("${gphdir}app_eval_panel", replace)
graph combine "${gphdir}app_eval1.gph" "${gphdir}app_eval2.gph" "${gphdir}app_eval3.gph" "${gphdir}app_eval0.gph", saving("${gphdir}app_eval_panelv2", replace)

qui tab age_grp, gen(iage)
qui tab race, gen(irace)
qui tab education, gen(iedu)
qui tab hhinc, gen(ihhinc)
rename he003s1 Medicare


reg own gender 
outreg2 using "${outdir}maintab", se ci noaster nolabel bdec(3) rdec(3) replace

reg own gender [aw=final]
outreg2 using "${outdir}apptab", se ci noaster nolabel bdec(3) rdec(3) replace

log using "${logdir}regs.log", replace
** FIRST NO WEIGHTING

for var iage* irace* hisplatino iedu* gender ihhinc* retired Medicare: tab X
 
su age iage* irace* hisplatino iedu* gender ihhinc* retired Medicare  if own==0
su age iage* irace* hisplatino iedu* gender ihhinc* retired Medicare if own==1

foreach Y in age hisplatino gender retired Medicare {
ttest `Y', by(own)
}

foreach Y in age age_grp race hisplatino education gender hhinc retired Medicare  {
tab `Y' own, chi
}


forvalues i = 0(1)3 {
reg eval`i' own, cluster(age_order)
outreg2 using "${outdir}maintab", se ci noaster nolabel bdec(3) rdec(3) append
reg eval`i' own i.age_order, cluster(age_order)  
outreg2 using "${outdir}maintab", se ci noaster nolabel bdec(3) rdec(3) append
reg eval`i' own i.age_order i.race hisplatino i.educa i.gen i.hhinc retired he002, cluster(age_order)  
outreg2 using "${outdir}maintab", se ci noaster nolabel bdec(3) rdec(3) append
 }

 
** HETEROGENEITY ESTIMATES **
gen own_male = own * gender
gen own_nocollege = own * (college==0)
gen nocollege = college ==0

gen age75plus = (age_grp==2 | age_grp==3)
gen own_age75plus = own * age75plus

foreach Z in "own_male gender" "own_nocollege nocollege" "own_age75plus age75plus" {
reg eval0 own `Z' i.age_order, cluster(age_order) 
outreg2 using "${outdir}heterog", se ci noaster nolabel bdec(3) rdec(3) append
}
 
** APPENDIX TABLES AND FIGUES
su age iage* irace* hisplatino iedu* gender ihhinc* retired Medicare [w=final_] if own==0
su age iage* irace* hisplatino iedu* gender ihhinc* retired Medicare  [w=final_] if own==1

foreach Y in age age_grp race hisplatino education gender hhinc retired Medicare  {
svy: tab `Y' own
}


forvalues i = 0(1)3 {
reg eval`i' own [pw=final_], cluster(age_order)
outreg2 using "${outdir}apptab", se ci noaster nolabel bdec(3) rdec(3) append
reg eval`i' own i.age_order [pw=final_], cluster(age_order)  
outreg2 using "${outdir}apptab", se ci noaster nolabel bdec(3) rdec(3) append
reg eval`i' own i.age_order i.race hisplatino i.educa i.gen i.hhinc retired he002 [pw=final_], cluster(age_order)  
outreg2 using "${outdir}apptab", se ci noaster nolabel bdec(3) rdec(3) append
 }
