* sec7_union3.do -- Section 7 of the paper: the manual's example (union3) read with the
* toolkit.  Produces the figures paper/fig_union3_curve.pdf and paper/fig_union3_mte.pdf
* and the reading sheet paper/reading_union3.txt; every number of the section is here.
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p esrdiag esrtest esrcurve esrmte _esreg_pwr esrreport {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}
* fonts of the paper (Computer Modern is not installed on Windows; Times is the closest)
cap graph set window fontface "Times New Roman"
cap graph set window fontfacesans "Times New Roman"
cap graph set window fontfaceserif "Times New Roman"
use "data/union3.dta", clear
summarize wage ln_wage union south black tenure age grade smsa if !missing(wage, union, tenure)
* ---- 0. the manual's reading: etregress, constrained model, wage in levels ----------------
etregress wage age grade smsa black tenure, treat(union = south black tenure) nolog
* ---- 1. the selection equation (probit gamma = two-step route) --------------------------
esreg ln_wage age grade smsa black tenure, select(union = south black tenure) method(twostep) nolog
esrdiag
esrtest, spec
esrtest, spec order(3)
* ---- 2. selection on gains ----------------------------------------------------------------
esrtest, pdid
esrtest, pdid nq(4)
esreg ln_wage age grade smsa black tenure, select(union = south black tenure) method(twostep) kappa(grade) nolog
esrtest, kappa
* ---- 3. the joint law: levels against logs --------------------------------------------
esreg wage age grade smsa black tenure, select(union = south black tenure) method(twostep) nolog
esrtest, normal
esreg ln_wage age grade smsa black tenure, select(union = south black tenure) method(twostep) nolog
esrtest, normal
* the FIML on the wage in levels (Table: the levels column) and its gamma contrast
esreg wage age grade smsa black tenure, select(union = south black tenure) nolog
esrtest, spec
* ---- 4. effects, profile, MTE: the two routes side by side ---------------------------------
esreg ln_wage age grade smsa black tenure, select(union = south black tenure) nolog store(fiml)
predict double P, pr
esrcurve, rank(P) nq(5)
esrmte, semipar
esreg ln_wage age grade smsa black tenure, select(union = south black tenure) method(twostep) nolog store(ts)
esrcurve, rank(P) nq(5) graph saving("paper/fig_union3_curve.gph", replace) ///
    title("Effect of union membership on ln(wage) by quintile of P(Z)")
graph export "paper/fig_union3_curve.pdf", replace
esrmte, semipar graph saving("paper/fig_union3_mte.gph", replace)
graph export "paper/fig_union3_mte.pdf", replace
matrix list r(mte)
matrix list r(mte_sp)
esrmte, est(fiml)
* ---- 5. the reading ----------------------------------------------------------------------
esrreport, export("paper/reading_union3.txt") replace
