*! esreg_examples 1.0.0  03oct2026  A. Araar (Universite Laval / PEP)
*! The examples of help esreg, run from their links.
*!   esreg_examples #          run example # in the command window
*!   esreg_examples #, db      open the dialog box of esreg filled in for example #
*!   esreg_examples #, do      open example # as a do-file in the Do-file Editor
*! The data in memory are never lost: the run keeps them (preserve) and gives
*! them back at the end, even after an error or a Break; the do-file does the
*! same; the dialog box, which needs the example data in memory, refuses to
*! replace data that have unsaved changes.
*! Examples 1-3 use the manual's union3 data, read with webuse (an internet
*! connection is needed); examples 4-6 use simulated samples that are ancillary
*! files of the package, read from the current folder (where "ssc install
*! esreg, all" or "net get esreg" copies them), else from the SSC archive, else
*! from GitHub (_esreg_exload); nothing is written.
program define esreg_examples
    version 16
    syntax anything(name=ex id="example number") [, DB DO NOEDIT]
    capture confirm integer number `ex'
    if _rc | !inrange(`ex', 1, 6) {
        di as err "esreg_examples: the examples are numbered 1 to 6 (see help esreg)"
        exit 198
    }
    if "`db'" != "" & "`do'" != "" {
        di as err "esreg_examples: db and do cannot be combined"
        exit 198
    }
    local U  "ln_wage age grade smsa black tenure, select(union = south black tenure)"
    local S  "income educ, select(treatment = educ i.region inst)"
    local n 0
    if `ex' == 1 {
        local title "The manual's example: the switching regression by full-information ML"
        local data "union3"
        local c1 "esreg `U'"
        local n 1
    }
    else if `ex' == 2 {
        local title "The two-step, the diagnostics of the selection equation, the tests and the reading"
        local data "union3"
        local c1 "esreg `U' method(twostep)"
        local c2 "esrdiag"
        local c3 "esrtest, normal"
        local c4 "esrtest, pdid"
        local c5 "esrreport"
        local n 5
    }
    else if `ex' == 3 {
        local title "The effect along the score and the marginal treatment effect"
        local data "union3"
        local c1 "esreg `U' method(twostep)"
        local c2 "predict double P, pr"
        local c3 "esrcurve, rank(P) nq(5) graph"
        local c4 "esrmte, semipar graph"
        local n 4
    }
    else if `ex' == 4 {
        local title "Selection on gains that varies with x, kappa(x), and its test"
        local data "esr_kx"
        local c1 "esreg y x, select(d = x z) method(twostep) kappa(x)"
        local c2 "esrtest, kappa"
        local n 2
    }
    else if `ex' == 5 {
        local title "A conditional mean not linear in u: the Hermite test, then the augmented two-step"
        local data "esr_nonlin"
        local c1 "esreg y x, select(d = x z) method(twostep)"
        local c2 "esrtest, normal"
        local c3 "esreg y x, select(d = x z) method(twostep) hermite(3 0)"
        local c4 "esrmte, semipar graph"
        local n 4
    }
    else if `ex' == 6 {
        local title "A survey design: the two-step with vce(svy), the likelihood with the svy prefix"
        local data "esr_wt"
        local c1 "svyset [pweight = wt], strata(region)"
        local c2 "esreg `S' method(twostep) vce(svy)"
        local c3 "svy: esreg `S'"
        local c4 "esreg"
        local n 4
    }

    * ---- as a do-file, in Stata's temporary folder ----
    if "`do'" != "" {
        local T = c(tmpdir)
        local T : subinstr local T "\" "/", all
        if substr("`T'", -1, 1) != "/" local T "`T'/"
        local fn "`T'esreg_example_`ex'.do"
        tempname fh
        file open `fh' using "`fn'", write text replace
        file write `fh' "* esreg, example `ex': `title'" _n
        file write `fh' "* Written by esreg_examples in Stata's temporary folder; save it elsewhere to keep it." _n
        file write `fh' "* preserve keeps the data in memory and gives them back when this do-file ends;" _n
        file write `fh' "* delete that line to keep working on the example data." _n
        file write `fh' "preserve" _n
        if "`data'" == "union3" {
            file write `fh' "* the manual's data, from Stata's web site" _n
            file write `fh' "webuse union3, clear" _n
        }
        else {
            file write `fh' "* the example data: the current folder (where ssc install esreg, all copies" _n
            file write `fh' "* them), else the SSC archive, else GitHub" _n
            file write `fh' "capture use `data', clear" _n
            file write `fh' `"if _rc capture use "http://fmwww.bc.edu/repec/bocode/e/`data'.dta", clear"' _n
            file write `fh' `"if _rc use "https://raw.githubusercontent.com/aabbdd12/esreg/main/examples/`data'.dta", clear"' _n
        }
        forvalues i = 1/`n' {
            file write `fh' `"`c`i''"' _n
        }
        file close `fh'
        if "`noedit'" == "" doedit "`fn'"
        di as txt "(example `ex' written to " as res `"`fn'"' as txt ")"
        exit
    }

    * ---- in the dialog box: needs the example data in memory ----
    if "`db'" != "" {
        * the example data loaded for a dialog box carry a mark (_dta[esreg_example],
        * or the mark of the examples of another package): they are replaced without
        * asking; the user's own data with changes are not
        local isex 0
        local cl : char _dta[]
        foreach c of local cl {
            if substr("`c'", -8, .) == "_example" & `"`: char _dta[`c']'"' == "1" local isex 1
        }
        if c(changed) & !`isex' {
            di as err "esreg_examples, db: the data in memory have changes not saved;"
            di as err "save them (or clear) first: the dialog box needs the example data in memory"
            exit 4
        }
        _esreg_exload `data'
        char _dta[esreg_example] "1"
        if `ex' == 6 qui svyset [pweight = wt], strata(region)
        di as txt "(example data `data' loaded for the dialog box)"
        if `ex' == 2 di as txt "(the dialog box fits the model; then type esrdiag, esrtest, normal, esrtest, pdid and esrreport)"
        if `ex' == 3 di as txt "(the dialog box fits the model; then type predict double P, pr, esrcurve, rank(P) graph and esrmte, semipar graph)"
        if `ex' == 4 di as txt "(the dialog box fits the model; then type esrtest, kappa)"
        if `ex' == 5 di as txt "(the dialog box fits the augmented two-step; run the example to see the Hermite test first)"
        if `ex' == 6 di as txt "(data svyset; the dialog box fits the two-step with vce(svy); svy: esreg is typed in the command window)"
        db esreg
        * Stata keeps the state of the dialog between two openings: every control
        * an example may set is first put back to its default
        .esreg_dlg.main.vn_dv.setvalue ""
        .esreg_dlg.main.vl_iv.setvalue ""
        .esreg_dlg.main.vn_tv.setvalue ""
        .esreg_dlg.main.vl_zv.setvalue ""
        .esreg_dlg.main.vl_hsig.setvalue ""
        .esreg_dlg.main.vl_hrho.setvalue ""
        .esreg_dlg.main.vl_kap.setvalue ""
        .esreg_dlg.main.ck_herm.setoff
        .esreg_dlg.main.cb_herm.setvalue "3"
        .esreg_dlg.main.cb_hreg.setvalue "both"
        .esreg_dlg.main.vn_hs.setvalue ""
        .esreg_dlg.main.ed_store.setvalue ""
        .esreg_dlg.byifin.ck_by.setoff
        .esreg_dlg.byifin.vl_by.setvalue ""
        .esreg_dlg.byifin.ed_if.setvalue ""
        .esreg_dlg.byifin.ed_in.setvalue ""
        .esreg_dlg.wgt.rb_wnone.seton
        .esreg_dlg.wgt.ed_wexp.setvalue ""
        .esreg_dlg.se.rb_oim.seton
        .esreg_dlg.se.vn_clu.setvalue ""
        .esreg_dlg.rpt.ck_noeff.setoff
        .esreg_dlg.rpt.ck_nosvy.setoff
        .esreg_dlg.max.ck_nolog.setoff
        .esreg_dlg.max.ck_diff.setoff
        .esreg_dlg.max.ck_iter.setoff
        .esreg_dlg.max.ed_iter.setvalue ""
        * the example
        if inlist(`ex', 1, 2, 3) {
            .esreg_dlg.main.vn_dv.setvalue "ln_wage"
            .esreg_dlg.main.vl_iv.setvalue "age grade smsa black tenure"
            .esreg_dlg.main.vn_tv.setvalue "union"
            .esreg_dlg.main.vl_zv.setvalue "south black tenure"
        }
        if inlist(`ex', 4, 5) {
            .esreg_dlg.main.vn_dv.setvalue "y"
            .esreg_dlg.main.vl_iv.setvalue "x"
            .esreg_dlg.main.vn_tv.setvalue "d"
            .esreg_dlg.main.vl_zv.setvalue "x z"
        }
        if `ex' == 6 {
            .esreg_dlg.main.vn_dv.setvalue "income"
            .esreg_dlg.main.vl_iv.setvalue "educ"
            .esreg_dlg.main.vn_tv.setvalue "treatment"
            .esreg_dlg.main.vl_zv.setvalue "educ i.region inst"
        }
        if `ex' == 1 {
            .esreg_dlg.main.rb_ml.seton
        }
        else {
            * the two-step and the controls that go with it (enabled explicitly, in
            * case the radio button does not run its script when set from here)
            .esreg_dlg.main.rb_2s.seton
            .esreg_dlg.main.tx_hsig.disable
            .esreg_dlg.main.vl_hsig.disable
            .esreg_dlg.main.tx_hrho.disable
            .esreg_dlg.main.vl_hrho.disable
            .esreg_dlg.main.tx_kap.enable
            .esreg_dlg.main.vl_kap.enable
            .esreg_dlg.main.ck_herm.enable
        }
        if `ex' == 4 {
            .esreg_dlg.main.vl_kap.setvalue "x"
        }
        if `ex' == 5 {
            .esreg_dlg.main.ck_herm.seton
            .esreg_dlg.main.cb_herm.enable
            .esreg_dlg.main.tx_herm.enable
            .esreg_dlg.main.cb_hreg.enable
            .esreg_dlg.main.cb_herm.setvalue "3"
            .esreg_dlg.main.cb_hreg.setvalue "treated"
        }
        if `ex' == 6 {
            .esreg_dlg.se.rb_svy.seton
        }
        exit
    }

    * ---- in the command window: the data in memory are kept ----
    preserve
    _esreg_exload `data'
    local src "`r(source)'"
    di as txt _n "{hline 78}" _n "esreg, example `ex': " as res "`title'" _n as txt "{hline 78}"
    di as txt "(example data `data', read from `src'; the data in memory come back at the end)"
    forvalues i = 1/`n' {
        di as txt _n `". `c`i''"'
        `c`i''
    }
end

* ============================================================================
* the data of the examples: union3 from Stata's web site; the simulated samples,
* ancillary files of the package, from the current folder (where "ssc install
* esreg, all" or "net get esreg" copies them), else from the SSC archive, else
* from GitHub; nothing is written to disk
program define _esreg_exload, rclass
    args f
    if "`f'" == "union3" {
        capture quietly webuse union3, clear
        if !_rc {
            return local source "Stata's web site (webuse)"
            exit
        }
        di as err "esreg_examples: union3 could not be read with webuse (is the internet reachable?)"
        exit 631
    }
    capture confirm file "`f'.dta"
    if !_rc {
        quietly use "`f'.dta", clear
        return local source "the current folder"
        exit
    }
    capture quietly use "http://fmwww.bc.edu/repec/bocode/e/`f'.dta", clear
    if !_rc {
        return local source "the SSC archive"
        exit
    }
    capture quietly use "https://raw.githubusercontent.com/aabbdd12/esreg/main/examples/`f'.dta", clear
    if !_rc {
        return local source "GitHub"
        exit
    }
    di as err "esreg_examples: `f'.dta is not in the current folder (`c(pwd)'),"
    di as err "  and it could not be read from the SSC archive or from GitHub."
    di as txt "  Copy the example data into the current folder with"
    di as txt `"  {stata "ssc install esreg, all replace"} (from SSC), or"'
    di as txt `"  {stata "net get esreg, from(https://raw.githubusercontent.com/aabbdd12/esreg/main)"} (from GitHub)."'
    exit 601
end
