*! esreg_engine 0.5.0  11sep2026  A. Araar
*! Shared engine of the esreg family (esreg, predict, esrcurve, esrmte,
*! esrdiag, esrtest).
*!   Mata  : individual effects, weighted aggregation, delta-method variance,
*!           stacked-moment variance of the two-step route, cell means by group.
*!   esreg_engine compiles the Mata code of esreg_mata.ado on demand.  The Stata
*!   helpers of the family are in their own files, as a command called from another
*!   ado must be findable by its own name: _esreg_getest.ado, _esreg_data.ado,
*!   _esreg_effects_post.ado.
*! The Mata code lives in esreg_mata.ado; esreg_engine compiles it on demand.

* ===========================================================================
* Stata-side helpers of the family
* ===========================================================================
cap program drop esreg_engine
program define esreg_engine
    version 16
    * are the Mata functions there?  (an ado-file autoloaded from a program does not
    * execute its mata: block, and r() cannot be relied on across versions: st_local)
    mata: st_local("esr_ok", strofreal(findexternal("_esreg_ddvars()") != NULL & ///
                   findexternal("_esreg_effects()") != NULL & findexternal("_esreg_hausman()") != NULL))
    if ("`esr_ok'" != "1") {
        cap findfile esreg_mata.ado
        if (_rc) {
            cap findfile esreg_mata.ado, path("stata")
        }
        if (_rc) {
            di as err "esreg_mata.ado not found: reinstall the esreg package"
            exit 601
        }
        local fn "`r(fn)'"
        cap noisily version `c(stata_version)': run "`fn'"
        mata: st_local("esr_ok", strofreal(findexternal("_esreg_ddvars()") != NULL & ///
                       findexternal("_esreg_effects()") != NULL & findexternal("_esreg_hausman()") != NULL))
        if ("`esr_ok'" != "1") {
            di as err "the Mata functions of esreg_mata.ado could not be compiled"
            exit 601
        }
    }
end

