*! _esreg_getest.ado 0.4.4  12sep2026  A. Araar
*! Internal to the esreg family: selects the estimation the post-commands work from
*! (est(name) > current e() if esreg > _esreg).  In its own file because a command
*! called from another ado must be findable by its own name on the adopath.

cap program drop _esreg_getest
* select the estimation to work from: est(name) > current e() if esreg > _esreg
* usage:  _esreg_getest [name]   (call after -_estimates hold-, see post-commands)
program define _esreg_getest
    version 16
    args name
    if ("`name'" != "") {
        cap estimates restore `name'
        if (_rc) {
            di as err "stored estimation `name' not found"
            exit 301
        }
    }
    else if ("`e(cmd)'" != "esreg") {
        cap estimates restore _esreg
        if (_rc) {
            di as err "no esreg estimation in memory; run esreg first or give est(name)"
            exit 301
        }
        di as txt "(using the last esreg estimation, _esreg)"
    }
    if ("`e(cmd)'" != "esreg") {
        di as err "`name' is not an esreg estimation"
        exit 301
    }
    * the data in memory must still be the estimation sample
    qui count if e(sample)
    if (r(N) != e(N)) {
        di as err "the data in memory do not match the estimation sample (e(sample) has " r(N) " obs, e(N) = " e(N) ")"
        exit 459
    }
end
