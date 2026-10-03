* bf2_setup.do <case> -- the data and the estimation of each case of the extended
* brute force (bf2.do, ifcode.do): the first 1,000 observations of
*   ts   : data/esr_sim.dta, esreg y x, select(d = x z) method(twostep)
*   fiml : data/esr_sim.dta, esreg y x, select(d = x z)                (FIML)
*   kx   : data/esr_kx.dta,  esreg y x, select(d = x z) method(twostep) kappa(x)
*   het  : data/esr_kx.dta,  esreg y x, select(d = x z) hetsigma(x) hetrho(x)
*   hm   : data/esr_nonlin.dta, esreg y x, select(d = x z) method(twostep) hermite(3)
*   hr   : data/esr_nonlin.dta, esreg y x, select(d = x z) method(twostep) hermite(3 0)
*          (the Hermite terms in the treated regime only)
* Sets the global BF2OPT (the options of esreg) and creates id and wt (= 1).
args case
* the folder that holds data/ and tests/ (the working directory)
local P = subinstr("`c(pwd)'", "\", "/", .)
if inlist("`case'", "ts", "fiml") use "`P'/data/esr_sim.dta", clear
else if inlist("`case'", "kx", "het") use "`P'/data/esr_kx.dta", clear
else if inlist("`case'", "hm", "hr") use "`P'/data/esr_nonlin.dta", clear
else {
    di as err "case: ts, fiml, kx, het, hm or hr"
    exit 198
}
keep in 1/1000
gen long id = _n
gen double wt = 1
if ("`case'" == "ts")   global BF2OPT "method(twostep)"
if ("`case'" == "fiml") global BF2OPT ""
if ("`case'" == "kx")   global BF2OPT "method(twostep) kappa(x)"
if ("`case'" == "het")  global BF2OPT "hetsigma(x) hetrho(x)"
if ("`case'" == "hm")   global BF2OPT "method(twostep) hermite(3)"
if ("`case'" == "hr")   global BF2OPT "method(twostep) hermite(3 0)"
