*! esreg_lf1 0.1.0  10sep2026  A. Araar
*! lf1 likelihood evaluator for esreg (Gaussian ESR, analytic score).
*! Equations: 1 y_1, 2 y_0, 3 selection, 4 lnsigma_1, 5 lnsigma_0,
*! 6 atanhrho_1, 7 atanhrho_0.  $ML_y1 = y, $ML_y2 = d.
cap program drop esreg_lf1
program define esreg_lf1
    version 16
    args todo b lnfj g1 g2 g3 g4 g5 g6 g7
    tempvar xb1 xb0 zg ls1 ls0 ar1 ar0
    mleval `xb1' = `b', eq(1)
    mleval `xb0' = `b', eq(2)
    mleval `zg'  = `b', eq(3)
    mleval `ls1' = `b', eq(4)
    mleval `ls0' = `b', eq(5)
    mleval `ar1' = `b', eq(6)
    mleval `ar0' = `b', eq(7)
    local y $ML_y1
    local d $ML_y2
    tempvar s1 s0 r1 r0 e1 e0 q1 q0 A1 A0 m1 m0
    quietly {
        gen double `s1' = exp(`ls1')
        gen double `s0' = exp(`ls0')
        gen double `r1' = tanh(`ar1')
        gen double `r0' = tanh(`ar0')
        gen double `e1' = (`y' - `xb1') / `s1'
        gen double `e0' = (`y' - `xb0') / `s0'
        gen double `q1' = sqrt(1 - `r1'^2)
        gen double `q0' = sqrt(1 - `r0'^2)
        gen double `A1' = (`zg' + `r1' * `e1') / `q1'
        gen double `A0' = (`zg' + `r0' * `e0') / `q0'
        replace `lnfj' = cond(`d' == 1, lnnormalden(`e1') - `ls1' + lnnormal(`A1'), ///
                                        lnnormalden(`e0') - `ls0' + lnnormal(-`A0'))
        if (`todo' == 0) exit
        gen double `m1' = exp(lnnormalden(`A1') - lnnormal(`A1'))
        gen double `m0' = exp(lnnormalden(`A0') - lnnormal(-`A0'))
        replace `g1' = cond(`d' == 1, (`e1' - `m1' * `r1' / `q1') / `s1', 0)
        replace `g2' = cond(`d' == 1, 0, (`e0' + `m0' * `r0' / `q0') / `s0')
        replace `g3' = cond(`d' == 1, `m1' / `q1', -`m0' / `q0')
        replace `g4' = cond(`d' == 1, `e1'^2 - 1 - `m1' * `r1' * `e1' / `q1', 0)
        replace `g5' = cond(`d' == 1, 0, `e0'^2 - 1 + `m0' * `r0' * `e0' / `q0')
        replace `g6' = cond(`d' == 1, `m1' * (`e1' + `r1' * `zg') / `q1', 0)
        replace `g7' = cond(`d' == 1, 0, -`m0' * (`e0' + `r0' * `zg') / `q0')
    }
end

