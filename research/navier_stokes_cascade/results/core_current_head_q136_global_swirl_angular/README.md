# Q136: regular-axis, finalized compact and exterior F0 angular assembly

Parent: Q135, 3b1cbefc91263d851ae0bf27d31f4bf5ffb9c434.
The fresh remote branch agreed with that head. Both parent workflows were
in progress on 4 November 2026. Q135 supplied the whole exterior swirl
bound; this checkpoint supplies the missing conditional regular-axis and
FINAL compact F0 angular bounds, then joins the SAME spatial cover.
The unchanged C^-20 setup checkpoint remains local to the leading core
and first collar. All original supports, roots, N, gauges and widths remain.
Historical frozen evidence is preserved.

## 1. Estimate the regular field at the axis

[collar_operator_envelope](../collar_operator_envelope/README.md), (3),(4),
bounds the actual regular F_a=F_r exp(J) as well as U_a. Its common complex
cap K(ell) controls BOTH fields and their first two X derivatives. It uses
the retained squared-flat activation, not an ordinary exponential cutoff.
As in [global_axial_angular_bound](../global_axial_angular_bound/README.md),
section 4, use ell_axis=tc/2<t/20 and 4 exp(ell_axis)<4.1. The identical
bounds for J and its derivatives and the zero-free core strip therefore
apply through the flat modulation ENTRY. There is no change of cutoff
support and no radial analyticity assertion.

Take the same delta=rho_axis/(800S) and K=K(ell_axis) from the explicit
collar formula (4), with the original regular core bound inside Xa.
At every real X in [0,axis_end], Cauchy on the same eta disks gives

    |d_eta^j F0| <= K delta^-j j!,
    A_m = K sum_(j=0)^m delta^-j.                         (1)

This repeats the analytic argument at arbitrary finite angular order,
not the activation or local growth calculation. In particular at X=0 it
estimates F0 directly. E0/sqrt(2X) is NOT a permissible axis estimate.
The actual strip, K expression and flat-entry identification remain
retained continuum premises; entering a test value of K does not prove them.
The binary bound for K<=2^bK is bK+ceilLog(m+1)+m ceilLog(1/delta).

## 2. Final compact field, with the full repair included

On the SAME compact cover [xmin,xmax], xmin>0, Q131 gives the modulated
raw E cap Pmod Bmod^j(j!)^3. Q132 (8),(9) gives the field correction cap
Prep Brep^j(j!)^3 B0, where B0=162*32600 is Q63's conservative value
bound for the ORIGINAL width-.1 unit-mass bumps. Disjoint fixed supports
justify that field cap; no bump is renormalized or moved here. Thus

    PE = Pmod+B0 Prep, B=max(Bmod,Brep),
    |d_eta^j E_final| <= PE B^j(j!)^3.

Since sqrt(2X) is independent of eta and
(2xmin)^-1/2 <= max(1,xmin^-1), we obtain

    PF = max(1,xmin^-1) PE,
    |d_eta^j F_final| <= PF B^j(j!)^3,
    C_m = PF sum_(j=0)^m B^j(j!)².                       (2)

The estimate covers the original gaps/flat extensions only if the supplied
Q131 caps cover the WHOLE compact interval, as required already by Q133.
It does not turn earlier PREMOD data into final-profile data. The API
rejects PREMOD and requires Q131 modulation and Q132 final repair stages.

Crucially the Q132 root still requires all five actual baseline deficits
and all five incoming modulation memories, including Pi and Pi_eta.
Their target normalization, inverse-lambda loss and original accepted-root
and positive-field gates remain. This checkpoint bounds that same field
conditionally, not a repair chosen independently of the pressure target.
No terminal equality is inferred from four rows, sample angles or the fact
that modulation ended. No incoming pressure is set to zero.

## 3. Join the intended norm and its first-coefficient interface

Q135 supplies O_m=Pouter sum_(j=0)^m 548^j j! for all X>=Xstar.
Require the SAME leading profile and actual endpoints satisfying

    0<xmin<=axis_end, xmin<=xmax, xmax>=Xstar>=100.         (3)

Then [0,axis_end], [xmin,xmax] and [Xstar,infinity) cover the full support,
and Q125's intended pointwise-X angular algebra norm satisfies

    sup_X sum_(j=0)^m sup_eta |d_eta^j F0(X,eta)|/j!
        <= F_m := max(A_m,C_m,O_m).                      (4)

Overlaps cause no problem. A gap does: the API rejects axis_end<xmin or
xmax<Xstar. The exterior wrapper retains its actual Xstar; a mismatched exterior entrance is rejected even if the other intervals overlap. Norm order and binary/exact conventions must agree. Q133's
stronger sum-of-global-suprema warning is retained: (4) is not a claim that
regional maxima control that stronger norm. Binary regional envelopes
join by maximum with no extra sum cost. Compact binary cost is
ceilLog(max(1,PF))+ceilLog(m+1)+m ceilLog(max(1,B))+2m ceilLog(m).

The paired interface now supplies the two conditional high-angular inputs
of Q125: F0 through m and the finalized Q133 U0 through m+2. The shift
belongs to U0, not F0; the API rejects reversed orders or a PREMOD U bound.
Actual common profile/cover binding and all source/core/collar constants
are still required to apply the first-discrepancy estimate.

## 4. Checks and next obligation

Exact tests cover the singular-axis countercontrol, compact E-to-F loss,
repair addition, factorial normalization, overlaps and gap rejection,
order/convention guards and high-order binary budgets. A manufactured
smooth field checks the finite jets; it is not the selected physical
profile. The audit independently reproduces Q135 and checks inherited
source hashes; fresh replay must match the new record byte for byte.
Passing these checks is not a continuum proof.

The next substantive obligation is the actual mixed radial/eta tuple on
Q121's WHOLE physical hull: F0 through m and U0 through m+1, with the
real cutoffs, inverse widths and all repair radial derivatives. A bound
on pure eta derivatives such as (4) cannot provide it. Actual state and
low-order cone gates and the unchanged physical scales must also be bound.
Later coefficients/reserved means, uniform physical pulse estimates,
complete covariance corrections with axial velocity/pressure, and the
viscous/cutoff/all-interaction residual remain open. No unconditional
pulse bound, residual improvement or blow-up is established.

Run test_audit.py and audit.py --verify-record evidence.json --output a
fresh path. Historical evidence remains frozen.
