# Q141: value-based mixed Gevrey recurrence on the Q72 short taper

Parent: Q140, `2cff8d39aca13598b4034e55d40fe3d8dd7fcac5`.
A fresh branch read on 10 November 2026 agreed with that head; its two new
workflows were still in progress with no failed step. Q140's frozen record,
Q72's exact short-taper map, Q67's full-memory state formulas and Q63's
fixed step were read first. Historical evidence, the cutoff, signed root,
pressure datum, frequency, supports, scales and `C_setup^-20` width are
unchanged. The setup activation conclusion remains local to the leading
core and first activation collar.

Q140 proved that Q72's published cap `exp(c T_m)` grows too quickly to be
compressed into fixed Gevrey-3 constants. Here the same short-taper fields
are differentiated at their actual values. This removes that artificial
loss and produces a conditional all-order mixed envelope for the four
PREMOD fields required by Q140. The selected order-zero value, reciprocal,
geometry and real gate caps remain premises; the helper does not certify
them merely by accepting numbers.

## 1. A closed two-variable Gevrey algebra

Write `G(P,A,B)` for the raw real-jet estimate

    |partial_y^i partial_eta^k f|
       <=P A^i B^k (i! k!)^3.                            (1)

The finite Leibniz calculation used in Q140 gives

    G(P,A,B) G(Q,C,D)
       ->G(PQ,2 max(A,C),2 max(B,D)).                     (2)

A derivative costs

    partial_y:G(P,A,B)->G(PA,8A,B),
    partial_eta:G(P,A,B)->G(PB,A,8B),                    (3)

because `(i+1)^3<=8^i` for positive `i`, with equality handled separately
at zero. Sums add prefactors and take maximum bases. If `f>=f_*>0`, Q140's
finite scalar composition gives

    1/f ->G(2/f_*,8A(1+2P/f_*),8B(1+2P/f_*)).            (4)

No complex slow continuation is inferred by (4).

If `q` obeys (1), `g(eta)` is a constant-in-`y` trace with angular base
`B_g`, and `0<=y<=ell`, then

    g(eta)+integral_0^y q(s,eta) ds
      ->G(P_g+ell P,A,max(B,B_g)).                        (5)

For `i>=1`, (5) is the fundamental theorem of calculus; for `i=0` it is
the interval bound. This is the step Q72's full norm estimate obscured.

Finally suppose `L` obeys `G(H,A,B)` and an order-zero outward cap
`exp(plus/minus L)<=E_L` is known on the actual real family. Finite
composition at that value gives

    exp(plus/minus L)->G(E_L,8A(1+H),8B(1+H)).            (6)

Only the value cap is exponentiated. The derivative norm is not. Thus the
bases in (6) are fixed independently of requested order.

## 2. Rebuild the same reference fields

Keep Q72's variables

    t=C_setup^-20, 0<=y<=2t, X=Xa exp(y).

Its core polydisks give caps `L0,U0` and radii `h0,b` for `log Phi` and
`U_nat`. Cauchy and the shift by one `y` derivative give

    (log Phi)_y ->G(L0/h0,max(1,2/h0),max(1,1/b)),
    (U_nat)_y   ->G(U0/h0,max(1,2/h0),max(1,1/b)).        (7)

The SAME Q63 step at width `t`, including its flat endpoints, obeys

    psi,kappa ->G(162,max(1,3260/t),1).                   (8)

Use (2), integrate by (5), and apply (6) to the exact Q72 definitions

    J_r,y=psi (log Phi)_y,
    F_r=F(Xa) exp(J_r),
    U_r=U(Xa)+integral_0^y psi (U_nat)_y ds.              (9)

The initial `F`, `1/F`, `U` traces use their original analytic angular
radius `b`. Equations (7)--(9) give fixed mixed Gevrey-3 envelopes for
`F_r,1/F_r,U_r`; no `T_m` or order-dependent exponential occurs.

From these fields, integrate all five *reference* moment rows from their
same values at `Xa`:

    M_y=XU, I_y=2X^2F, J_y=2X^2FU,
    S_y=XU^2-X^2F^2, Pi_y=XF^2.                          (10)

Absolute envelopes replace the minus sign by a sum. The implementation
then evaluates the unchanged identities

    W=1-X^-1(2D eta M+(1-eta^2)M_eta),
    NQ=(1-h)I-D eta I_eta-(1-eta^2)J_eta+2(h-D)eta J,
    p1r=[-XW+NQ/(2XF)]/L,
    nsr=[-WU+X^-1{D(M-eta M_eta)+4h eta S
                  -(1-eta^2)S_eta}
              +4Aphys eta Pi-(1-eta^2)Pi_eta]/L.         (11)

Thus both reference sources have fixed mixed envelopes. All five incoming
traces remain: `I,J` occur in `p1r`, and both `Pi` and `Pi_eta` occur in
`nsr`. A local velocity bound cannot replace them.

## 3. Rebuild the actual fields and complete p2

Use the exact prescribed slopes, not the reference fields themselves:

    (log F_a)_y=-kappa p1r/2,
    (U_a)_y=-kappa X nsr/2.                              (12)

Apply (2),(5),(6) again from the original core traces at `Xa`. This gives
fixed envelopes for `F_a,1/F_a,U_a`. Reintegrate (10) with these ACTUAL
fields and the original incoming traces; the reference moments are not
copied. With `E=sqrt(2X)F_a`, construct

    a=1-2 E_y/E, b=2 U_y/E,
    w=X-2D eta M-(1-eta^2)M_eta,
    p2=[-wU+D(M-eta M_eta)+4h eta S-(1-eta^2)S_eta
        +X{4Aphys eta Pi-(1-eta^2)Pi_eta}]/(L E).         (13)

The positive algebra uses absolute values, so the signs in (11)--(13) do
not change the envelope. They remain written here to identify the original
map. Equations (1)--(13) yield explicit common `H,Ay,Be` for the actual
short-taper `a,b,p2,E`. Q141 passes them to Q140, then Q139's fixed-origin
primitive and Q138's original `Ny` pullback.

This is a conditional analytical estimate. Required inputs are the Q72
order-zero caps for `X,X^-1,sqrt(2X)`, the five traces, the two exponential
values in (9),(12), and the original real floors/value/frequency gates.
They are independent of derivative order and follow from the same selected
profile if their stated inequalities are checked. The code records them as
premises rather than manufacturing a numerical certification.

## 4. Checks and remaining obligation

Exact controls test the value-based exponential majorant, fixed binary
growth at order 100005, independent reintegration of pressure memory, the
nonzero `Pi_eta` response, and the Q140/Q139/Q138 binding. The old
`exp(cT_m)` lower budget is shown to exceed the new fixed-base budget in a
manufactured moderate case. Frozen provenance is replayed and a separate
CLI replay must reproduce the record byte for byte. These controls and CI
are not a continuum proof.

Next verify the actual order-zero Q72 caps and real gates, then compress
Q69--Q71/Q73/Q127 into common regional four-field constants and take their
maximum on the original compact cover. After that, join Q138/Q137 with
Q121 on the whole hull, including the entire axis primitive. Coefficient
and reserved-mean supports, complete wave-covariance corrections including
axial velocity and pressure, uniform actual physical pulse estimates, and
the complete viscous/cutoff/all-interaction residual remain unproved. No
unconditional pulse bound, complete residual improvement or blow-up is
claimed.
