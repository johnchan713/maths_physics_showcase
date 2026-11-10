# Q140: four PREMOD mixed fields bind the six scalar loop inputs

Parent: Q139, `90f910a0eb4ebbbb40bad8b6e831d222c8672539`.
A fresh branch read on 10 November 2026 agreed with that head. Its field
workflow 41771514554 and research workflow 41771514582 both completed
successfully. Q139's frozen evidence, Q65/Q66/Q128, and the actual regional
PREMOD estimates Q69--Q73/Q127 were read first. Historical evidence, the
selected cutoff, signed root, pressure datum, frequency, supports, scales
and `C_setup^-20` width are unchanged. The setup activation result remains
local to the leading core and first activation collar.

This checkpoint removes four of Q139's six independent mixed scalar inputs:
once genuine mixed Gevrey-3 bounds for the actual PREMOD fields
`a,b,p2,E` are supplied on a region, the original cutoff algebra gives
`t_s`, the same signed root and the augmented variance with explicit common
constants. It also identifies why the published Q72 short-taper cap cannot
be used directly for that supply. This is a correction to the proof route,
not a claim that the actual short-taper profile is non-Gevrey.

## 1. Finite two-variable Gevrey composition

At every real slow center suppose the SAME four PREMOD fields satisfy

    sum_(f in {a,b,p2,E}) |partial_y^i partial_eta^k f|/(i! k!)
        <= H Ay^i Be^k (i! k!)^2,                         (1)

with `H,Ay,Be>=1`. The pressure is the complete Q67/Q128 `p2`, including
the incoming `Pi`, `Pi_eta` and all five primitive traces. Assume the
original real gates

    0<a_*<=a<=A_*, |b|<=B_*, delta>0,
    |mu|,|p2|,|mu p2|<=G, G>=2, Eexp>=exp(G), 0<d0<=1.   (2)

No slow complex extension is assumed. The following finite composition
lemma will be used repeatedly. Let `F=(F_l)` obey (1), and at its real
value let an outer scalar map obey

    |partial^beta g|/beta! <= Q C^n (n!)^2, n=|beta|>=1. (3)

Then its real finite slow jets obey

    |partial_y^i partial_eta^k g(F)|
      <= Q [8Ay(1+CH)]^i [8Be(1+CH)]^k (i!k!)^3.         (4)

Indeed, for an ordered decomposition of `(i,k)` into `n` nonzero pairs,

    n! product_l (i_l! k_l!) <= (i+k)! <=2^(i+k)i!k!.

There are at most `2^(i+k) binom(i+k-1,n-1)` such decompositions. After
squaring the factorial inequality, summing over `n`, and using the sum in
(1), the coefficient is bounded by
`8^(i+k)(1+CH)^(i+k)(i!k!)^2`; restoring `i!k!` gives (4). This argument
uses only the downward-closed finite jet set requested at a given order.

Products use the similarly finite bound

    P A^i B^k(i!k!)^3 * Q C^i D^k(i!k!)^3
      -> P Q [2 max(A,C)]^i [2 max(B,D)]^k(i!k!)^3.      (5)

The arrows in (4),(5) mean bounds after the exact chain/Leibniz rules;
they are not multiplication of the displayed derivative values.

## 2. The same cutoff and signed root

On the scalar disk

    r=min(1/4,a_*/4),

the maps `t_s=-b/a` and `v0=a+b^2/a` have the safe caps

    T=2(B_*+1)/a_*,
    V=A_*+1+2(B_*+1)^2/a_*.                            (6)

Cauchy plus (4) gives both maps the common bases

    A0=8Ay(1+H/r), B0=8Be(1+H/r).                       (7)

Use exactly Q66's cutoff

    s=K(8(v0-2)/delta-1) sqrt((2+delta/2-v0)/a).        (8)

Where it is nonzero, the original real gap is `g=delta/4`. On
`q=min(1/4,a_*/4,g/4)` the analytic radicand factor has cap
`Cr=2 max(1,(2+delta/2)/a_*)`. Q66 gives
`|K^(n)|/n! <=729*8^n(n!)^2`. Apply (4) separately to the radicand and
to `K`, then (5). The implementation records

    Ms=729 Cr,
    As=2 max( 8A0(1+64V/delta),
              8 max(Ay,A0)(1+(H+V)/q) ),                (9)

and the analogous `Bs`. At both flat support boundaries the finite chain
rule agrees with the zero extension, so (9) holds without evaluating the
square root where its radicand may be negative. It does not assert an
analytic slow extension of `K`.

Q65's unchanged signed scalar inverse supplies `a_Psi,sigma_Psi>0` and

    |partial_p^u partial_s^v Psi|/(u!v!)
        <= a_Psi (2/sigma_Psi)^(u+v), u+v>=1.            (10)

Applying (4) to `(p2,s)` gives explicit mixed bases for the SAME accepted
`mu`, including zero pressure and zero variance. Its real zeroth-order cap
remains `G`; no division by `p2`, `s`, or `mu` occurs.

Finally set

    v=v0+a s^2.                                          (11)

Two applications of (5) close the augmented variance. Taking the sum of
the six component prefactors and the maximum of their two bases yields
explicit `Hsix,Asix,Bsix` such that

    sum_(lambda in {a,t_s,p2,mu,E,v})
      |partial_y^i partial_eta^k lambda|/(i!k!)
        <= Hsix Asix^i Bsix^k (i!k!)^2.                  (12)

Thus (12) is exactly Q139 equation (1). The helper immediately passes it
to Q139's fixed-origin primitive and Q138's original `Ny` pullback. This
is conditional on (1),(2); accepting numerical caps does not certify them
on the actual profile.

## 3. Precise obstruction in the inherited short-taper bound

The regional literature does not yet provide (1) with order-independent
constants. In Q72 the normalized order-`m` step norm is bounded by

    T_m=1+sum_(j=1)^m 162*3260^j*(j!)^2/t^j,
    t=C_setup^-20,                                      (13)

and the actual-field cap is then formed as `F0 exp(c T_m)` for a fixed
positive `c`. For fixed `H0,A0`, compression of this *cap sequence* into
`H0 A0^m(m!)^3` would require

    c*162*3260^m*(m!)^2/t^m
       <= log H0 + m log A0 + 3 log(m!).                 (14)

The left side divided by the right side tends to infinity. Therefore for
every fixed `H0,A0`, (14) fails at all sufficiently large `m`. Q72's
published finite-order upper bound cannot itself discharge Q139's uniform
mixed Gevrey hypothesis. A huge finite cap is not the same as a common
all-order Gevrey cap.

This does **not** show that the actual Q72 field violates a Gevrey bound.
The loss comes from exponentiating the entire order-`m` norm. The needed
repair is to differentiate `F_y=F (log F)_y` at the actual value and use a
positive Gevrey recurrence, as Q69/Q70 already do in their regions. That
recurrence must retain the original five incoming moments and both pressure
terms before the regional maximum is taken.

## 4. Checks, scope and next obligation

Exact controls test the factorial/decomposition inequality, the same
augmented variance identity, Q139/Q138 stage binding, binary bounds at
order 100005 and a moderate explicit separation of (13) from a fixed
Gevrey-3 budget. Frozen parent evidence and every source binding are
replayed; a separate CLI replay must match the record byte for byte.
These manufactured algebra checks and passing CI are not a continuum proof.

Next derive the value-based differentiated recurrence on Q72's short
taper, assemble one actual four-field `(H,Ay,Be)` from Q69--Q73/Q127 on
every original compact subregion, and verify the real gates (2). Then
apply this checkpoint, Q139, Q138 and Q137 and join Q121's whole hull,
including the entire axis primitive. Coefficient/reserved-mean supports,
complete wave-covariance corrections including axial velocity and pressure,
uniform actual physical pulse size/derivative/support bounds, and the full
viscous/cutoff/all-interaction residual remain open. No unconditional pulse
bound, complete residual improvement or blow-up is claimed.
