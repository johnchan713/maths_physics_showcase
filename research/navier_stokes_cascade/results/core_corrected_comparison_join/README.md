# Corrected diagonal field: pressure recovery and global comparison

Parent `40c2c5e56222ef49da8be912eeccc856bf49ae59` conditionally joins the
corrected diagonal field to zero initial data and one smooth compact force.
This checkpoint applies the ordinary expanding-ball comparison argument to
that **same** velocity, complete fixed-gauge pressure and force.  It does not
substitute an older localized family.

The conclusion is conditional on the retained field, force and positive-core
contracts.  The pressure and radius arguments below are analytical; the
exact controls only audit their exponents, signs, provenance and scope.

## 1. Difference equation and essential energy hypothesis

Fix `T<1`.  The parent supplies, conditionally,

    u_t+(u dot grad)u=Delta u-grad p+f,      div u=0,        (1)

where `u,p` are smooth and supported in one compact set `K` on `[0,T]`,
`u(0)=0`, and `f` is the single smooth compact global force.  The pressure
`p=cp0` retains the corrected field's complete fixed gauge.

Let `v,P` be any smooth solution on all space with the same zero datum and
the same force, and assume

    sup_(0<=t<=T) ||v(t)||_2 < infinity                     (2)

for each fixed `T<1`.  No global derivative or dissipation bound and no
pressure decay hypothesis is added.  Put

    w=v-u,       pi=P-p,       g=v tensor v-u tensor u.

Then `w` is uniformly `L2` and `g` uniformly `L1` on the fixed interval,
and the common force cancels exactly once:

    w_t+div g=Delta w-grad pi,       div w=0,
    -Delta pi=partial_i partial_j g_ij.                     (3)

Condition (2) is essential.  Without it, `v=t e1,P=-x1` and `u=p=0`
have the same zero datum and zero force but differ by a nonzero harmonic
pressure gradient.  Their energy grows like the volume of a ball.

## 2. Recover the pressure gradient without pressure decay

Let

    pi_star=sum_(i,j) R_i R_j g_ij,                         (4)

where the multiplier is `-xi_i xi_j/|xi|^2`.  Since `g` is `L1`, its
Fourier transform is bounded and `pi_star` lies uniformly in `H^-2`.
The strict threshold is `s>3/2`; equality is insufficient.

For a compact smooth time test `a`, (3) gives

    integral a grad pi dt
      =Delta integral a w dt+integral a' w dt
        -div integral a g dt.                              (5)

The right side belongs to `H^-3`.  Subtracting the corresponding average
of `grad pi_star` gives a harmonic `H^-3` distribution.  Its Fourier
transform is a weighted `L2` function supported at the single point zero,
so it vanishes.  Therefore

    grad pi=grad pi_star                                   (6)

as a space-time distribution.  This derives the needed temperedness of the
averaged gradient and does not assume that `P` itself is tempered.

Individually the pressure equations remain forced:

    -Delta p=DD(u tensor u)-div f,
    -Delta P=DD(v tensor v)-div f.                         (7)

Neither pressure may be replaced by an unforced Riesz formula.  The common
`div f` has already cancelled in (3) and must not be subtracted again.

## 3. Compact mixed pressure and the noncompact quadratic part

Write

    g=w tensor w+h,       h=w tensor u+u tensor w.          (8)

Because `u` is supported in `K subset B_L`, so is `h`.  On the cutoff
annulus `R<=|x|<=2R`, with `R>=max(1,2L)`, the double-Riesz kernel and
`|x-y|>=R/2` give

    |sum R_i R_j h_ij(x)| <= C H_T R^-3,                   (9)

where `H_T=sum ||h_ij||_1`.  Pairing with `R^-1 w` over the annulus and
using the global `L2` bound yields an `O(R^-5/2)` pressure flux.

No cancellation of the leading compact moment is assumed.  Divergence-free
compact fields can produce an anisotropic mixed moment, so isotropy does not
follow from incompressibility or support.  It is also unnecessary for (9).

For the quadratic part use `chi_R=phi(x/R)^8`, with `phi=1` on the unit
ball and zero outside radius two.  Set

    E_R=integral chi_R |w|^2,
    D_R^2=integral chi_R |grad w|^2,
    B_R=||phi_R^4 w||_6.                                  (10)

Local Sobolev gives

    B_R<=C(D_R+R^-1 W_T).                                 (11)

Split the localized quadratic pressure into its main double-Riesz term and
commutator.  The commutator kernel obeys

    C |z|^-3 min(|z|/R,1).                                (12)

At exponent `4/3`, its radial norm satisfies

    R^-4/3 integral_0^R r^-2/3 dr
      +integral_R^infinity r^-2 dr=4/R,                   (13)

so the `L^(4/3)` norm is `O(R^-3/4)`.  Both near cancellation and the far
cap are necessary.  This is not an `L1`-to-`L1` pressure estimate.

The resulting quadratic pressure and transport fluxes are

    C_T R^-1 B_R^(3/2),       C_T R^-7/4 B_R^(3/4).        (14)

The eighth cutoff power is sufficient for this pairing; seven is not.

## 4. Expanding-ball inequality

Pair the velocity form of (3) with `chi_R w`.  The exact identity is

    (1/2)E_R'+D_R^2
      =-integral chi_R (w dot grad)u dot w
       +(1/2)integral |w|^2 Delta chi_R
       +(1/2)integral |w|^2 v dot grad chi_R
       +integral pi w dot grad chi_R.                      (15)

For `R>=2L`, `u=0` and hence `v=w` on the cutoff annulus.  The strain term
is bounded by `||grad u||_infinity E_R`; the Laplacian term is `O(R^-2)`.
For `R^-b D_R^a`, `0<a<2`, Young absorption leaves radius power

    2b/(2-a).                                              (16)

The two powers in (14) give respectively `4` and `14/5`; their Sobolev
offsets and the compact mixed-pressure flux decay like `R^-5/2`.  The
coarse uniform-`L2` Laplacian term is the limiting `R^-2` term.  Therefore

    (1/2)E_R'+(1/2)D_R^2
      <=||grad u||_infinity E_R+C_T R^-2,      E_R(0)=0.   (17)

The `R^-2` rate is the scope of this coarse ledger, not a universal optimal
rate for each fixed solution.

## 5. Limits and conditional exclusion

First fix `T<1`.  With
`L_T=sup_[0,T] ||grad u||_infinity`, Gronwall gives

    E_R(t)<=2 C_T T exp(2 T L_T) R^-2.                     (18)

Let `R` tend to infinity.  Each fixed ball eventually lies where
`chi_R=1`, so `w=0` there.  Continuity gives `v=u` on `[0,T]`; since `T`
was arbitrary, equality holds for every `t<1`.  No monotone cutoff family
and no constants uniform as `T` approaches one are needed.

Only after this radius limit approach the parent's protected core path.
Conditionally, the localization plateau preserves

    |u_theta(x_tau,1-tau)| >= (e0/2) tau^-A,               (19)

with positive `e0,A`.  These points converge to `(0,1)`, while a global
smooth competitor is bounded on one compact neighborhood of that point.
Thus the same-force global smooth bounded-energy competitor is excluded
**if** all retained field, force and positive-core premises hold.

This closes the downstream comparison implication for the corrected field.
It is not an unconditional admissible-force construction or blow-up proof.
No future velocity or pressure is constructed.  The next obligation is an
independent foundation audit, beginning with identification of the exact
core pressure family with the completed global backward pressure, followed
by the retained background and finite defining/existence/inverse/gap
inequalities at the unchanged parameters.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
