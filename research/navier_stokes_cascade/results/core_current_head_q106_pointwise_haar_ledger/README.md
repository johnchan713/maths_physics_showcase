# Q106: exact pointwise Haar rows and a moment-blind covariance direction

Parent: `595c9222fcc5c389dfe41a370870975ca26afd3c` (Q105).
Research CI 41494913919 and field-comparison CI 41494912987 were both
in progress when inspected on 2 November 2026. The remote head was Q105.
This checkpoint follows its pointwise Haar-row obligation. It supplies an
exact differential ledger and a precise missing norm inequality, not a
new fine epsilon estimate or a canonical-wave counterexample. Frozen
historical evidence, fields, root, accepted cover and width C_setup^-20
are unchanged.

## 1. Average the current pre-fast state, including products

Use precisely Q105's pre-fast state O: initial fast/slow plus first-real
fast/slow, together with the TOTAL old-plus-sourced complete wave W.
The new Q101 fast F and Q102 slow S have not yet been applied here.
Let a bar mean auxiliary Haar averaging, X=bar beta, U=bar v, Zm=bar gamma,
and Qab=bar(a b) for the FULL means. Qab is generally NOT bar a times bar b.
The original auxiliary-independent b,V,G, rho and moving chi are fixed.
Label parameters are fixed as in Q105. Set De=partial_R+e/R and
Dz=epsilon partial_Z. The exact Haar rows are

    Htheta = -epsilon partial_T U + D2 At + Dz Bt
                    -epsilon(partial_R^2+R^-1 partial_R-R^-2+Dz^2) U,
    Hz = -epsilon partial_T Zm + D1 Az + Dz Bz
                    -epsilon(D1 partial_R+Dz^2) Zm,                 (1)

    At=b U+V X+Q(beta,v)+bar Wrt-Sigma_theta,
    Bt=G U+V Zm+Q(gamma,v)+bar Wzt,
    Az=b Zm+G X+Q(beta,gamma)+bar Wrz-Sigma_z,
    Bz=2G Zm+Q(gamma,gamma)+bar Wzz+bar pm.                         (2)

Stress symbols here denote the Haar of the SAME assigned radial fluxes;
if their original construction is auxiliary-independent these are just
the original fluxes. There is no new freely chosen correcting stress.
The full graph contribution to Dr integrates to zero on the auxiliary
periodic torus. The R^-1 connections and both viscous pieces survive.
Time averaging removes c Nt but not -epsilon partial_T. All identities
require the inherited smoothness and common domain; closed cover and
physical derivative transfers are still separate obligations.

Only the fast means have zero Haar, not their correlations. Even when
X,U,Zm consist only of the two slow histories, Q(beta,v), Q(beta,gamma),
Q(gamma,v) and Q(gamma,gamma) must contain ALL four incoming mean stages.
Averaging products of Haar fields would omit precisely the fast/old
correlations retained in Q102. `rows` requires these averaged products
explicitly. It is a symbolic identity helper, NOT a certified native
provider or an outgoing fine-row callback.

## 2. Complete pressure survives averaging

Put

    S0=2b X+Q(beta,beta)+bar Wrr,
    H0=b Zm+G X+Q(beta,gamma)+bar Wrz,
    C0=2V U+Q(v,v)+bar Wtt,
    r0=epsilon partial_T X-S0/R-Dz H0+C0/R+epsilon Dz^2 X,
    K0=-S0+epsilon D1 X.                                          (3)

Translations on the full auxiliary torus preserve Haar measure. Therefore
averaging the SAME characteristic pressure T0 gives the radial moving
operator L0 f=int_a^R f dr-chi(R,Z,T) int_a^b f dr. The complete pressure
increment relative to the retained balanced base is

    bar pm=K0+L0(r0-rho Ppre),  Ppre=int_a^b r0 dr.                (4)

The equality for Ppre uses the whole compact K0 trace contract. With seam
or endpoint defects those terms must remain, exactly as in Q98. Ppre is
the Q99 pre-fast scalar, not the final Q103 scalar. Pi1/Pi2 and all earlier
pressure histories are retained by the original base subtraction and
Q105's accumulated state; no exterior/base residual is erased by (4).
In particular neither -S0 nor epsilon D1 X can be dropped from (4).
Dz in Bz differentiates the ENTIRE moving primitive, including chi_Z and
rho_Z Ppre+rho partial_Z Ppre. Integrating globally or knowing Ppre alone
does not bound that pointwise function. `radial_split` checks the local
identity bar gr=r0+partial_R K0; it does not evaluate a native L0 norm.

## 3. The exact missing inequality and full join

For each output n and a common proved norm/pole p, define H_n as the norm
of the two ACTUAL expressions (1), with (2)--(4) substituted. The missing
fine Haar estimate at requested exponent a is precisely

    sup_(accepted pulses, 0<epsilon<=epsilon0)
          epsilon^(-a) ||(Htheta,Hz)||_(n,p) < infinity.            (5)

The constant must be uniform on the original physical cover, with every
required multiplier, pressure, product and trace jet supplied. Q103's
three scalar coefficients do not establish (5). Nor is bounding each
flux separately at its coarse order a proof of the cancellations in (1).
A sufficient coefficient ledger, retaining the epsilon factors, is

    Ht_n <= epsilon(n+1) U_(n+1)
          +(n+1) At_(n+1)+2 r_n At_n
          +epsilon(n+1) Bt_(n+1)
          +epsilon[(n+1)(n+2) U_(n+2)+r_n(n+1) U_(n+1)
                     +r2_n U_n+epsilon^2(n+1)(n+2) U_(n+2)],
    Hz_n <= epsilon(n+1) Zm_(n+1)
          +(n+1) Az_(n+1)+r_n Az_n+epsilon(n+1) Bz_(n+1)
          +epsilon[(n+1)(n+2) Zm_(n+2)+r_n(n+1) Zm_(n+1)
                     +epsilon^2(n+1)(n+2) Zm_(n+2)].              (6)

Here r_n and r2_n are the factorial-through norms of R^-1 and R^-2 on
the SAME positive hull. These formulas follow from Leibniz and ordinary
coefficient derivative shifts; pole embeddings/products must be paid
before using them. They are analytical sufficient estimates, not evaluated
uniform constants. `envelope` implements (6) on nonnegative rational coefficient sequences,
requiring field jets through n+2 and flux jets through n+1. These are
conditional inputs, not a native factory or a proof of (5). `rows` computes
the exact symbolic expressions independently of that norm arithmetic.

Once Ht_n,Hz_n are actually supplied, Q105's exact identity gives FULL rows

    final_theta_n <= Ht_n+RFtheta_n+RStheta_n,
    final_z_n <= Hz_n+tail_n+RFz_n+RSz_n.                          (7)

There is no extra Qcirc factor on full RF or RS. The original A1 tail
still acts on Qcirc Epre_z; it is not replaced by the Haar row. `join`
implements only this positive arithmetic, not a current factory or a
proof of the missing input. Thus Q105's projected estimates alone still
do not close the full cone. Q104 radial A0 and scalar gates are unchanged.

## 4. All three global moments can miss a positive covariance direction

A sharper comparison than Q105's arbitrary zero-moment row is available.
On [1,2] take any nonzero h in C_c^infinity((1,2)), f=h^2, and zero
mean fields. For example h=exp(-1/((R-1)(2-R))) inside (1,2), with
zero extension, is flat at both endpoints: each fixed derivative is an
exponential times a rational function with finite-order endpoint poles,
and the exponential dominates every such pole. Consider the smooth
compact positive semidefinite covariance

    W=f [[1,1,0],[1,1,0],[0,0,1]]  in (r,theta,z) coordinates.    (8)

Its constant matrix has eigenvalues 0,1,2. It is the exact auxiliary
covariance of sqrt(2) h (cos y,cos y,sin y). This comparison velocity is
NOT asserted divergence-free or canonical and is not a counterexample
to the actual constrained construction. Its role is to test what follows
from positivity and moments ALONE, without those missing extra constraints.

For the fixed zero assigned stress in this comparison, (3)--(4) give
r0=0, K0=bar pm=-f. The exact Q102 covariance moment changes, including
pressure, are ALL zero:

    delta P=int (Wtt-Wrr)/R dr=0,
    delta Jtheta=int R Wzt dr=0,
    delta Jz=int R Wzz dr-(1/2) int R(Wrr+Wtt) dr=0.              (9)

Yet Htheta=D2 f=f'+2f/R is not zero; Hz=0 and
int R^2 Htheta dr=[R^2 f]_1^2=0. Scaling (8) by any positive constant
preserves all three zeros and positivity while scaling the pointwise row.
Consequently no inequality controlling this pointwise row solely by
P,Jtheta,Jz, with zero right-hand side at zero moments, follows on the
comparison class. The omitted axial pressure term would instead give a
spurious nonzero delta Jz; keeping it is essential to this control.
For the actual construction the assigned stress may cancel part of W;
that cancellation has to be proved in At,Az and Bz in (2), with its real
support. It cannot be inferred from (9). This is a rigorously documented
information gap, not a lower bound on the canonical residual.

## 5. Checks and next obligation

The executable comparison uses h=(R-1)^3(2-R)^3 on [1,2]. Its zero
extension is only C^2 (and f is C^5); it is a finite exact diagnostic,
not the globally smooth witness above. The analytical identities (8)--(9)
hold for every smooth compact h without evaluating its integrals.

Exact Fourier/polynomial checks compare full graph rows against (1),
detect replacing averaged products by products of averages, retain local
pressure memory and the moving primitive derivative, verify the radial
split and all three moment zeros for (8), and check its PSD spectrum and
Fourier realization. Frozen replay binds this checkpoint to Q105 and its
untouched recursive evidence. Tests are exact diagnostics of these finite
identities, not a continuum proof of native inputs or (5).

Next bind the ACTUAL averaged products and stress-subtracted fluxes in
(2) to the current target/root and obtain (5), or identify a canonical
uncancelled term. A second scalar reset cannot substitute for this work.
Then join (7) with Q104's radial gates and the inner/outer estimates.
Uniform physical-scale pulse estimates, reserved-mean supports, full
viscous/cutoff/interaction residual improvement, infinite realization and
blow-up remain unproved. No new fine Haar grade is claimed here.
