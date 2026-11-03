# Q121: all-order radial tuple bound with the axis primitive

Parent: `3cf9920bdda832a8de0715724f75653f735a274c` (Q120).
On 3 November 2026 the freshly fetched remote agreed with that parent.
Q120 research CI 41608567542 and field CI 41608567714 were in progress.
Q118 field CI 41594347819 had passed; its research CI remained in progress.
Current notes, Q120 frozen evidence and its dependencies were inspected.
The original core_growth_activation result remains a leading-core and first
activation-collar estimate; width C_setup^-20 and historical evidence are
unchanged. No actual scale, coefficient, pressure gauge or cutoff is altered.

Q120 leaves a full tuple bound M_(n,m) in its existing-scale inequality.
core_background_radial_reconstruction supplied an explicit recipe only for
m=0,1,2. Here is a closed estimate at EVERY fixed m, from raw finalized F,U
jets. It includes the axis primitive and the extra pure-angular derivative.
Actual profile inputs and all existing scale/selector conditions remain
premises. This is an analytical implication, not a measured profile bound.

## 1. Same finalized fields and whole primitive interval

Use 0<h<1/100, n>=0, D=1/2-h, c=D+2nh, and

    R=sqrt(2X), L=1-2h eta^2, d=1-eta^2,
    M(X,eta)=integral_0^X U(Y,eta)dY,
    V=[2eta X U-2c eta M-d M_eta]/L,
    E=R F, W=V/R, J=2eta M/(LR).                     (1)

All are the SAME globally finalized coefficient n, including its full
pressure and exact moment repair. A velocity plateau does not reset M.
No pressure support identity follows from (1). Assume the existing positive
annulus Rlo^2/2<=X<=Xhi, Rhi^2>=2Xhi, |eta|<=1. Let Fraw bound all RAW
mixed X,eta derivatives of F through total order m on this annulus. Let
Uraw bound all raw mixed derivatives of U through order m+1 on the WHOLE
interval [0,Xhi], including the established extension beyond Xv if used.
The assumptions require C^m F and C^(m+1) U, not analytic profiles.
Domain coverage of the band chart is still the old condition
Rlo<=Rmin/2, Xhi>=2Rmax^2. Native and band radii are not interchangeable.

Put B=max(1,Xhi), ell=1-2h. For all r,k needed below,

    M_(0,k) = integral_0^X partial_eta^k U,
    M_(r,k) = partial_X^(r-1) partial_eta^k U, r>=1.

Thus raw derivatives of M and M_eta through order m are bounded by B Uraw.
For M_eta the r=0,k=m entry needs the pure (m+1)st angular derivative of U;
for r>=1 it needs at most total order m of U. This explains the one angular
loss. It cannot be inferred from U jets restricted to the observation
annulus. F needs no order m+1 input.

## 2. Closed positive real-jet majorant

Introduce nonnegative formal increments x,z. Collapse both to t only AFTER
forming normalized Taylor coefficients. The raw input hypotheses majorize
the finite jets of F,U,M,M_eta by, respectively,

    Fraw exp(2t), Uraw exp(2t), B Uraw exp(2t),
    B Uraw exp(2t),                                 (2)

through order m. Only finite coefficients enter; exp(2t) is a majorizing
series and does not assert convergence of the profiles' Taylor expansions.
Set

    tau=min(1,Rlo^2)/8.                              (3)

Then tau<=1/8 and 2tau/Rlo^2<=1/4. The actual geometric factors have
positive coefficient majorants, on evaluating at tau, satisfying

    Rhat<=2Rhi, Ihat<=2/Rlo, Lhat<=2/ell,
    etahat<=2, dhat<=2, Xhat<=Xhi+1<=2B,
    exp(2tau)<2.                                   (4)

For Rhat use the absolute binomial coefficients of sqrt(2X0+2x), with
constant Rhi and j>=1 coefficients evaluated at Rlo. Their sum is
Rhi+Rlo[1-sqrt(1-2tau/Rlo^2)]<=Rhi+Rlo<=2Rhi.
For Ihat, the binomial sum is Rlo^-1(1-2tau/Rlo^2)^(-1/2)<2/Rlo.
For Lhat use [ell-4ht-2ht^2]^-1; at t<=1/8 and h<1/100,
4ht+2ht^2<ell/2, giving (4). etahat=1+t, dhat=1+2t+t^2,
Xhat=Xhi+t. These are bounds at every real center in the retained hull,
including eta endpoints by the one-sided/extension convention already used.

Combining (1)--(4) gives positive majorant sums

    Ehat <=4Rhi Fraw, Uhat<=2Uraw,
    Vhat <=(40+16c)B Uraw/ell,
    What <=(80+32c)B Uraw/(ell Rlo),
    Jhat <=32B Uraw/(ell Rlo).

In particular their componentwise maximum is bounded by the convenient sum

    T=4Rhi Fraw + [2+(112+32c)B/(ell Rlo)]Uraw.      (5)

If r+k<=m, positivity bounds the (r,k) normalized coefficient by
T tau^(-r-k). Multiplying by r!k! and using tau<=1 and r!k!<=m! proves

    max_(component in E,U,W,J; r+k<=m)
      |partial_X^r partial_eta^k component|
       <= m! tau^(-m) T =: Mclosed_(n,m).           (6)

This estimates every fixed order, including m=0, uniformly over the same
hull whenever the actual Fraw,Uraw inputs are uniform. It introduces no
physical q power and makes no uniform-in-order bound claim. The old
anisotropic low-order estimates can be sharper and remain valid unchanged.
The historical factorial C_eta norm is NOT the raw norm in this theorem:
a row |partial_X^r U|_(C_eta^(m+1-r))<=u_r only gives k!u_r for the raw
(r,k) jet. These factorials must be included before using (6).

## 3. Join to Q120 without evaluating a giant factorial

For m>=1, m!<=m^m. A padded binary upper bound for (6) is

    Btuple = ceil(log2 max(1,T))
             + m ceil(log2 m) + m ceil(log2(1/tau)),        (7)

with the middle term zero for m=0. No factorial or high power need be
evaluated to compute (7). Feeding Btuple into Q120 gives the explicit
sufficient retained-scale condition, for 0<=m<=n,

    nh ell_n >= n+1+b_(m+1)
                   +m ceil(log2(1+2nh))+Btuple,            (8)

where ell_n<=log2(a_n) is a certified lower bound for the SAME existing
a_n and b is Q120's implemented Borel Euler cutoff binary budget.
Actual cutoff identity is still required. The finite head n<m in Q119's
background norm uses (6) directly; it is not covered by a fictitious m<=n
row. Its actual low-index high-jet profile bounds remain necessary.

The function diagonal_from_profiles only composes these two sufficient
inequalities. The supplied profile constants and log2 interval are not
certified by the function. The order-100010 diagnostic is manufactured,
not an actual native source or infinite-schedule certificate. A failed
sufficient inequality is not a counterexample to the actual field.
The original symbolic schedule prescription may admit constants of this
size, but this checkpoint does not reselect a_n or claim that an earlier
choice absorbs them. Other weighted-stress, support and selector gates
are untouched. Complete pressure memory is retained through the finalized
F,U inputs; its validity, support and moment identities are not derived
from a velocity norm.

## 4. Checks and limits

Exact rational positive-series checks compare the finite bivariate
reconstruction through order eight with (5)--(6), retaining inverse-radius
and L derivatives. Independent symbolic polynomial profiles and their
axis-to-X primitives check actual mixed jets through order four, including
eta=+-1 and coefficient n>0. Plateau controls with U=0 locally and
M=mu eta^k/k! give partial_eta^(k-1)W(X,0)=-mu/R for k>=1, so incoming
primitive memory is not lost. These controls are norm/identity tests,
not actual PDE coefficients or continuum proof by sampling.
The audit replays Q120 and its dependency/source hashes; new evidence is
immutable. Independent replay must match byte-for-byte. The CI workflow
runs this suite and replay in addition to its inherited checks.

Next prove the actual whole-interval finalized F,U raw bounds needed by
(6), on the required hull and at each coefficient/order, including finite
heads. Verify their binary budgets in (8) for the retained schedule, and
carry these actual constants through Q118/Q119 into the high-pole source
and full-pressure compact-modulation/exact-repair cone join. Pressure and
moment/support identities, whole-wave mean corrections including axial
velocity and pressure, uniform actual pulses over every physical scale,
and the complete viscous/cutoff/interaction residual remain open.
No unconditional pulse bound, full residual improvement or blow-up follows.
