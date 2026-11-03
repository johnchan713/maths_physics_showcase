# Q134: real endpoint heat jets and the same compensation root

Q133 joins the finalized angular U0 estimates; Q125 still needs the whole
F0 input. Q76 already provides finite-order heat tables and Q78 retains
the real Gamma formula at the poles. This checkpoint replaces their
order-dependent heat table by a closed real Gevrey-2 envelope, integrates
the COMPLETE infinite heat discrepancy, and propagates it through the
SAME accepted three-row relative heat root. It does not repeat the local
core/activation checkpoint or assert the whole F0 bound.

All conclusions below are analytical implications of the original reference
schedule, real Gamma representation, attenuation and accepted heat root in
[Q76](../core_current_head_q76_heat_compensation_jets/README.md),
[root foundation](../root_foundation/README.md), and
[global swirl envelope](../global_swirl_envelope/README.md). Actual continuum
profile binding and all physical pulse gates remain separate obligations.
No historical frozen record is changed.

## 1. Real endpoint estimate; no complex Gamma disk

Let 0<h<=1/100, X>=100, eta in [-1,1], Z=2(1-eta²)/X and

    H(Z) = integral_0^infty exp(-v) v^h (1+Zv)^(-h) dv / Gamma(1+h).

Differentiation at nonnegative Z gives
|H^(a)(Z)|/a! <= (h)_a(1+h)_a/a! <= h(a+1)! for a>=1.
Also |H(Z)-1|<=4h/X. For the finite angular jet at eta, the absolute
increment is bounded coefficientwise by J(q)/X, J(q)=4q+2q².
Only a<=m contributes at positive degree m. Since X^(-a)<=X^-1 and
J(1/16)=33/128<1/2, the positive series sum_(a>=1) J(q)^a has coefficient
at degree m less than 16^m. Thus the normalized angular coefficient is
at most (h/X)(m+1)!16^m. Multiplying by m! and using m+1<=2^m yields,
including m=0 by padding,

    |d_eta^m(H-1)| <= (4h/X) 32^m (m!)².                  (1)

This uses a finite chain rule for each m, uniformly up to one-sided pole
limits. It never evaluates H at Z<0. The absolute Z-Taylor coefficients
at zero are A_m=(h)_m(1+h)_m/m!, with
A_(m+1)/A_m=(m+h)(m+1+h)/(m+1) tending to infinity. For h>0 their
radius is zero. A common complex Gamma disk cannot be substituted here;
Q126's axial analytic disk did not supply such a disk.

For real Gevrey-2 products, raw caps P A^m(m!)² and Q B^m(m!)² give
PQ(A+B)^m(m!)² by Leibniz and 1/binom(m,j)<=1. Put delta=chi(H-1),
with the original radial cutoff 0<=chi<=1 independent of eta. Then
P_heat=(1+delta)²-1 obeys

    |d_eta^m P_heat| <= (16h/X)64^m(m!)²,                (2)

since 2(4h/X)+(4h/X)² <=16h/X. The value bounds do not assume
analyticity of chi or of the slow angular variable.

## 2. Full tail, correct signs and polynomial normalizers

On the ENTIRE edited future tail X=XK*x, x>=1, the same angularly uniform
pre-heat reference satisfies |Eref|<=2eK*x^(-1/2-h), XK>=Xstar>=100.
Its terminal taper is included in this estimate. With increments defined
as new minus reference, integrate the original X densities
E²/(2X), -E²/2 and sqrt(2X)E. From (1),(2), for every m>=0,

    |d_eta^m Delta Cp| <= (16h eK²/XK) 64^m(m!)²,
    |d_eta^m Delta S | <= (32h eK²)    64^m(m!)²,
    |d_eta^m Delta I | <= (8sqrt(2)eK sqrt(XK))64^m(m!)². (3)

The respective majorant integrals are 1/(2+2h), 1/(1+2h), and 1/h.
These integrable derivative caps justify differentiation of the full
renormalized-I difference for each fixed h>0 and each finite order.
In particular h integral_1^infty x^(-1-h) dx=1: the last row has no
extra h factor. Setting h=0 before infinity, or replacing the complete
future by a fixed finite radius, loses this contribution.

The three target rows are (-Delta Cp, Delta S, -Delta I), divided by
(e_star² f², Xstar e_star² f², sqrt(2)Xstar^(3/2)e_star f), f=(1+eta²)^-1.
The signs and half-factors remain those of the original three-row map.
The normalized polynomial coefficient sums of f^-1 and f^-2 at any
eta in [-1,1] are at most 5 and 25. Multiplication costs these constants
without increasing base 64: (m-j)!/m!<=1 and 64^(-j)<=1.
Using rH=eK sqrt(XK)/(e_star sqrt(Xstar))<=h^4 and XK>=Xstar gives

    row raw prefactors <= (400h^9, 800h^9, 40h^4)/Xstar,
    |d_eta^m z| <= Z L^m(m!)²,
    Z=800g, L=64, g=h^4/Xstar<=10^-10.                  (4)

These are complete terminal DIFFERENCES. The fixed incoming pressure is
common, not bounded by (3). Locally Pi(X)=Pi(Xstar)+integral E²/(2x)dx
with the actually edited E. The compensation restores totals only when
combined with the FUTURE heat edit. It neither erases Pi_eta nor resets
any of the genuine five incoming traces; M,J being unaffected by this
three-row relative edit does not make their incoming values zero.

## 3. Closed bound on the same accepted real root

Use Q76's conservative beta=50, q=100 and retained |c|<=6400g for
Bc+Q(c,c)=z. B,Q are constant in eta on this specific leading relative
patch. The pointwise Jacobian inverse is bounded by gamma=100 because
2 beta q (6400g) <=.0064<1/2. This is the original accepted root,
not the five-row compact repair and not a new all-order contraction.

For d_m=sup|c^(m)|/m!, differentiation gives

    d_m <= gamma[Z L^m m! + q sum_(j=1)^(m-1) d_j d_(m-j)].

Set v_m=d_m/m!. Since j!(m-j)!/m!<=1, these coefficients are dominated
by the analytic MAJORANT V=a Lx/(1-Lx)+b V², where a=100 Z=80000g,
b=10000. On |x|<=1/(8L), the forcing is <=a/7. The ball |V|<=a
is invariant and contractive: ab<=.08<1/8, so its image is <=a/7+a/8<a
and its Lipschitz constant <=1/4. Cauchy bounds the majorant coefficients,
not the physical slow-angle root. Therefore

    |c|<=6400g,
    |d_eta^m c| <=80000g*512^m*(m!)²,  m>=1.             (5)

The constant base is valid on the whole stated h,Xstar range; no
infinite family of C^k acceptance gates was added. Original positive
field, support and profile-existence premises are still required.

On the original reserved patch 0<=t=log(X/Xstar)<=20,
E=e_star f exp((-1/2-lambda)t)(1+sum_(j=1)^3 c_j beta_j(t)), U=0.
The beta_j are the ORIGINAL width-.3 unnormalized sigma' bumps, with
|beta_j|<=64. With |f^(m)|<=2*2^m*m! from the original radius-1/2 disks,
(5) yields the conservative correction bound

    |d_eta^m(E-Eref)| <=30720000 e_star g*514^m*(m!)².    (6)

The sum of three bump values was retained; disjointness could sharpen
it. Since sqrt(2X)>=1 and is independent of eta, (6) also bounds the
F0 correction on THIS patch. It gives no radial derivative or moving
support gain. Q76's finite mixed budgets and original supports remain.

## 4. Verification and remaining inequality

The checks use exact rational endpoint chains, the infinite moment
integrals, polynomial normalization, weighted quadratic recurrences,
a manufactured scalar root, and source/parent-evidence replay. They
are checks of the stated formulas, not a continuum proof. The frozen
record binds Q133 and the retained Q76/Q78/root/Gamma sources.

The next obligation is to assemble all exterior F0 stages (including
the axial pulse swirl, later angular edit, terminal heat multiplier and
taper) with the same reference and accepted roots, then join the actual
axis and finalized compact covers. For Q125 the still-needed inequality is

    sup_(0<=X<=Xvel) sum_(j=0)^m sup_eta |d_eta^j F0(X,eta)|/j! <= F_m,

on the SAME whole support used with Q133's U0 bound through m+2.
At the axis F0 must be estimated as the regular field, not by dividing
an E bound by sqrt(2X). Moreover Q121 still needs mixed X,eta derivatives
through m for F0 and m+1 for U0 on the whole physical hull; pure angular
bounds cannot supply them. No unconditional cone/pulse bound, residual
improvement, scale-uniform construction or blow-up is established.

Run test_audit.py and audit.py --verify-record evidence.json --output a
fresh path. Historical evidence remains frozen.
