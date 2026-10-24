# Symbolic first real signed covariance and fresh pressure-change join

Parent: `084a424c17fc348ce6ca98426ac26b71e39e0dc6`.
On 24 October 2026 its symbolic signed-root CI step passed; the full
run remained in progress without failures. The exact first real signed
covariance identity and its rates already exist in
`real_signed_covariance_update`. This checkpoint substitutes the recent
explicit primary/signed norm recipes, with symbolic constants, into all
six covariance entries and the complete pressure/scalar change formulas.
It enforces the NEW pressure-domain requirement. The actual projector
uniform seminorm and original background/root/value/support premises
remain OPEN; no physical family, root, band, width or gauge is changed.

## 1. Use the same actual roots, phase and common coefficient chart

Let A,B be the original principal and normalized full-curl coefficient
jets, and L,S the new signed principal and normalized full-curl jets
from the SAME primary U-chart root. For through-order n define

    a_n=16*188^n*n! sum_(j<=n) A_j((25/16)Sscale),

and identically b_n,l_n,s_n from B,L,S. These dominate raw derivative
maxima on the retained matched coarsest common chart. The factor16,
188^n, S_i<=(25/16)Sscale, normalized Haar measure and at-most-one
closed wave label are the existing atlas constants. No Eedge unweighting
is inserted; all coefficients still have sqrt(zeta) P times their edge
powers. For both original and signed physical normalization the amplitude
ratio factor is <=8: with Aexp=1/2+h and alpha=2-3kappa, all four
exponents Aexp-h b, b in {1/2,1-kappa,alpha-1/2,alpha-kappa}, lie
between0 and3/4. Their derivatives use the same fixed linear charts.

The current Xi is a PHYSICAL covariance target. For r=Q*/Qi>=1,
its per-label normalized coefficient has additional factor
r^(h*alpha-2*Aexp)=r^(-1-3h*kappa)<=1. The inverse slow coordinate
changes from a label to the common chart have derivative multipliers
<=1, and Xi is auxiliary independent. Also S*<=Si and every norm
polynomial has positive coefficients. Thus the common source bound
legitimately supplies the native per-label target input before the
atlas conversion; it is not silently copied between incompatible
physical normalizations.

For every label the OLD and NEW actual complete waves have phase-locked
real cosine principal and sine curl coefficients. Angular averaging
therefore gives the EXACT tensor identity

    Delta W=Lambda+Ecov,
    Lambda=1/2(A tensor L+L tensor A),
    Ecov=1/2(Bold tensor Snew+Snew tensor Bold
                               +L tensor L+Snew tensor Snew). (1)

Bold and Snew here denote the physical sine coefficients, with orders
1-kappa and alpha-kappa. All cutoff curl terms remain in them. The two
selected Haar entries of Lambda match Xi by the original fixed-root,
squared-partition identity; its other four entries need not vanish.
P^2<=1 can be dropped in upper bounds, but the full zeta is retained.

## 2. Explicit all-entry bounds and chart/pole bookkeeping

For the raw common-U tensor derivatives a safe sufficient coefficient is

    Lambda_n <= 2^(n+1) a_n l_n,
    Ecov_n <= 2^n(2 b_n s_n+l_n^2+s_n^2).             (2)

These intentionally overestimate angular halves and include Leibniz
costs. Apply the derived SAME-WEIGHT U-to-R tensor transfer Tr_n ONCE
to both expressions. Both then bound the historical raw maximum and
factorial through-order R norms. Use common wave pole

    q_n=max(6n+9,6n+d+6), d=target quotient allowance,

so every covariance term has full zeta and pole2q_n. This is a
finite requested-order estimate, not a uniform pole in all derivatives.
The three error orders in (1) are alpha+1-2kappa, 2alpha-1 and
2alpha-2kappa. With alpha=2-3kappa their minimum is nu=3-6kappa;
the signed square is the limiting, NONZERO term. The complete change
and Lambda have only order alpha. The historical real-phase assumption
is essential; a later cosine-quadrature error has a cross with L which
survives averaging, so this estimate is not extended to later complex
states.

## 3. Recompute the complete mean pressure and all scalar changes

Mean VELOCITIES are held fixed during this wave substep; no old mean
is reset. Let D_ab=bar(Delta W_ab). The full changes are

    dg=-D1 Drr-epsilon partial_Z Drz+Dtt/R,
    dP=int dg, dp=I0(dg-rho dP),
    dJtheta=int R^2 Dtheta,z,
    dJz=int R Dzz-(1/2)int R^2 dg.

These are changes to the full existing pressure, which still contains
Pi1/Pi2 and old mean/wave components. They are not replacement by a
leading pressure. Invoke the EXISTING `pressure_target_budget` with
the derived complete Delta coefficients at n,n+1,n+2 and error at n.
Its axis-pressure response keeps Wzz, Wtheta,z, rho dP and all required
derivatives. The whole integrated pressure contribution is

    epsilon partial_Z(dJz+c_rho dP),

including c_rho,Z dP. The code also ADDS these scalar changes to the
old order2 scalar coefficients to bound the actual updated defects.
They now have order alpha, not a claimed retained S2 order. The updated
projected targets and Haar tangential rows have the historical order
min(nu,alpha+1,3)=nu; the full pointwise pressure/rows can still have
the weaker alpha-kappa order. This is no full residual improvement.

## 4. Reject reuse of an insufficient pressure proof domain

Pressure through n+1 needs covariance through N=n+2, which needs the
signed target through n+3. The NEW pressure source pole must cover

    2q_N=2max(6N+9,6N+d+6).                          (3)

The API's pressure output n supplies TARGET coefficients through n;
recovering Haar ROW jets through n needs this target recipe at n+1.
Declare the finite maximum accordingly rather than differentiating a
bound which has no next input jet.

The old initial Haar pressure map was supplied at d=2p. It is generally
insufficient for (3). The operator callback MUST declare coverage of
the new pole and order and supply rho,Kp,Kmom, the original oriented
target primitive, projector bump and whole c_rho bounds there. A stale
old-domain declaration is rejected even if its numerical entries look
large. These declarations remain mathematical lower hypotheses, not
automatic certificates of the actual maps. The earlier exact symbolic
edge kernels give a conditional route at any fixed new pole, keeping
the original cutoff gap and supports; do not introduce a new S threshold
or infer a new-domain estimate from the old one.

## 5. Evidence and next obligation

Exact real angular tensor controls retain every signed square and curl
term and detect the later-complex obstruction. New substitutions verify
the explicit coefficient and pressure recipes, the enlarged domain,
all current scalar additions and moving pressure memory. Marked
manufactured jets keep diagnostics moderate; no giant continuum field
is evaluated. Historical frozen evidence is preserved and replayed.

Next supply the actual projector-family seminorm and the new-domain
fixed pressure/projector/primitive estimates, then apply fresh fast
corrections to the CURRENT post-signed residual and slow corrections
to CURRENT recomputed defects. Preserve signed squares and later complex
cross terms before the first sourced inverse. The original background,
phase, value/cone and support assumptions, uniform physical-scale
realization, all viscosity/cutoff/interactions and complete residual
remain unresolved. Unconditional pulse bounds, infinite realization
and blow-up are not proved.
