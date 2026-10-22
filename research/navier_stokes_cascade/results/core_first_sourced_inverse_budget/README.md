# First forced inverse: a Gaussian-adapted ambient jet budget

Parent: `a0bce249f120368f8437962aabaecb05de293264`.
The parent derives the weighted initial Haar rows feeding the original-root
signed target. This checkpoint makes the subsequent forced ODE's norm
recurrence explicit, including a moving-constraint issue that a direct
ambient differentiation must handle. It does NOT verify the current
nonzero-angular source inequality after all fresh mean updates.

The new estimate is conditional on the original exact core/background,
full-normal floor, tangent propagator and coefficient/source inputs below.
No profile, root, carrier, cutoff, accepted band, radius or width C^-20
changes. Historical evidence is preserved. The frozen controls are exact
manufactured algebra and finite polynomial bounds, not continuum estimates
of the selected physical fields.

## 1. The tangent propagator is not an ambient derivative estimate

Work on one common lifted path s in [a,b], with frozen discrete label,
epsilon, k and nonzero input harmonic m. Let x denote the remaining
continuous parameters. The entrance a is fixed in these coordinates;
moving entrance coordinates require their endpoint chain terms first.
Put N=|n|^2, h=n/N, lambda_m=m^2 epsilon k^2 N, and let

    A0=-K+h(n^T K-n'^T),  A=A0-lambda_m I,
    H=h n^T,  Pi=I-H,  t'=A t-Pi f,  t(a,x)=0.          (1)

A prime means path differentiation. This is exactly the projected
forced equation of the first-sourced foundation, and n.t=0 is preserved.
The physical pressure is still

    pi=i[(n^T K-n'^T)t+n.f]/(km N).                      (2)

No component of n is inverted; only the retained floor N>=c_N>0 is used.
For an x derivative, n.partial_x t=-partial_x n.t generally is NONZERO.
Consequently the tangent homogeneous propagator cannot simply be applied
to an arbitrary differentiated ambient equation.

There is a second pitfall. The natural ambient extension A has
(n.t)'=-lambda_m(n.t). That is not sufficient to preserve the pulse
weight P(s,x)>0. On a decaying portion with P=e^-2s and lambda_m=1,
the normalized normal response is e^s. Dropping P or replacing its
inverse by a polynomial S constant is unjustified. This is an obstruction
to that shortcut, not a refutation of the historical moving-frame proof.

Let p=P'/P, with all requested mixed coefficient jets bounded as inputs.
Use instead the ANALYTICAL AUXILIARY extension

    Ahat=A+(lambda_m+p)H.                                (3)

It agrees with A on tangent vectors, so changes no actual solution in
(1), source, pressure or operation. Its normal component z=n.t obeys
z'=p z. The weighted normal propagator is therefore exactly P(s)/P(w).
This extension is used only to estimate derivatives.

## 2. A polynomial ambient bound retaining the entire pulse weight

Assume the inherited tangent propagator, in the Euclidean physical
tangent norm (include frame conversion constants if starting in a frame),
obeys, for a<=w<=s<=b,

    ||V_T(s,w)|| <= T(S) P(s)/P(w),  b-a<=L(S).          (4)

T and L are nonnegative coefficient polynomials in S=ell^2, with T>=1.
Let H0,N0 bound |h|,|n| on the WHOLE path, and set

    q=A0 h-h',  Q0=||A0|| H0+||h'||,
    Kamb=T(1+L Q0 N0)+H0 N0.                            (5)

The identity n.q=0 follows by differentiating n.h=1. For arbitrary
ambient initial vector v, decompose t=u+h z, u tangent. Then

    z(s)=P(s)z(w)/P(w),
    u'=A u+q z,  u(w)=Pi(w)v.

Orthogonality gives |Pi(w)v|<=|v| and |z(w)|<=N0|v|.
Apply (4) to the tangent inhomogeneity q z. The profile ratios cancel
INSIDE the Duhamel integral. This proves

    ||Vhat(s,w)|| <= Kamb(S) P(s)/P(w).                  (6)

There is no exponential Gronwall constant, inverse Gaussian floor, or
new normal/angular component floor. A finite pole in edge distance and
sqrt(zeta) are constant along this lifted path and pass through (6).
Their parameter derivatives are handled on the actual unnormalized
fields, not by differentiating an assumed bounded reciprocal edge weight.

## 3. Explicit raw mixed-jet majorants

Use Euclidean vector/operator norms and the maximum of every raw mixed
derivative of total order j. For positive polynomial sequences set

    (X*Y)_j=sum_(i=0)^j binom(j,i) X_i Y_(j-i).          (7)

This bounds multivariate Leibniz because the sum of multibinomial
coefficients with |alpha|=i is binom(j,i). No factorial norm is silently
reused. Full-path inputs are normal jets n_j, moving normal jets n'_j,
K_j, p_j, and source jets F_j. The moving-normal input at order j uses
the normal at order j+1. All actual chart/pullback jet conversions must
be included in these inputs. The bound applies to each frozen input m;
the first family is m=+/-1,+/-2.

For r=1/N use

    D_j=(n*n)_j,  r_0=c_N^-1,
    r_j=c_N^-1 sum_(i=1)^j binom(j,i) D_i r_(j-i),
    h_j=(n*r)_j,  v_j=(n*K)_j+n'_j,
    A0_j=K_j+(h*v)_j,  H_j=(h*n)_j,
    lambda_j=4m^2 D_j,
    Ahat_j=A0_j+lambda_j+((lambda+p)*H)_j.               (8)

The last damping bound uses epsilon k^2<=4, fixed per label. If that
factor varies in continuous coordinates its FULL jets must instead be
included in lambda_j. The normal-square inverse, lift and adapted
connection are thus finite positive polynomial recipes, not separate
unverified inverse/connection norms. For Pi use Pi_0=1 (orthogonal),
Pi_j=H_j for j>=1, and G_j=(Pi*F)_j.

Assume the actual source satisfies the precise remaining inequality

    |partial^alpha f_m(s,x)|
      <= epsilon^B sqrt(zeta(x)) delta(x)^-J P(s,x) F_j(S),
      |alpha|=j, all s in [a,b], all accepted labels/scales. (9)

Use one common finite pole J covering the requested input jets. Fixed
unweighted coefficient polynomials in (8) must be uniform on this domain;
if they have poles, add those finite poles to the output allowance.
Here B=1-kappa is the proposed first-source order, conditional on (9).

Let U_j bound all order-j mixed velocity jets with the same weight and
epsilon^B. A sufficient recurrence is

    U_0=Kamb L G_0,
    E_j=Kamb L [G_j+sum_(i=1)^j binom(j,i) Ahat_i U_(j-i)],
    V_j=G_(j-1)+sum_(i=0)^(j-1) binom(j-1,i) Ahat_i U_(j-1-i),
    U_j=E_j+V_j, j>=1.                                 (10)

For parameter-only derivatives use differentiated (3), zero entrance
data and (6); lower mixed jets give E_j. For a derivative containing
at least one path derivative, differentiate t'=Ahat t-Pi f only j-1
times to obtain V_j. Induction proves (10) for every fixed finite order.
It does not require differentiated solutions to be tangent. Each output
uses finitely many coefficient jets and source jets through j; no new
starting band depending on derivative order is selected here.

Pressure and potential coefficient polynomials, normalized at order
epsilon^(B+1/2), are explicitly

    PiPress_j=|m|^-1 [((r*v)*U)_j+(h*F)_j],
    Pot_j=|m|^-1 (h*U)_j.                               (11)

These follow from (2), i(n cross t)/(km N), and k^-1<=sqrt(epsilon).
The UNPROJECTED F enters pressure. For f=n g, Pi f=0 and the zero-entry
velocity is zero, but pi=i g/(km) generally is not. A projected-source
norm alone is therefore insufficient for the pressure estimate.

## 4. Support, cutoff and what has not been closed

If f vanishes on a WHOLE slow/transverse parameter path, zero entrance
and uniqueness give zero t and pi there. Pointwise source support at a
single path location is insufficient. Endpoint zero values of f do not
force zero uncut exit velocity. Preserve the SAME actual cutoff on
potential AND pressure, then curl and smoothly extend; cutoff jet
products are the convolution (7). Full coefficient curl introduces its
retained single graph loss, so its order is B+1/2-kappa. Full velocity
remains W_B, conditional on (9) and the original support/cutoff lemmas.

The exact additive tail (1-psi)f+T_psi remains outside the forwarded
regular source. This note gives no new complete cutoff/viscous residual
norm and does not delete flat A0/A1 velocities from nonlinear products.
Complete base pressure, Pi1/Pi2, current moving pressure moments, and
fresh fast/slow updates remain governed by their existing full formulas.
The local pressure (2) is not a replacement for that mean pressure memory.

The next obligation is to supply F_j in (9) for the ACTUAL current
primary plus original-root signed wave after its prescribed fresh means,
retaining full viscosity, longitudinal curl cancellation, mean-wave phase
advection, cutoff velocities and whole-path supports. Neither the sharper
initial Haar rate nor (10) supplies that inequality. Also substitute the
actual n,K,p jets and the core/frame tangent bound into (8),(5). The
historical qualitative source/inverse foundations remain conditional;
this checkpoint does not independently certify their lower premises.
All later complex states, every physical scale, complete residual closure,
infinite realization and blow-up remain open.

## 5. Reproduction

`python test_audit.py -v` runs the new controls. `python audit.py
--verify-record evidence.json --output /path/to/fresh.json` reproduces
the frozen record and checks the declared historical bytes without
overwriting them. Exact moving-normal jets test (1)--(3) and the lift
decomposition; a separate finite polynomial solution bounds ALL mixed
jets over a full unit box by coefficient sums, not a mesh. Pure-normal
forcing, a non-tangent parameter derivative, conjugate input harmonics,
normal-rate adaptation, cutoff derivatives, finite input orders and
nonzero uncut exit are negative controls. These fixtures do not evaluate
the actual giant pulse or prove the analytical inputs of (4),(9).
