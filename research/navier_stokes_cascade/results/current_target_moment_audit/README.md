# Current signed targets: exact moment map and a quantitative cancellation gate

Parent: `08caced6f10105d98faf59af4e251af26c3a2791`.
Status: `current-target-full-weight-identities-with-higher-order-cancellation-gate`.
The parent's new checkpoint CI step passed on 18 October 2026; its complete
workflow was still running. Historical sources and frozen records are unchanged.

This supplies a direct audit of the CURRENT targets left open by the
[source/edge checkpoint](../source_envelope_edge_audit/README.md). It proves
full-zeta functional membership, exact compactness and both integrated mean
identities for each finite current state satisfying the stated field hypotheses.
It also isolates an additional ORDER requirement: field/support/moment bounds
alone do not imply arbitrarily improving signed-target orders. The older
finite-cycle argument remains conditional on those cancellations; its protocol
and exponent tables are not used as proof that a current target satisfies them.
No actual enormous field is sampled, no historical first band is certified,
and no unconditional pulse, complete residual or blow-up theorem is claimed.

## 1. Current fields, not historical values

Use the same accepted chart, epsilon=Q^h, fixed width C_setup^-20 and complete
base (b,V,G,p_B). The corrected initial tail selector, including the separate
ell^6>=8*Cfreeze/gamma0 condition and the frozen margin in Lreq, is retained.
There is no new domain, root, carrier, cutoff or realization here.

Write the ACTUAL mean increment (beta,v,gamma), with
D_1 beta+epsilon*partial_Z gamma=0, and the complete angular covariance
W_ab=<w_a w_b>_theta. A bar below is auxiliary Haar average, not angular
average. All six covariance entries are retained. Assume compact smooth
radial support in a fixed positive shell, the reviewed slow supports, all
fixed coefficient jets with full zeta weights for means, W and the retained
base stress Sigma (with finite inverse-edge powers), the reviewed base
coefficient bounds, and the two
constraints int R^2 bar(v)=int R bar(gamma)=0. Mean maps may fill the torus.
The complete A0/A1 velocities are part of beta,v,gamma, not discarded errors.
Flat base and designated cutoff residual ledgers remain separate.

Put B=bar(beta), v0=bar(v), g0=bar(gamma), Qab=bar(Wab), and

    F = 2bB + bar(beta^2) + Qrr,
    H = b*g0 + G*B + bar(beta*gamma) + Qrz,
    U = 2V*v0 + bar(v^2) + Qtheta,theta,
    K = G*v0 + V*g0 + bar(gamma*v) + Qtheta,z,
    Zf = 2G*g0 + bar(gamma^2) + Qzz,
    tau = epsilon*partial_T + epsilon^3*partial_ZZ.               (1)

Products are averaged AFTER multiplication: bar(beta*gamma) is not B*g0.
Haar invariance removes every graph derivative and the fast time derivative,
including those of actual curl remainders. The current radial source is

    gbar = tau*B - (partial_R+1/R)F - epsilon*partial_Z H
                   + U/R + epsilon*(partial_RR+R^-1 partial_R-R^-2)B.
                                                                  (2)

The signs follow the conservative equations in compact_mean_corrections,
not a swirl-only model. The base may have nonzero radial flow and axial
velocity. Equation (2) is an identity, not an exponent estimate.

## 2. A closed current moment map, including radial time and axial flux

Let rho be the same moving unit-integral pressure bump, and
c_rho=(1/2)int R^2 rho. Radial integration by parts gives exactly

    P = tau int B + int (U-F)/R - epsilon*partial_Z int H,
    Qg = int R^2 gbar
       = tau int R^2 B + int R(F+U) - epsilon*partial_Z int R^2 H,
    Jtheta = int R^2 K,
    Jz = int R Zf - Qg/2.                                       (3)

Both radial vector-viscosity integrals vanish. All endpoint terms vanish
by actual compact smooth extension, not by Gaussian smallness alone. No
radial graph derivative appears in (3); thus moments cost slow derivatives,
not a hidden epsilon^-kappa graph loss. For every fixed jet they are bounded
by finite sums of current field jet products, with finite shell constants.
In particular P is NOT just the initial covariance anisotropy: it also
contains tau int B, both mean squares, radial base transport and the FULL H.

The Haar pressure is the ordinary compact primitive

    pbar(R)=int_left^R (gbar-rho*P),
    int R pbar = -Qg/2+c_rho*P.                                 (4)

Its input has zero total by definition, so the two-edge estimate of the
parent applies. This is not an inference from smoothness alone. The actual
unaveraged pressure remains T0(gr-rho P), with its nonzero A0 term; averaging
does not license its removal from later nonlinear products.

Define the radial fluxes and axial fluxes

    Atheta = bar(b*v+beta*V+beta*v+Wrtheta-Sigma_theta)
                         - epsilon*(partial_R v0-v0/R),
    Az = bar(b*gamma+beta*G+beta*gamma+Wrz-Sigma_z)
                         - epsilon*partial_R g0,
    Htheta = K,  Hz = Zf+pbar,
    ftheta = epsilon*partial_Z Htheta-tau*v0,
    fz = epsilon*partial_Z Hz-tau*g0.                            (5)

Here b,V,G,Sigma are auxiliary-independent. Direct differentiation, including
the cylindrical viscous connections, gives

    bar(Etheta)=D2^0 Atheta+ftheta,
    bar(Ez)=D1^0 Az+fz,           De^0=partial_R+e/R.             (6)

Integrating with R^2 and R respectively, and using the two velocity moments,
proves FOR THIS CURRENT STATE

    int R^2 bar(Etheta)=epsilon*partial_Z Jtheta,
    int R bar(Ez)=epsilon*partial_Z(Jz+c_rho*P).                 (7)

Without the velocity constraints retain -tau Mtheta and -tau Mz in (7).
The derivative acts on the WHOLE moving product: c_rho,Z*P need not vanish.
Wrr, Wtheta,theta, Wrz, Wtheta,z, Wzz enter (1)--(5); the sixth entry Wrtheta
enters Atheta. None may be replaced by just the two matched principal entries.
Pi1/Pi2 cancel only on subtraction of the SAME complete base balance; (3)--(7)
do not assert that p_B equals leading p0 or that its axial derivative is zero.

## 3. An explicit signed-target formula and all-fixed-jet inequality

For auxiliary-independent data let Ie f=R^-e int_left^R r^e f(r)dr.
Let the same interior sigma_theta,sigma_z have integrals int R^2 sigma_theta
=int R sigma_z=1. Set

    Mtheta=epsilon*partial_Z Jtheta,
    Mz=epsilon*partial_Z(Jz+c_rho*P),
    htheta=ftheta-sigma_theta*Mtheta,
    hz=fz-sigma_z*Mz.

By (6),(7), int R^e he=0. Therefore the ACTUAL signed targets have the
equivalent exact representations

    Xi_theta=-Atheta-I2 htheta,
    Xi_z=-Az-I1 hz.                                             (8)

Indeed Ie De^0 Ae=Ae for compact Ae, and the original definition is
-Te[bar(Ee)-sigma_e Me]; here Te=Ie because its total moment is zero.
This uses no historical target, protocol stamp or principal-covariance
substitution. The radial primitives cannot be replaced by a one-edge
primitive without the moment condition.

Let [.]_(m) denote a finite full-zeta coefficient seminorm through order m,
allowing the finite inverse-edge powers required by those derivatives. The
parent's oriented primitive estimate and Leibniz give

    [Xi_e]_m <= C_m( [Ae]_m + [he]_m )
      <= C_m( [Ae]_m + epsilon*[He]_(m+1)
                    + epsilon*[ue]_(m+2)
                    + [sigma_e Me]_m ),                         (9)

where ue=v0 or g0 and epsilon<=1. Constants retain the fixed shell, bump,
width and all derivative/edge costs; no uniform-in-m or stage bound is claimed.
All derivatives of c_rho*P and of moving bumps are inside these seminorms.
Products of full-zeta fields stay full-zeta because zeta^2<=zeta, and complete
wave products have zeta because each wave has sqrt(zeta). Reserved bumps
are on the fixed INTERIOR potential hull where zeta has a positive floor.
The compact mean primitives preserve this class at both edges. Equations
(1)--(9) thus prove full-zeta functional membership and smooth extension
of the current Xi at every finite expression with these field hypotheses.
They do not prove its desired epsilon^alpha order by themselves.

For the sharper order the precise remaining estimate is

    [Ae+Ie he]_m <= C_(j,m) epsilon^alpha_j,                     (10)

with the same full-zeta weights, at EVERY current signed substep. A sufficient
but stronger condition is that both terms separately have that order.
Testing (10) on old targets, assuming a small radial flux alone, or replacing
bar(beta*gamma) by a product of means does not establish it.

## 4. Why the higher-order gate is real, not a support issue

There is a smooth continuum control with all those supports and both velocity
moments, but with target of order only mu+1. This is NOT the proposed actual
canonical pulse and does not refute a further cancellation theorem.

Take a fixed nonzero C-infinity bump psi in the shell and f=R^-2 psi', so
int R^2 f=0. Set beta=gamma=w=b=G=Sigma=0 and v=epsilon^mu a(T)f(R),
independent of Z and the torus, with a smooth compact temporal bump. Choose
any fixed base swirl V=V(R). Full current
radial pressure, P and Jz are still reconstructed from 2Vv+v^2; their axial
derivatives vanish, and Jtheta=0. At an interior T0 choose a(T0)=0,
a'(T0)=1. Then (8) gives EXACTLY

    Xi_theta(T0)=epsilon^(mu+1) I2 f
                =epsilon^(mu+1) psi/R^2,       Xi_z=0.          (11)

All fixed full-zeta seminorms are finite by choosing psi on an interior
compact hull. No boundary or pressure-tail omission is used. Hence a bound
with alpha>mu+1 and a constant independent of epsilon would require
epsilon^(mu+1-alpha)<=C as epsilon tends to zero, which is false.
The actual preserved means have mu=1-kappa, so this necessary raw-state
ceiling is 2-kappa. This is NOT a sufficient bound for arbitrary raw states:
unmatched leading radial covariance fluxes can already be only order 1.
The canonical cycle's desired alpha_j=B_j+1/2 exceeds the control's ceiling
already at j=1 (alpha_1=5/2-3kappa). The finite-cycle residual cancellations
may remove (11), but must be shown for the actual recomputed state via (10).
This control neither satisfies nor is claimed to satisfy that extra residual
induction hypothesis. It corrects a possible inference, not frozen evidence.

## 5. Validation and next obligation

`target_algebra.py` independently computes the unaveraged graph equations
for rational Laurent-Fourier/polynomial current data, Haar averages them,
and checks (3)--(8), exact divergence, both moments, compact pressure, and
the moving-bump derivative. A second current-state perturbation checks
the full linear plus quadratic radial-source and moment increments; no
old mean products are silently removed. Polynomial controls have finite
endpoint regularity; the all-order arguments above use actual C-infinity
supports. The symbolic fixtures are not actual pulse fields or continuum
supremum certificates. Frozen source digests are checked transitively.

Next establish (10) through the canonical sourced/signed/fast/slow sequence
using the actual current residual decomposition, including its full pressure,
all six covariance entries, zero-Haar correlations, reserved support hull and
exact additive tail reclassification. Until then an exponent table alone
does not close the higher-order target loop. Uniform actual physical pulse
and complete residual/forcing obligations remain conditional. Run:

```sh
python test_audit.py -v
python audit.py --verify-record evidence.json --output fresh.json
```

Fresh outputs only; never overwrite historical evidence.
