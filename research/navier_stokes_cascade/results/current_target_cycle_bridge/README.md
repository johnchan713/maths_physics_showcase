# Exact current-target cancellation through signed and mean updates

Parent: `f370f9155a7286d613ad4e00ba5e22a3ddf43b0d`.
Status: `conditional-current-target-cancellation-cycle-bridge`.
On 18 October 2026 the parent's checkpoint CI step passed; the complete
workflow was still running. Historical sources and frozen evidence are intact.

The parent isolated the missing inequality [Ae+Ie he]_m<=C_jm epsilon^alpha_j.
Here it is derived by EXACT CURRENT-STATE DIFFERENCES through the prescribed
generic cycle, relative to the reviewed wave/mean map estimates and an accepted
initial canonical input. Raw field membership alone is still insufficient.
The parent's smooth swirl control does not satisfy the cancellation recurrence
below, and remains valid. No actual giant field or first band is evaluated.
This closes a conditional target-estimate loop, not all analytical premises,
the complete residual/forcing theorem, or an unconditional blow-up proof.

## 1. A linear target map applied to nonlinear CURRENT equations

Retain epsilon=Q^h, kappa=10^-5, all fixed jets, full zeta weights with finite
inverse-edge powers, the actual width C_setup^-20, complete pressure memory,
the same fixed canonical roots/frame, and the corrected initial accepted-tail
selector (including ell^6>=8*Cfreeze/gamma0 and the frozen cone margin in Lreq).
All constants retain physical radius, fixed lambda, jet and finite-stage costs.
There is no jet- or stage-dependent first-band choice and no new positive root.

Let e=2 for the angular row and e=1 for the axial row, and define

    Ie f = R^-e int_left^R r^e f dr,
    Qe f = f-sigma_e int R^e f,        int R^e sigma_e=1,
    Xe(Ee) = -Ie Qe bar(Ee).                                  (1)

The parent proves that these are the actual compact signed targets because
both velocity moments vanish. Its pressure/flux identities give the current
integrals epsilon*partial_Z Jtheta and
epsilon*partial_Z(Jz+c_rho P), respectively. Average products AFTER multiplication.
On compact full-zeta inputs the oriented primitive estimate bounds (1) on
every fixed jet. Qe is a bounded finite-rank operator because its bump has
the same fixed interior support hull. No trace at a shell edge is divided by.

This map is linear in Ee, not in the fields. Therefore for ANY actual update

    Xnew-Xold=-Ie Qe bar(Delta Ee).                            (2)

The cancellation is in (2). Bounding new raw terms separately loses it, as
the parent's countercontrol demonstrates. Conversely assuming a desired bound
for Enew instead of proving its difference would be circular.

Moving compensation is important. With Dz=epsilon*partial_Z,

    Qe Dz H = Dz(Qe H) + epsilon*(partial_Z sigma_e)*int R^e H. (3)

The last term is generally nonzero; it has one extra epsilon, not an excuse
to erase it. Likewise all axial pressure derivatives act on the WHOLE c_rho P.
We use Qe Dz, never commute them silently.

## 2. Exact covariance response, then the actual signed cancellation

During a wave update hold mean velocities fixed. Let D_ab=bar(Delta W_ab)
be the change of the COMPLETE angular covariance, not its principal two
entries. Its exact radial source and pressure changes are

    dg=-(partial_R+1/R)Drr-Dz Drz+Dtheta,theta/R,
    dP=int dg, dp=I0(dg-rho*dP),
    dJtheta=int R^2 Dtheta,z,
    dJz=int R Dzz-(1/2)int R^2 dg.                            (4)

The input of dp has zero total. All fixed full-zeta jets are therefore
controlled with the original edge weights. After Haar averaging no graph
loss remains: D in M_a implies dg,dp,dP,dJ in order a. Unaveraged A0 pressure
remains in the equations; (4) does not remove it from later products.
The exact target changes from the conservative equations are

    dXtheta=-Drtheta-I2 Q2 Dz Dtheta,z,
    dXz=-Drz-I1 Q1 Dz(Dzz+dp).                                (5)

Both pressure reconstruction and moment compensation use NEW current data.
All six covariance entries enter (4),(5). In particular dp is not supplied
by the two selected principal fluxes. Equations (4),(5) retain the moving
bump terms in (3) and in partial_Z(c_rho*dP).

For the actual signed complete curl V=L+R about the fixed canonical principal
wave W0c, write U-W0c=E. The complete tensor identity, with angular averaging,
is exactly

    Delta W = Lambda + Bil(E,L)+Bil(U,R)+Cov(V),
    Lambda=Bil(W0c,L),    Bil(a,b)=<a tensor b+b tensor a>theta. (6)

The fixed pair inverse gives bar(Lambda_rtheta,Lambda_rz)=Xold from the
CURRENT target, including the angular half and the factor two. The other
entries of Lambda need not vanish. Define Derror=bar(Delta W-Lambda).
Substitute (6) into (5), rather than merely estimating (1):

    Xtheta,new=-Derror_rtheta-I2 Q2 Dz Dtheta,z,
    Xz,new=-Derror_rz-I1 Q1 Dz(Dzz+dp).                        (7)

This proves the parent cancellation inequality at this substep. If Xold
has order C, U in W_1/2, E in W_rho, L in W_(C-1/2), and R in W_(C-kappa),
the full covariance change is order C, while its error is order

    nu=min(rho+C-1/2, C+1/2-kappa, 2C-1).
    Xnew in M_min(nu,C+1).                                   (8)

Products of complete sqrt(zeta) waves have full zeta. Leibniz and the exact
current inverse give (8) for EVERY fixed jet, with finite polynomial S and
inverse-edge costs. Generic complex old/new crosses and all signed squares
are retained; the real initial cosine/sine shortcut is not used here.

A sourced complete wave of order B changes the full covariance at order
C=B+1/2 (since B>=1/2), including its self-square. Equations (4),(5) give
dX in M_C. Thus it preserves a pre-existing order-C target and supplies the
CURRENT target to (6),(7). It is not a historical inverse response.

## 3. Fast means: the linear Haar response vanishes EXACTLY

Keep W fixed, and let the ACTUAL fast increment d=(db,dv,dgamm) have zero
Haar mean, tangential order H and radial order H+1. Its curl/A1 correction
has zero Haar mean too; replacing it by a desired principal axial velocity
is neither needed nor allowed. The existing mean map also preserves both
weighted velocity moments and divergence.

All base-linear, time and viscous Haar changes vanish exactly. The remaining
averaged radial source and tangential fluxes are

    dg_fast= -D1^0 bar(2beta*db+db^2)
             -Dz bar(beta*dgamm+gamma*db+db*dgamm)
             +bar(2v*dv+dv^2)/R,
    dFtheta=bar(beta*dv+db*v+db*dv),
    dFz=bar(beta*dgamm+db*gamma+db*dgamm),
    dHtheta=bar(gamma*dv+dgamm*v+dgamm*dv),
    dHz=bar(2gamma*dgamm+dgamm^2)+dp_fast,                     (9)

where dp_fast=I0(dg_fast-rho*int dg_fast). Apply (2) to the radial divergences
of dFe plus Dz dHe. With old tangential means M_mu, radial M_(mu+1),

    dP,dJ in S_min(mu+H,2H),
    dX_fast in M_min(mu+H+1,2H+1).                            (10)

The extra one in the target response follows from the radial mean factor
in dFe and Dz in dHe, not from pretending that the old means are zero.
Old/new fast correlations in (9) generally survive Haar averaging. In
particular the defects CHANGE; the slow five-map must read the POSTFAST
targets. With H=C-kappa and mu=1-kappa, stale defects can leave order
C+1-2kappa, insufficient for the desired postslow C+1 scalar gain.

## 4. Slow means: use increments, not a raw-state smallness assertion

The reserved five-map at current defects of order C has tangential increment
d in M_C and radial increment db in M_(C+1). It is auxiliary-independent,
compact on the fixed I4 POTENTIAL HULL, and preserves the two moments. Actual
b is O(epsilon); the base V,G and their fixed coefficient jets have their
reviewed bounds. Include every slow-time, vector-viscous and pressure term.

The exact current source change is

    dg = -t*db-D1(2b*db+2beta*db+db^2)
         -Dz(b*dgamm+G*db+beta*dgamm+gamma*db+db*dgamm)
         +(2V*dv+2v*dv+dv^2)/R
         +epsilon*(Delta0-R^-2)db.                           (11)

The tangential equations change by radial divergence of the complete new
fluxes, Dz of their complete axial fluxes (including dp), time, and viscosity.
Here dg,dp have leading order C because of 2V*dv/R. That pressure is not
order C+1, but its axial momentum contribution Dz dp IS order C+1.
All resulting target changes obey

    dX_slow in M_min(C+1,mu+C+1,2C+1)=M_(C+1).                (12)

For example b*dv, V*db, epsilon*partial_R dv, epsilon*partial_T dv,
epsilon*partial_Z(V*dgamm), and Dz dp each have order C+1. Old/new slow
products have order mu+C+1, and new squares have order 2C+1. Haar averaging
turns graph derivatives into ordinary radial derivatives. No claim of a
uniform-in-lambda inverse is made; its fixed 1/lambda cost remains.
Nonlinear postslow defects are not zero. Their earlier C+1 estimate still
requires the fresh current linear rows, not a historical five-map target.

## 5. One accepted seed, then every finite target cycle

The input bundle is the canonical finite-cycle state, not JUST a target:
regular source W_(C-1/2), full tangential mean rows M_C and defects S_C,
current target M_C, old wave U in W_1/2 with U-W0c in W_(1-kappa), total
tangential/radial means M_mu/M_(mu+1), exact divergence and both moments.
The reviewed sourced/signed maps change full unaveraged rows at order
H=C-kappa because graph derivatives cost kappa, while their Haar response
is precisely (5)--(8). Thus the CURRENT fast map has order H; it is not
assigned a high order merely because the target is compact.

There is also no need to assume the desired POSTFAST full mean bound.
From (1), bar(Ee)=-De^0 Xe+sigma_e Me. The signed target estimate (8) and
current defects S_C give bar(Ee) in M_min(nu,C+1). The fast inverse cancels
the complementary E^circ exactly up to the retained actual A1 term. Its
remaining full mean orders are H+1, H+1-kappa, H+1-2kappa (viscosity),
mu+H+1-kappa, H+1 (pressure), and 2H+1-kappa. The smallest linear row is
C+1-3kappa, which exceeds Cnext by 1/2-2kappa>0. The compact-mean interior
A0/A1 estimates permit every fixed finite order while retaining the actual
velocities and remainders. Slow full row changes have orders C+1,
mu+C+1-kappa,2C+1, also above Cnext. Hence the corresponding next FULL
mean-row bound follows; it is not used to assert its own Haar target gain.
The complete nonzero-angular source bound still uses the separately
reviewed residual operators and is not proved by these mean identities alone.

Let rho=mu=1-kappa and C0=2-2kappa. For C>=C0, the sufficient target orders
through the ACTUAL sequence are

| Current operation/group | Target order |
|---|---:|
| Sourced full covariance change | C |
| Signed old-error cross | C+1/2-kappa |
| Signed full-curl cross | C+1/2-kappa |
| Signed square | 2C-1 |
| Signed axial covariance and pressure response | C+1 |
| Fast old/new means | C+2-2kappa |
| Fast square | 2C+1-2kappa |
| Slow linear / old means / square | C+1 / C+2-kappa / 2C+1 |

Each order minus Cnext=C+1/2-kappa is affine in C, with nonnegative slope
and nonnegative value at C0. The first two signed error groups tie the
bound; no arbitrary larger gain is asserted. This proves the whole-half-line
implication X in M_C -> Xnext in M_(C+1/2-kappa) for the prescribed current
maps, relative to their retained operator and source bounds. The target
estimate is now a CONSEQUENCE of (2),(7),(9),(11), not an extra per-stage
assumption. Constants may depend on the finite stage and fixed jet.

For seed routing, begin with the accepted canonical initial-mean estimate
alpha0=2-3kappa from canonical_initial_mean_join. Its complete real primary
has the reviewed real curl-square identity. The real signed error gives
3-6kappa; the real fast/slow target responses have orders 4-5kappa and
3-3kappa or better. The first sourced solve has B=1-kappa and covariance
order C1=3/2-kappa. Generic signed errors then all tie at 2-2kappa=C0;
the following fast/slow target responses are at least 7/2-3kappa and
5/2-kappa, respectively. Thus (7)--(12) route that SAME accepted canonical
input to the required target seed. No historical root, field, or target
value is transplanted. This does not independently prove the retained
initial canonical wave and mean estimates or their background premises.

All fixed-jet functional hypotheses come from the parent and source/edge
audit. The accepted tail is unchanged across this induction, including all
physical scale costs. This is NOT a bound uniform in stage, an evaluation
of actual constants, or permission to use an uncontrolled infinite diagonal.

## 6. Exact tail routing and remaining limits

Linear wave terms, including the pressure/curl cutoff commutator and
(1-psi)f, have nonzero angular frequency and hence zero angular mean in
the cylindrical component equations. They do not add a missing term to
(4)--(12). Quadratic transport uses COMPLETE cut velocities, so all their
angular means are already in the full covariance. Mean A0/A1 remainders
are retained: linear zero-Haar pieces vanish only in this averaged target
map, while products of actual velocities stay in (9),(11). The flat BASE
residual belongs to the separately subtracted full base balance, not these
increment targets. Pi1/Pi2 cancel only against that SAME complete base.

No arbitrary physically flat tail is normalized or forwarded through a
sourced inverse. This angular-mean routing does not by itself validate the
ENTIRE nonzero-angular residual ledger or its pressure-based reclassification.
That remains the next obligation: join the now explicit target recurrence
to all sourced/mean residual increments, full viscosity, cutoff tails and
interactions with actual support and physical-jet estimates on every accepted
scale. Retained core/background/initial-map hypotheses, historical first-band
and old-force acceptance, and unconditional conclusions are not recertified.

## 7. Checks

`target_cycle.py` uses exact manufactured current conservative equations,
not protocol flags, to compare independently recomputed targets after a full
six-entry covariance change, an actual zero-Haar curl increment and a slow
curl increment. Its fixed stress is engineered to give a specified current
target, solely as a finite tensor fixture; it is not an actual canonical
wave or a certified positive covariance realization. It detects nonzero axial
pressure, moving compensation and old/new mean correlations. The all-jet,
all-scale continuum argument is (1)--(12), conditional on the named analytical
inputs; polynomial endpoint tests and affine budgets are not continuum proofs.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json`.
Fresh outputs only; historical frozen evidence is never overwritten.
