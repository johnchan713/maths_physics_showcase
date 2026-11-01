# Q81: exact exterior cancellation in the principal covariance comparison

Parent: `a8b1efc0a7a3582209bd44621dc5b471fee704c3` (Q80).
Q2, Q25 and Q26 already join the leading cone, canonical columns and
rule-defined later-tail domain, relative to their stated primitive estimates.
Q79's enlarged upper propagator bound does not extend that positive-column
domain. We retain it. The new result here is support localization of Q61's
COMPARISON ERROR, using Q80's exact positive-order velocity support.

This concerns the principal tensor after angular averaging. It does not
construct a new physical wave, replace a root, or remove exceptional pressure
memory. All estimates below retain the original same-field background,
closed-label and canonical inverse hypotheses. Historical records stay frozen.

## 1. The common target and fixed slow parameter are indispensable

Let V=Xvel be the positive-order velocity support endpoint. Q80 gives
F=F0, G=G0 for X>=V, with every available jet equal at the boundary as
well. The superscript zero denotes the leading-background comparison of
Q61, keeping the SAME actual carriers, epsilon, k, L, reference rate P,
cutoffs, frames, slow partition, Haar maps and target Tbeta. In Q59 this
unstripped target is the prescribed scaled T0; Q61 explicitly holds it
fixed. No argument below replaces any full pressure-dependent target by a
leading-pressure approximation.

Fix a slow point x with X(x)>V. Auxiliary path time v varies with x fixed:
X(x) is constant on the whole pulse path. The exact phase formulas give
n=n0, m=m0 and K=K0 for every v. The complete ambient generators and
canonical entrance vectors are consequently identical, including moving-
normal pressure and axial normal motion. Uniqueness of the finite-interval
ODE gives t=t0 on the whole path. Multiplying by the unchanged cutoffs,
and integrating over the unchanged auxiliary Haar fiber, gives

    h_sigma=h_sigma^0, Hbar=Hbar^0.                     (1)

This is equality of the entire column integral, not a conclusion from a
pointwise sampled direction. The Haar fiber holds the slow variables fixed;
a radial average mixing different X would not justify (1).

Where the actual column is invertible, (1) makes the reference column
invertible too. With the SAME target and s0,

    z=Hbar^-1 Tbeta=z0,
    P_ab=eta_beta^2 z_sigma h_sigma,a h_sigma,b/s0
        =P_ab^0.                                      (2)

No positive square root of z0 is needed. On any common smooth comparison
domain, identities (1)--(2) hold on the exterior open set, and their jets
extend as zero difference to X=V by smoothness. Thus

    supp D^I(P-P0) is contained in {X<=V}               (3)

for every fixed available derivative. A whole chart with r^2/8>=V has
identically zero difference, regardless of its fixed-band epsilon size.
For a crossing chart (3) is a support theorem, not permission to insert a
new discontinuous cutoff into the wave.

## 2. The low-order axial error has the same input support

Q58's diagonal weights do not change. Its shear weights are

    w3=2R G_U,
    w4=(2C-R^2)F_U-4CD F, C=s mu, D=R_U/R.             (4)

The cutoff-shape moment C is common to both weights, with its full slow
motion retained. It is not the physical pressure itself. Q80 and (4) imply
Delta w3=Delta w4=0 with all jets for X>=V. Hence the exact replacement
error in the principal weighted integrand is

    wP-w0 P0=(w-w0)P+w0(P-P0).                         (5)

Both terms in (5), and their slow derivatives, have the support (3).
This includes the derivative of C; equality does not require C_Z=0.
The same statement holds after the purely auxiliary linear projections
used in the moment, since those do not mix X. The Q58/Q61 direct O4
cost can therefore use correction and product seminorms on the intersection
with X<=V, rather than the entire exterior hull. The HIGH-order leading
comparison integral still uses the whole original hull.

For the positive original radial coordinate U=R^d, d>0, at fixed s the
intersection of [Ulo,Uhi] with that support has length

    Leff=max(0,min(Uhi,(2sV)^(d/2))-Ulo).               (6)

A direct integral of a bounded error is at most Leff times its supremum.
Any nonconstant kernel still incurs its own existing bound. Formula (6)
is a sharpening of the input integration-length factor only; it supplies
no extra epsilon power. When differentiating with respect to slow variables,
use the globally smooth zero-extended integrand before taking the norm.
All its jets vanish at the moving endpoint X=V, so no endpoint distribution
or omitted boundary flux is created. We do not differentiate the piecewise
majorant Leff as though it were the actual field.

## 3. Two restrictions prevent an unsupported global improvement

First, Q61's successful comparison branch supplies a common reference
normal, propagator and column inverse on the full hull. Equations (3)--(6)
apply there. If its fixed-band smallness gates fail, its fallback estimates
the TRUE tensor directly. Exterior equality alone does not construct a
smooth reference inverse on the rest of that hull, nor license replacing
that fallback by a compact difference integral. No new band threshold or
patched inverse is introduced. On wholly exterior charts reference equality
is valid wherever the original actual inverse is defined.

Second, a shifted radial inverse or a backward pressure primitive need not
preserve radial support. For example a smooth bump f supported in (0,V)
with integral one has integral_0^X f=1 for X>V. Localization of the INPUT
error in (5) does not show that a later A0/A1 velocity or pressure output
vanishes there. Its complete kernel and moment conditions must be retained.
The reserved pressure still has Pi1=K1(eta)-mu^2/(4X) and Pi2=K2(eta).
Replacing the common target in (2) by a different one would instead give
z-z0=Hbar^-1(Tbeta-Tbeta^0), generally nonzero even when the columns agree.

## 4. Result and next obligation

The comparison error in Q61 and the shear-weight error in Q58 are confined
to the positive-order velocity support, on their legitimate comparison
domain. This removes unnecessary exterior inputs from their low-order O4
remainder estimate; it does not remove the leading high-order exterior
integral, the failed-comparison branch, or either exceptional pressure tail.
No unconditional pulse bound, improved complete residual, or blow-up follows.

Next evaluate the remaining high-order LEADING principal moment and its
radial pressure counterpart using the established profile, path, target and
column bounds on the certified domain. Keep the full shifted inverse and
complete-wave covariance/axial-pressure mean terms. The inherited primitive
contracts and all viscosity, cutoff and interaction residuals remain.

Reproduce with `python test_audit.py -v` and `python audit.py
--verify-record evidence.json --output FRESH.json`. Exact manufactured
controls test the common-target cancellation, its failure when targets
change, moving-boundary flatness at the tested orders, and nonlocal support
failure. Finite tests do not establish the continuum premises.
