# Explicit background spectral and same-tilt cone certificate

Parent: `165b75afee625f8ec597e117e4986452907e7a3c`.
Status: `explicit-background-value-gap-certificate-relative-to-leading-and-weighted-jets`.
On 18 October 2026 the parent's new residual-ledger CI step passed; the
full workflow was still running. Historical sources and frozen evidence
remain unchanged. No unconditional pulse, residual or blow-up theorem follows.

The dependency trace found that `canonical_finite_cycle_join`, section 2,
already constructs the canonical initial route. Repeating that seed table
would not discharge a new premise. The next genuine estimate here replaces
the unspecified neighborhood radius and Lipschitz constant in
`weighted_background_edges`, section 5, by explicit rational expressions.
It supplies a sufficient value-domain certificate for the actual summed
background, GIVEN the reviewed leading and weighted-error estimates.
It does not independently prove the infinite core, regional accepted-root
inequalities, those weighted-error estimates, or actual numerical constants.

## 1. The input and a denominator-free growth discriminant

Use normalized band variables and write

    g=(R partial_R F, partial_R G), N=g/|g|, K=(-N_z,N_theta),
    D=|g|^2+2F g_theta, t=-g_theta.

Here g is the background SHEAR vector, not the wave's phase normal n.
The leading quantities have subscript 0. The existing leading-profile
certificate gives fixed positive f,t_*,d_* and finite Fbar,Gbar such that

    f<=F0<=Fbar, -g0_theta>=t_*, |g0|<=Gbar, D0>=d_*.

These constants need not be representable in floating point. In the notation
of `leading_frame_foundation`, one sufficient choice is

    f=f_*, Fbar=f^a, t_*=alpha*f, d_*=alpha*nu*f^2,
    Gbar=Fbar*(A_b+B_b),

since g0=F0(-a,b) and D0=F0^2*a*(v-2). Apply the stated positive
s^(-1-h) chart factors to f,Fbar,t_*,Gbar and its square to d_* BEFORE
using one uniform certificate over 1/4<=q/Q<=4. Cone ratios are invariant
under this common positive shear/swirl scaling. The existing positive
tilt u is FIXED, not chosen again by this estimate.

Let n0=T0/|T0| be the smooth extended leading stress direction. It remains
defined at the flat edges even though T0 itself vanishes there. Require
the reviewed leading orientation and cone margins

    H0=-N0.n0>=omega>0,
    sqrt((1+u^2)*c0^2)*|K0.n0| <= (1-gamma)*u*H0,
    0<omega,gamma<=1,
    c0^2=D0/(2F0*(-g0_theta)).                              (1)

The EXACT shear growth identities are

    lambda_shear^2=-2F*g_theta*D/|g|^2,
    c_shear^2=D/(-2F*g_theta).                              (2)

Thus positivity of F and -g_theta alone is insufficient: D must be
positive. Equations (2) concern background algebra, not redefinitions of
the construction's already frozen representative frame or carrier.

By the actual weighted-background and chart estimates, on one domain

    max(|F-F0|, |g-g0|, |r|) <= Cval*q^(2h),
    r=(That-T0)/|T0|,  That=|T0|*(n0+r),  Cval>=1.           (3)

The vector norms are Euclidean. The quotient in (3) is its smooth factored
extension at both stress edges, not division of numerical zero stresses.
The weighted-stress foundation supplies its specific smooth-edge argument.
A generic absolute stress bound would NOT supply this input. This is a
value/first-profile-jet input with finite Cval, not all high jets made small.

## 2. Explicit perturbation radius, with all denominators retained

Define the following finite positive constants from the input bounds:

    LD=4Gbar+2Fbar+3, LP=Fbar+Gbar+1,
    Dbar=Gbar^2+2Fbar*Gbar, cbar2=Dbar/(2f*t_*),
    Lc=2LD/(f*t_*)+2Dbar*LP/(f^2*t_*^2),
    LH=1+2/t_*,
    Lphi=3u^2*LH+(1+u^2)*(4Lc+3cbar2*LH),
    mphi=u^2*omega^2*gamma*(2-gamma),
    r_* = min(1/2,f/2,t_*/2,d_*/(2LD),
                       omega/(2LH),mphi/(2Lphi)).             (4)

If the maximum error delta in (3) is <=r_*, then F>=f/2,
t>=t_*/2 and |g|<=Gbar+1. Expanding D exactly gives

    Delta D=2g0.Delta g+|Delta g|^2
                 +2F0 Delta g_theta+2Delta F g0_theta
                 +2Delta F Delta g_theta.

For delta<=1 this has absolute value <=LD*delta, so D>=d_*/2.
Similarly |Ft-F0*t0|<=LP*delta. Since Ft>=f*t_*/4 and
F0*t0>=f*t_*, the exact quotient difference

    c^2-c0^2=Delta D/(2Ft)-D0*Delta(Ft)/(2Ft*F0*t0)

has absolute value <=Lc*delta. The leading |g0|>=t_* implies
|N-N0|<=2delta/t_* (the usual normalization difference inequality);
K has the same bound. No division by a possibly disappearing angular
carrier, output wave frequency, target stress, or pressure is used here.

Write H=-N.(n0+r), J=K.(n0+r). Then

    |H-H0|, |J-J0| <= LH*delta,
    H>=omega/2, |H|,|J|<=2.

The polynomial cone functional avoids differentiating its square root:

    Phi=u^2 H^2-(1+u^2)c^2 J^2,
    Phi0>=mphi,
    |Phi-Phi0|<=Lphi*delta.                                (5)

Indeed the H and J square differences are at most 3LH*delta;
|Delta(c^2 J^2)|<=4Lc*delta+3cbar2*LH*delta. Consequently
Phi>=mphi/2. This proves (4)--(5) on the WHOLE closed domain from the
displayed input inequalities, not by scanning parameters.

## 3. Explicit actual-background floors and the SAME tilt

Equations (2)--(5) imply

    F>=f/2, -g_theta>=t_*/2, D>=d_*/2,
    -N_theta>=t_*/[2(Gbar+1)],
    lambda_shear^2>=f*t_*d_*/[4(Gbar+1)^2],
    c_shear^2>=d_*/[4(Fbar+1)(Gbar+1)].                     (6)

For the actual normalized stress direction nB=(n0+r)/|n0+r|,
the factor norm is between 1/2 and 3/2. Thus -N.nB>=omega/4,
and the stress magnitude is at least |T0|/2 in the active INTERIOR.
There is no nonzero stress-magnitude floor claimed at its zero edges.
The smooth direction extends to those closed edges.

For the same u, (5) and H^2<=4 give

    ((1+u^2)c^2 J^2)/(u^2 H^2)
        <=1-omega^2*gamma*(2-gamma)/8.

Using sqrt(1-x)<=1-x/2 therefore gives the explicit cone margin

    gamma_new=omega^2*gamma*(2-gamma)/16>0.                   (7)

Normalizing n0+r cancels from the ratio. This is a possibly weaker
SUFFICIENT certified margin, not an assertion that the actual optimal
margin equals (7). In particular positive spectral growth alone does not
imply a stress cone for the prescribed u.

Once these value floors hold, every fixed jet of the smooth quotients,
positive square roots and normalization maps is bounded by the chain rule
using the already supplied fixed background jets. Higher derivatives change
constants, NOT r_* or the starting domain. This is real-domain smoothness;
it does not assert complex-parameter analyticity of a flat modulation or
of the heat endpoint. The ACTUAL moving pulse frame is still the prescribed
one with its nprime and Bprime terms; its inverse estimates are not replaced
by the constant-reference diagnostic below. Equations (6)--(7) only supply
the value-gap interface used by the reviewed frame perturbation argument.

## 4. The precise remaining value gate and historical-domain limitation

The new explicit sufficient gate is

    Cval*q^(2h)<=r_*.                                       (8)

It holds for all q<=qcap=min(qexisting,(r_*/Cval)^(1/(2h))).
For a full dyadic chart q<=4Q, Q=2^-ell, a sufficient threshold is

    ell >= 2+log2(Cval/r_*)/(2h),                            (9)

rounded upward and combined with EVERY previous geometry, normal, cone,
root and corrected-freezing requirement. The extra 2 is the actual 4Q
chart cost; it must not be omitted. This does not retune h,N,lambda,
the width C_setup^-20, leading roots, primary roots, carriers or Borel
coefficients. It restricts the already defined background to a later tail.

The original qexisting condition obtained from an unspecified continuity
constant need NOT imply our more conservative sufficient (8). No numerical
Cval or actual huge leading constants are evaluated here, and failure of
this sufficient test would not prove failure of the actual cone. For a
historically fixed first band, (8) is the precise unverified inequality;
it is not certified by manufactured constants. A later distinct realized
field/force must use its OWN controlled schedules and be rejoined as such,
not identified with the historical field by restricting a parameter label.

Full Pi1/Pi2, the source-normal and moving-normal wave pressure, actual
A0/A1 velocities and the complete residual ledger remain untouched. A
spectral/cone value certificate is not a pressure identity or full PDE
residual estimate. The previous distinction between temporal Gaussian
tails and regular spatial-localizer derivatives remains in force.

## 5. Validation, limits and next obligation

Exact rational perturbation controls independently evaluate every term of
the quotient and polynomial-cone bounds. Reference-frame controls vary
rational F,c0 and Pythagorean s,sqrt(1+s^2), checking the complete 2x2
diagonalization. These are manufactured constant reference frames, NOT
the selected huge profile or its moving-frame evolution.
Near-neutral examples keep F>0 and g_theta<0, but an arbitrarily small
absolute error changes the growth sign; another exact example has strict
growth and fails the same-tilt cone. Relative-stress omission controls
detect orientation reversal by a tiny absolute stress change at a tiny
stress. These are obstructions to invalid inferences, not counterexamples
to the actual construction. Exact selector controls use a permitted
diagnostic h=1/200; they do not select the construction's h or first band.

Historical byte bindings are verified recursively, separately from fresh
test/evidence replays of the relevant foundations. Tests are not continuum
proofs. Analytical implications (2)--(9) use the displayed leading and
weighted inputs and do not discharge them by declaring a flag true.

Next independently audit the actual retained regional small-root/value
inequalities and normalized leading/weighted-error constants feeding (8),
or use a rigorously controlled later-tail existential selection with all
those inputs established. Do not repeat the already completed canonical
seed route or promote this one value-gap certificate to an unconditional
initial-field, full residual, single-field/force or blow-up theorem.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json`.
Never overwrite historical evidence.
