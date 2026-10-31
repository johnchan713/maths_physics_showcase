# Q62: reconstruct the unstripped axial-test target from velocity and terminal memory

Q59–61 use the actual leading stress T₀ linearly, before any square root
or division by the edge weight. Its high seminorm need not be supplied
independently of the high leading velocity seminorm. Here is a bound
on the entire original stress shell, with explicit terminal data. It
is conditional on the same finalized leading fields, exact five global
identities and terminal heat join. The pressure gauge cannot be reset
at either edge. No bound on a normalized stress direction follows.

## Remove the field reciprocals before estimating

Use profile coordinates X,η, r=√(2X), E=E₀, U=U₀, H=rE,
d=1−η², L=1−2hη², A=1/2+h, D=1/2−h. Define
w=X−2DηM−d Mη and V=(1−h)I−Dη Iη−d Jη+2(h−D)ηJ.
With subscripts denoting differentiation, the two **unweighted** entries are

    r Tθ = (−E w+V/r)/L − E + 2X E_X,
    r Tz = [−w U+D(M−ηMη)+4hηS−d Sη
                  +X(4AηΠ−d Πη)]/L + 2X U_X.                  (1)

Indeed T₀=(E/r)(p₁−a,p₂+b), a=1−2X E_X/E, b=2X U_X/E.
Substitute the full forward formulas (6)–(7) of
[the activation recheck](../core_activation_volterra_recheck/README.md):
XW=w and (E/H)V=V/r. In particular both pressure terms survive.
The only reciprocals in (1) are geometric; they remain bounded on
X≥Xa>0, |η|≤1. Small swirl and the flat stress at Xb cost no field
inverse in this estimate. The exact global identities below are needed
for the reconstruction, not for the algebra of (1).

## The heat endpoint supplies the whole incoming trace

Set B=Xb, the actual upper stress endpoint. For X≥B, U=0,
E=c X^(-A) f(z), z=2(1−η²)/X, k=√2 c. The heat Gamma integral
in the [global moment recheck](../core_global_pressure_moment_recheck/README.md)
gives f(0)=1 and, for z≥0,

    |f^(j)(z)|/j! ≤ a_j=(h)_j(1+h)_j/j!,   a₀=1.

Here (h)_j is rising factorial. For an angular jet use the sum of
absolute normalized Taylor coefficients through order l at a real point.
No convergent Taylor series at the heat endpoint is assumed. For X≥B
the norm of f−1 is at most C_l/X, where

    C_l=2a₁+Σ_(j=1)^l a_j 6^j B^(1−j),    G_l=1+C_l/B.

The value follows by the mean value theorem. At positive degree the
increment of z has coefficient sum ≤6/X (linear and quadratic in
the angular increment), so finite Taylor composition proves the sum.
All bounds include one-sided jets at η=±1. Integrating termwise is
justified by the resulting integrable majorants, for the fixed h>0.

Use the [signed terminal identities](../global_moment_stress/README.md),
not I(∞)=0. M(B)=J(B)=0, while the angular trace norms satisfy

    ||I(B)|| ≤ |k| [B^(1−h)/(1−h)+C_l B^(-h)/h],
    ||S(B)|| ≤ c² G_l² B^(-2h)/(4h),
    ||Π(B)|| ≤ c² G_l² B^(-1−2h)/(2(1+2h)).                 (2)

Indeed I(B)=k B^(1−h)/(1−h)−∫_B^∞ kX^(-h)(f−1)dX;
S(B)=∫_B^∞ E²/2 and Π(B)=−∫_B^∞ E²/(2X), with dX understood.
For l=0 the subtraction in I is already necessary for integrability.
The constants need not be uniform as h decreases to zero.

## A common bound on the finite hull

On [Xa,B] reconstruct each primitive P from P(B)−∫_X^B g, taking
the five densities in order M,I,J,S,Π:

    U,   rE,   rEU,   U²−E²/2,   E²/(2X).                  (3)

In particular (3) integrates through all intervening modulation and
repair supports. A locally zero velocity does not set a primitive to
zero. Let g_ij bound the absolute (i,j) Taylor coefficient in (X,η)
of its density uniformly. If t_l is the relevant bound (2) or zero,
a coefficientwise bound for P is

    P_0j ≤ t_l+(B−Xa) g_0j      (j≤l),
    P_ij ≤ g_(i−1,j)/i          (i≥1).                      (4)

These are bounds on magnitudes. The first row uses a single angular
coefficient of the endpoint trace, the second just the fundamental
theorem of calculus. In particular high radial orders do not propagate
to the trace. To evaluate (1), multiply absolute Taylor tables by
ordinary convolution; differentiate via (i+1) times entry (i+1,j)
or (j+1) times entry (i,j+1). Use r and its reciprocal from the
binomial coefficients of (1+ΔX/X)^(±1/2), with prefactors bounded by
√(2B), 1/√(2Xa) and X⁻¹≤Xa⁻¹. Use the geometric expansion for 1/X;
for 1/L use value bound 1/(1−2h) and incremental majorant
(4h z+2h z²)/(1−2h) in its finite geometric series (z is the
formal η increment). Apply signs only to the exact formulas.
This gives explicit finite positive sums from raw leading E,U
derivatives (divide by i!j! before multiplication). No supremum of the
stress itself, complex heat strip or unknown axis-pressure jet is used.

For Q57's slow restriction let d₀=n+99997, e=n+1. In the chart
Tβ=s^(-1-h) T₀(X,η), X=R²/(2s), η=Z s^(-D), each derivative on
η uses a Z or T slot. Radial changes R↔U_chart do not change this.
Thus the outer target table needs i+j≤d₀, j≤e. This restriction also
survives Leibniz in Q59's slow column solve: subindices of a requested
derivative never add axial or time slots. In (1) a differentiated
primitive can add one η slot, hence (2) needs only l=n+2. Local
velocity enters through degree ≤d₀+1=n+99998, with angular order
≤e+1. The chart factors are the original finite Taylor compositions
in Q58 and the target budget; the same smooth zero extension applies
at both stress edges. A coarser isotropic velocity bound through
n+99998 as in Q61 therefore suffices simultaneously for this target
and the high generator inputs on their respective domains. This
angular saving is asserted for the target and its incoming memory;
the entire ODE is not assigned that restricted budget here.

Next enclose the leading velocity derivatives themselves, region by
region, and the low-order physical value/propagator and correction
bounds on the accepted bands. Q62 is a conditional target estimate,
not a numerical enclosure of the finalized velocity or an independent
proof of the exact repair. The radial pressure tail, uniform physical
pulse and full residual estimates remain open.

Reproduce with test_audit.py and audit.py --verify-record evidence.json
--output /path/to/fresh.json. The exact manufactured checks test the
algebra, derivative normalization and terminal majorant.
