# Q29: explicit weighted constants for the shifted compact maps

Parent: `d218ebe560bc54586a7eac536b7267b981153fd1` (Q28).
This checkpoint combines the existing two-edge kernel theorem and raw
compact-map jet estimate into explicit weighted T_e/A_e constants at each
pole in Q28's archive. It keeps the logarithmic-to-characteristic Jacobian
and the actual cutoff plateau gap visible. No field, cutoff, h, width,
pressure gauge or historical evidence changes. The result is analytical
and conditional on the stated actual fixed-profile data, not a numerical
certificate of the selected profiles or complete source.

## 1. Coordinates and fixed data

Let x=log(X/Xa), 0<x<Lx, R=sqrt(2sXa) exp(x/2), U=R^d,
s in [1/4,4], d>0 fixed, and

    zeta=exp(-a/x^2-4/(Lx-x)^2), a=t1^2,
    delta=min(1,x,Lx-x), w_p=zeta delta^-p.

The source coefficient derivatives are in (U,Z,T,y) with M and all label
parameters frozen. Suppose each input derivative through m is bounded by
A w_p. This is the weighted INPUT hypothesis; the theorem does not infer
it from an unweighted bound. Fixed native/log chart conversions are applied
before this hypothesis is used.

Retain the original moving cutoff chi(U,Z,T). Require chi=0 for x<eta
and chi=1 for x>Lx-eta, with ONE eta>0 uniform in the slow chart.
Its transition and every positive-order derivative are then supported in
[eta,Lx-eta]. These are actual moving support statements. They follow
for a fixed transition profile with that support, but not from pointwise
compactness of an arbitrary band-dependent cutoff family. The precise
remaining support test, if the actual profile data are not supplied, is

    inf_over_bands,slow_parameters, supp(d chi)
                    min(x,Lx-x) >= eta>0.              (1)

Also keep the original zero/one plateau values. A transition support bound
alone does not specify them. No new cutoff is chosen to satisfy (1).

Let f_m,g_m,k_m,chi_m,(1-chi)_m,(chi_U)_m be raw through-order multiplier
bounds for f_e=R^(e-d+1)/d, g_e=R^-e, k_e=d R^(d-1-e) and the actual
cutoff. The last input uses chi through m+1. On the positive fixed hull
these are finite fixed-profile/chart bounds. Define

    Jmax >= sup dU/dx = sup d U/2,
    r=min(1,Lx/4,eta),
    E_p >= sup zeta delta^-p,
    Z_eta=exp((a+4)/eta^2).                            (2)

For example the fixed bound Umax=(8Xa)^(d/2)exp(d Lx/2) gives
Jmax=d Umax/2 on the whole scale chart. This constant has no running-band
or M dependence. The existing Eedge(p), which bounds sqrt(zeta)delta^-p,
is also a valid E_p. Giant exponentials remain symbolic finite constants.

Use the signed-target edge provider's three-branch constant C_p with its
proof collar reduced to r. For either directed logarithmic integral,
it bounds its weight ratio by C_p on the corresponding endpoint collar
and on [r,Lx-r]. It is NOT a bound for a left primitive near the right
edge without cancellation. The cutoff is what excludes that region below.

## 2. Weighted estimate without differentiating the large shift

Keep exactly

    I_- = int_(a_shift<0) (f_e f)(U+a_shift,y+M a_shift vr) da_shift,
    I_+ = int_(a_shift>0) (f_e f)(U+a_shift,y+M a_shift vr) da_shift,
    T_e f=g_e[(1-chi)I_- -chi I_+],
    A_e f=k_e chi_U int_all (f_e f)(U+a_shift,y+M a_shift vr) da_shift.

For each coefficient derivative, differentiate at FIXED a_shift. The
translation is independent of every differentiated variable. No M factor
appears. Smooth zero extension removes moving support boundary terms.
One must not first change to x and then silently discard derivatives of
its moving endpoints or Jacobian. The change of variable is used only
AFTER differentiation, to estimate the absolute integral.

Each differentiated source obeys the same assumed A w_p at its input
radius. The torus shift changes no radial weight. In absolute integrals,
da_shift=dU'=J(x')dx', bounded by Jmax dx'. Apply the two-edge theorem
to the left branch only where a derivative of 1-chi is nonzero (or its
zeroth-order value is nonzero); that region excludes the right collar.
Apply its backward counterpart to the right branch, which excludes the
left collar. Leibniz with four factors costs at most 2^(3m). Hence

    |D^I T_e f| <= A w_p KT_(e,m,p), |I|<=m,
    KT_(e,m,p)=2^(3m) g_m[(1-chi)_m+chi_m] f_m Jmax C_p. (3)

There is no additional Lx factor in (3): C_p already bounds the directed
logarithmic integral. Nor may the characteristic length Lu be substituted
for Jmax without a separate inequality. The unweighted compact estimate's
Lu and this weighted estimate's Jmax C_p are alternative integral bounds,
not identical factors to multiply indiscriminately.

For A_e, the whole integral is at most A f_m Jmax Lx E_p. Its output
and all differentiated output multipliers lie in the actual transition
support. There zeta^-1<=Z_eta and w_p^-1=zeta^-1 delta^p<=Z_eta.
Consequently

    |D^I A_e f| <= A w_p KA_(e,m,p),
    KA_(e,m,p)=2^(3m) k_m(chi_U)_m f_m Jmax Lx E_p Z_eta. (4)

This retains the full shifted moment. Haar zero does not make it zero.
For auxiliary-independent input with its exact weighted radial moment
zero, the A_e output vanishes by the separate exact identity; that special
case is unnecessary for (4). T_e can fill radial gaps; A_e can fill the
cutoff transition. Both retain smooth zero extension and the slow support
projection. The floor in (4) is never applied on a rectangular union of
moving radial supports, which can reach outside the shell.

## 3. Bind the two archives and pressure

Apply (3)--(4) independently to e=0,1,2 as needed, using their actual
multipliers. Q28's fast poles are {2p,2r}, and its pressure/moment poles
are {2p,4p,2r,4r}, where r=max(2p,p_signed). Retain each distinct constant;
their sum supplies a common padded constant without any monotonicity
assumption in pole. Constants are also padded over the finite requested
orders if a through-order provider is needed. No numerical enumeration
of the full giant Q14 request graph is required to prove finite existence.

For pressure use e=0 and the FULL compensated input gr-rho Pdef. The
fixed interior rho has its own weighted profile bound from Q4. Its moving
Pdef and the entire derivative of c_rho Pdef remain. In the axial fast
velocity use gd-A1 gd with (4), not gd alone. The temporal Fourier bound
from Q28 then composes with these constants at the original m+4/m+5
requests, with no new source derivative loss or q-power loss.

This removes unspecified weighted compact operator constants relative to
actual multiplier jets, plateau gap, positive hull and input weights.
The slow five-response profile constants and its whole potential-hull
support, full background/weighted stress schedule and canonical premises
remain separate. Next bind those slow-response data and (1) to the exact
fixed choices, then complete the first-source input certificate. Full
residual improvement, localization and blow-up remain unproved.

Tests check the coordinate Jacobian, directed-kernel domain, support-floor
logic, shifted derivatives and distinct archive padding on exact fixtures.
They are not continuum certification of selected profile suprema. Frozen
Q28 evidence and its dependency graph are replayed. Run `python
test_audit.py -v` and `python audit.py --verify-record evidence.json
--output FRESH.json` here.
