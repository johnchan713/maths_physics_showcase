# Q23: finite-band completion of the uniform upper propagator

Parent: `30725005ddb9aaa93737238451f757fc473df5c0` (Q22).
Q22's comparison with the original selector constants is sufficient, not
necessary for a uniform UPPER propagator bound. Here a finite-head/tail
argument removes that comparison from this particular conclusion. It does
not remove it from an assertion about the original selector's numerical
bounds, primary cone, Haar-column invertibility or signed root.

## 1. Keep the field and separate an analytical threshold

Use the same fixed h>0, radius r0=1/6057072, enlarged positive radial
hulls, leading spectral floors and normalized Borel schedule. Q22 derives
band-independent Db_req and Mz_req from the finite-head background theorem.
Insert these as NEW ESTIMATE constants into the frame formulas, keeping
all actual fields, phases, lengths, entrances and cutoffs unchanged. Call
the resulting constants CE*, Cd*, Kframe* and let

    Lreq=max(1, normal_acceptance*, Riccati_acceptance*).

The Riccati threshold is harmless but unnecessary for the upper propagator.
It is included only to reuse the complete existing sufficient value gate.
Neither these constants nor this analysis threshold replace an input to
the construction or identify a different selected pulse.

The original geometry gives S=ell^2, 2r0 S<=L<2r0 Tg S<=S. The last
inequality follows from Tg=4+sqrt(2)<6 and the fixed radius. Also
k=ceil(epsilon^-1/2), epsilon=2^(-h ell), so

    S^2(epsilon+epsilon^2+k^-1)
       <=3 ell^4 exp(-h log(2) ell/2)
       <=360 (4/h)^5 / ell.                            (1)

Indeed log(2)>=1/2 and exp(x)>=x^5/120. A sufficient integer threshold is

    ell0=ceil(max(16, 360(4/h)^5, Lreq/(2r0))).          (2)

Since ell>=1, (2) gives 2r0 ell^2>=Lreq. This is a finite exact-real
recipe at the original h, not h=1/200 and not a reselected accepted band.
For all existing bands ell>=ell0 the original frame proof now applies
with the starred estimate constants and gives

    ||V_T,1(v,w)|| <= Ttail P(v)/P(w),
    Ttail=Kframe* exp(CE*+Cd*).                         (3)

No Q22 original-constant dominance is used in (3). The old constants and
all historical evidence remain unchanged. All fixed background and
leading-floor hypotheses of the frame provider are retained.

## 2. The finite earlier set has an independent normal floor

Consider only existing accepted bands with ell<ell0. This set is finite,
even though ell0 can be enormous. On each band there are finitely many
slow labels, torus charts and deck copies on the closed active hulls.
Continuous parameters range over the same compact enlarged charts.
Frozen nearest-NONZERO rounding gives p=j/k, j a nonzero integer. Thus

    |n_tan|>=|n_theta|=|p|/R>=1/(k Rmax)>0.             (4)

This bound was inadequate as a uniform asymptotic floor; it is adequate
on a FINITE set. Let khead be the maximum k on that set, Rmax the maximum
of its finitely many fixed hull bounds, and set eta=1/(khead Rmax).
Then |n|>=eta. No positive lower bound on another normal component,
actual-reference comparison or canonical-root identity is required.
If the head is empty, omit its constants and use (3) alone.

Let K0, N1, d0 and p0 bound ||Kbg||, |n'|, d=epsilon k^2|n|^2 and
|P'/P| on this finite union; Lhead bounds its lengths. These are finite
by smoothness, the complete background bounds, (4) and frozen-label
semantics. The pressure-eliminated m=1 ambient generator is

    A1=-Kbg+n(n^T Kbg-nprime^T)/|n|^2-d I,
    ||A1||<=2K0+N1/eta+d0=:Bhead.                     (5)

It preserves tangency because nprime^T+n^T A1=-d n^T.
For y=t/P the generator is A1-(P'/P)I. Gronwall on EVERY subinterval
of the finite head therefore gives

    ||V_T,1(v,w)||<=Thead P(v)/P(w),
    Thead=exp[Lhead(Bhead+p0)].                        (6)

This uses the logarithmic Gaussian rate, not a minimum of P on all bands.
The finite-head constant may be extremely large; it is independent of the
running band. Compactness is used ONLY after the band set is finite, fixing
the quantifier error in Q3. One cannot use (6) on the entire infinite tail.

## 3. Uniform conclusion, harmonics and derivative scope

Set T=max(Ttail,Thead) (or T=Ttail for an empty head). Then on every
original accepted band, every lift and 0<=w<=v<=L,

    ||V_T,1(v,w)||<=T P(v)/P(w).                        (7)

The estimate is a constant, hence polynomial in S. It proves existence
of a uniform upper bound relative to the retained background/geometry
and leading value hypotheses; no numerical materialization of ell0 or
finite-head suprema is needed for that analytical conclusion.

For all nonzero integer m, A_m=A1-(m^2-1)d I, so exactly

    V_T,m(v,w)=exp[-(m^2-1) integral_w^v d] V_T,1(v,w).

This retains nonnegative additional ACTUAL damping without comparing
m^2 d to m^2 dref. The uniform damping floor can be the minimum of the
starred tail floor and eta^2, since epsilon k^2>=1. The complete normal
reciprocal is uniformly bounded by the maximum of their corresponding
reciprocals. Differentiated inverse recipes still need their own full
coefficient/source jets and adapted normal connection; (7) is not a
shortcut that drops those terms. Finite-head fixed-order coefficient
jet suprema are finite by the same argument, separately at each order.
No constant uniform in derivative order is asserted.

In particular the unprojected normal source remains in
 pi=i(n.Kbg.t-nprime.t+n.f)/(km|n|^2).
The mean pressure primitives, Pi1/Pi2 memory, moving c_rho P moment,
axial A1 velocity and all covariance entries are unaffected.

## 4. What this does NOT repair

An upper ODE bound and n!=0 do not establish a positive canonical pulse
component, invertible covariance column matrix H, or a cone margin.
For example a constant two-dimensional rotation generator on [0,pi/2]
has propagator norm one while its first component starting at e1 reaches
zero. Compactness cannot turn that component into a positive floor.
Thus a finite-head completion of any original-root lower bound would
require an independent proof of nonvanishing/positivity on that head.
It is not obtained by applying (6) backwards or by enlarging T.

Removed for the UPPER tangent estimate: Q22's sufficient dominance tests
Db_req<=Db_original and Mz_req<=Mz_original. Not removed: all retained
actual-background hypotheses, root/cone/column lower estimates, source
norms, exact same-field identification and complete residual obligations.
Next audit the canonical Haar-column/root lower bounds on the original
accepted domain; they are the next possible obstruction when composing
(7) into the first signed-prefix source estimate. No unconditional pulse,
complete residual improvement, endpoint theorem or blow-up is claimed.

## 5. Reproduction

The code retains exponentials as exact symbolic recipes. Rational controls
verify the tail threshold and rounding floor, the two-sided frame/head
budget composition, the moving-normal constraint identity and harmonic
scalar damping. The rotation control prevents treating an upper estimate
as a positive component floor. They are manufactured algebraic controls,
not evaluations of the actual giant finite head. Parent frozen evidence
and source hashes are replayed.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json` in this directory.
