# Q25: direct canonical cone estimates imply the three root margins

Parent: `e3059ea9875b457244397e68a0217b4d2fcd6f44` (Q24).
The canonical direction theorem already supplies more structure than arbitrary
column enclosures. This checkpoint derives Q24's determinant margins directly
from that structure on the certified direct-cone tail. No actual column
integral needs numerical evaluation there. Earlier bands not proved to belong
to that tail retain Q24's obligation; no field or selector is changed.

## 1. Use the canonical averaged directions and their positive masses

In the orthonormal basis e1=-N, e2=K of the target plane, the normalized
canonical columns C=H/s0 have the form

    c_plus = j_plus [A(1+a_plus)e1+u(1+b_plus)e2],
    c_minus= j_minus[A(1+a_minus)e1+u(-1+b_minus)e2],    (1)

where A=|c0|sqrt(1+u^2)>0, u>0, and jmin<=j_sigma<=jmax.
The j_sigma are the positive radial-square Haar masses from
`canonical_primary_foundation` (5), with the original angular half,
transverse area and temporal cutoff retained. They are not equalized.
The direction theorem bounds the Euclidean error before division by A,u
by Cdir/sqrt(L). Hence

    |a_sigma|,|b_sigma|<=delta,
    delta=(Cdir/sqrt(L))*max(1/A,1/u).                 (2)

For T/zeta=Vhat(A e1+u r e2), the retained target cone says
Vhat>=vstar>0 and |r|<=1-gamma, 0<gamma<=1.
Vhat can diverge at a flat edge; only its lower bound is used. The same
statement holds for tau=T/|T| with Vtau in place of Vhat. Since |tau|=1,

    Vtau>=1/(A+u).                                    (3)

This follows from 1<=A Vtau+u|r|Vtau<=(A+u)Vtau.
The normalized target direction extends to the flat shell edges by Q2.
These are the original canonical columns and leading target; no historical
entrance family or full curl covariance is substituted.

## 2. Uniform oriented determinant and numerator bounds

Orient all determinants so that the ideal plus/minus column pair has
positive determinant. Divide the column determinant by j_plus j_minus A u.
Expansion of (1) gives the lower and upper bounds

    2-4delta-2delta^2 <= D/(j_plus j_minus A u)
                      <=2+4delta+2delta^2.             (4)

For the plus Cramer numerator with target tau, division by
j_minus A u Vtau gives

    1-b_minus+r(1+a_minus)
       >= gamma-delta(2-gamma).                       (5)

The minus numerator divided by j_plus A u Vtau is
1+b_plus-r(1+a_plus) and has the SAME lower bound. These formulas
retain orientation; independently taking absolute determinants would not
prove positivity.

The original direct-cone gate is delta<=gamma/4. Since delta<=1/4,

    D>= (7/8)jmin^2 A u,
    D<= (25/8)j_plus j_minus A u,
    Aplus,Aminus>=gamma jmin A u/[2(A+u)]>0.           (6)

Thus (6) explicitly supplies all three Q24 margins on this domain. The
bounds are uniform in the whole closed chart when the original direction,
mass and target-cone estimates are uniform there. They do not follow from
sampled directions or positive pulse mass alone.

## 3. Recover the original K0 and beta without losing mass conditioning

Using Vhat in (5) and cancelling the OTHER column mass against (4) gives

    b_sigma=(C^-1 T/zeta)_sigma
       >=4 gamma Vhat/(25 j_sigma)
       >=4 gamma vstar/(25 jmax)
       >=gamma vstar/(8 jmax).                        (7)

This proves the beta used by `core_primary_root_budget` (1). A generic
bound that first replaces both column masses by jmin or jmax would lose
this useful cancellation, though it could still yield a weaker floor.

For the inverse, factor C=O diag(A,u) B diag(j_plus,j_minus), where O
is the orthogonal basis matrix and B has the coordinates in (1).
Each row of B^-1 has absolute sum at most
2(1+delta)/(2-4delta-2delta^2)<=20/7. Multiplying by diag(A,u)^-1
costs 1/min(A,u), while multiplication by O^T in the induced infinity
norm costs at most sqrt(2)<2. Consequently

    ||C^-1||_infinity <=(40/7)/(jmin min(A,u))
                      <6/(jmin min(A,u)).             (8)

The basis conversion in (8) is essential: a Euclidean rotation need not
have infinity norm one. The physical inverse H^-1 additionally costs
s0^-1; (8) is for C=H/s0, exactly the root budget's normalized input.
No uniformity in varying radius, tilt or edge width is asserted.

## 4. Which accepted domain has now been checked

The canonical direct-cone condition is

    L >= [4 Cdir max(1/cmin,1/u)/gamma]^2,              (9)

in addition to its normal, reference-frame Riccati and mass gates.
Because A>=cmin, (9) implies (2) with delta<=gamma/4. The
`core_primary_tail_recheck` rule takes Lreq at least this direct length
and chooses ell_*>=ceil(sqrt(Lreq/(2r0))). Its original geometry gives
L>=2r0 ell^2. Therefore ALL ell>=ell_* obey (9), as an exact-real
inequality for the SAME budgets. No decimal evaluation of the gigantic
constants is required for this implication, and no new start band is chosen.
The other rows of that rule still supply their separate normal/slow gates.

This is a composition of the retained analytical providers, conditional on
their actual same-profile constants and domain hypotheses. It proves (6)--(8)
on their rule-defined certified direct-cone tail. It does NOT prove that
an arbitrary historical minimum band, an earlier localized field, or
Q23's finite earlier set lies inside that tail. The exact remaining domain
comparison is ell_min>=ell_* (or direct verification of the original
value gates on every remaining band). If ell_min is already DEFINED by
that rule for the same object, this comparison is tautological; one must
not confuse it with an unknown historical numerical floor.

No comparison with inherited-family roots is needed. A stronger historical
root-comparison length may remain in Lreq but is not a premise of (7).
Q24's finite-head enclosure route remains available outside this domain.

## 5. Scope and next obligation

Established: explicit determinant margins and the original root lower and
inverse constants on the canonical direct-cone tail, including physical
Haar normalization and the infinity-norm basis cost. Not established: new
actual leading/frame input enclosures or membership of every historical
band in that tail. Upper propagator estimates alone still cannot supply
these cone hypotheses.

Next resolve which minimum-band/domain the current same-field source chain
actually uses, then compose the established canonical root jets and typed
whole-path source bounds on that domain. Do not infer all-band source
bounds by silently restricting a previously localized field. Full six-entry
curl covariance, axial mean velocity, Pi1/Pi2, moving pressure moments,
A0/A1, viscosity and cutoff interactions remain in their existing ledgers.
No unconditional residual improvement, endpoint or blow-up result follows.

Exact rational controls cover all error-box corners, both Cramer numerators,
mass cancellation, orientation, the nontrivial orthogonal basis cost and
selector arithmetic. They test the analytical formulas, not the actual
continuum hypotheses. Historical evidence is replayed unchanged.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json` here.
