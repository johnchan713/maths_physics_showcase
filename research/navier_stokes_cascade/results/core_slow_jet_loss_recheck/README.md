# Exact slow-map potential moments and the remaining target-jet inequality

Parent: `51d9ae01a20fa3d2d3fe00f6ae2f8aa6f021a833`. The preceding checkpoint
supplies uniform reserved geometry and compensator jets. This checkpoint
identifies exact NONZERO coefficients multiplying higher incoming target
jets in the prescribed axial slow inverse. It gives a precise obstruction
to estimating this isolated repair remainder in the same finite-order
target norm, rather than silently treating a small epsilon as a contraction.

This is not a refutation of the actual canonical sequence. Its current
targets depend on the incoming field; they are not freely prescribed while
that field stays fixed. The retained smooth finite-cycle estimates already
allow higher input jets and stage-dependent constants. They are preserved.
No unconditional pulse bound, complete residual improvement or blow-up is
claimed. Historical evidence and the common tail, width and tilt are unchanged.

## 1. The prescribed two-bump inverse has a logarithmic potential moment

Use exactly the geometric axial block of `compact_mean_foundation`. Let
g=exp(d)>1, phi0(x)=eta0(x), phi1(x)=g^-1 eta0(x/g), with a fixed
nonnegative nonzero smooth compact eta0 and disjoint copies. The entire
potential hull lies in I4. Write mu_p=integral x^p eta0(x) dx and
q_lambda=g^(-2lambda), lambda>0. The retained coefficient a(eta) is positive.
The zero-flux and current angular rows give

    gamma_q=s1(phi1-g phi0),
    s1=Jtheta_q/[a(eta) mu_(1-2lambda) g(1-q_lambda)],
    integral x gamma_q=0, integral x^2 V_q gamma_q=-Jtheta_q.

Signs agree with the CURRENT five-map, not a historical target. Set
Psi_q=x^-1 integral_left^x y gamma_q(y) dy. It is compact because its
weighted source total is zero. Integration by parts, with no boundary term,
gives the exact identities

    integral Psi_q = -integral x log(x) gamma_q(x) dx,
    integral x^2 Psi_q = -(1/2)integral x^3 gamma_q(x) dx.

Under a geometric dilation the logarithmic moment gains g mu_1 log g;
the two ordinary log moments cancel. Consequently

    integral Psi_q = -H Jtheta_q,
    H=mu_1 log g/[a mu_(1-2lambda)(1-q_lambda)] >0,
    integral x^2 Psi_q = -H2 Jtheta_q,
    H2=mu_3(g^2-1)/[2a mu_(1-2lambda)(1-q_lambda)] >0.       (1)

In particular zero axial flux does NOT imply zero unweighted potential
moment. The potential fills the radial gap between the bumps; dropping
that gap loses the log g contribution. For fixed profiles,

    mu_1/[2a mu_(1-2lambda)lambda] <= H
      <= mu_1 exp(2d lambda_bar)/[2a mu_(1-2lambda)lambda]

when 0<lambda<=lambda_bar, and lambda*a*H tends to 1/2 as lambda tends
to zero after factoring out a. Thus this moment, not just the pointwise
inverse amplitude, has the unavoidable fixed inverse-tilt cost. No uniform
bound for the construction as lambda varies is inferred.

## 2. Transfer to the SAME Q-chart before differentiating

Keep A=1/2+h, s=q/Q, x=R/sqrt(s). The physical row normalization is
Jtheta_q=s^(2A-3/2)Jtheta, and the mean velocity acquires s^-A. Equivalently
the Q-chart power coefficient is a_eff=a s^(lambda-h), and a normalized
dilated bump has moment mu_p^R=s^(p/2)mu_p. Either calculation gives

    integral Psi_R = -K Jtheta,       K=s^h H,
    integral R^2 Psi_R = -K2 Jtheta,  K2=s^(1+h) H2.          (2)

All s and a dependences remain inside differentiation. Omitting the
s^h factor would contradict the physical five rows. The exact slow radial
increment is B=-epsilon partial_Z Psi_R, so

    integral B=epsilon partial_Z(K Jtheta),
    integral R^2 B=epsilon partial_Z(K2 Jtheta).             (3)

Take the linear response in the axial target direction (P,Jtheta,Jz)
=(0,J,0), so the angular increment u is zero. Hold the complete six-entry
wave covariance fixed as required during this slow operation. In the
exact recomputed remainder map, tau=epsilon partial_T+epsilon^3 partial_ZZ
acts on the moments (3). Its contributions to the new pressure defect and
weighted radial remainder are EXACTLY

    epsilon^2 partial_T partial_Z(K J)+epsilon^4 partial_Z^3(K J),
    epsilon^2 partial_T partial_Z(K2 J)+epsilon^4 partial_Z^3(K2 J).  (4)

Radial vector viscosity has zero plain and R^2-weighted integrals, whereas
axial viscosity does not. All old means and the actual b=O(epsilon) are
retained. At linear order their transport terms contain at most two Z
derivatives of J; they cannot remove the third-Z-derivative coefficient
in (4). A jet with J=J_Z=J_ZZ=0 isolates that coefficient at a slow point,
even with moving s,a and nonzero old means. For the three defects the
principal third-Z-derivative response vector is therefore

    (Pnew,Jtheta_new,Jz_new):
       epsilon^4 (K,0,-K2/2) partial_Z^3 J.                (5)

Here "principal" refers only to the highest number of Z derivatives of
this incoming target in the LINEARIZED remainder. It is not a statement
that the full finite remainder equals (5). Slow time, transport, lower
target jets, squares and all pressure compensation still remain.

## 3. Whole moving pressure can cancel one principal coefficient

For the SAME physical compensator, c_rho=s*c_hat. The complete axial
residual moment is epsilon partial_Z(Jz_new+c_rho Pnew). Its principal
fourth-Z-derivative coefficient from (5) is

    epsilon^5 [c_rho K-K2/2].                              (6)

It vanishes precisely if

    c_hat = mu_3(g^2-1)/(4 mu_1 log g).                     (7)

No assertion is made that the actual chosen rhohat satisfies (7). This is
an explicit coefficient that must be checked before claiming generic
fourth-order loss or a special cancellation in the full axial moment.
The pressure-defect coefficient K in (5) stays strictly positive either way.
Also (7) does not erase the whole pressure contribution. If K2=2c_rho K,
then its axial-viscous part becomes

    epsilon^5 partial_Z[c_rho partial_Z^3 F
                                      -partial_Z^3(c_rho F)], F=KJ,
    c_rho F'''-(c_rho F)'''=-3c_rho' F''-3c_rho'' F'-c_rho''' F.

The moving-compensator commutator remains. Pi1/Pi2 and their physical
cutoff derivatives stay in the incoming complete base equations; fixing W
or subtracting the same base balance does not replace pressure by leading p0.

## 4. Precise norm obstruction and the admissible higher-jet budget

Use the normalized slow coefficient norms of `core_dynamic_mean_recheck`.
Leibniz and three derivative shifts give the sufficient bound

    |epsilon^4 partial_Z^3(KJ)|_m
      <=epsilon^4 (m+1)(m+2)(m+3) |K|_(m+3)|J|_(m+3).       (8)

For the whole axial-viscous moment a corresponding bound is

    epsilon^5 (m+1)(m+2)(m+3)(m+4)
      [|c_rho|_(m+1)|K|_(m+4)+|K2|_(m+4)/2]|J|_(m+4).    (9)

The preceding reserved-coordinate certificate supplies the s and c_rho
costs in these expressions; the a^-1 jets, fixed profiles and actual current
target jets must still be controlled. For proposed tolerances T_m the
specific remaining sufficient inequalities are that (8),(9), TOGETHER
with the retained lower-order transport and nonlinear budgets, fit T_m.
They are not replaced by an exponent table or by low-order smallness.

At fixed epsilon>0 the isolated target-to-pressure remainder map with
incoming coefficients frozen is not locally Lipschitz from C^m to C^m.
Indeed choose smooth slow cutoffs equal to one near an interior point and
J_n=n^-m chi sin(n(Z-Z0)+theta). Their C^m norms stay bounded, but the
m-th Z derivative of its linear response has leading size epsilon^4 K n^3;
lower-order terms are O(n^2). Choose the phase to realize the leading term.
Thus a same-order operator bound would require

    epsilon^4 K n^3 <= constant + O(n^2),                  (10)

which fails as n tends to infinity. Taking an infinitesimal amplitude first
isolates the linear response for each n; the argument does not assume that
large finite corrections preserve any field budget. This is a continuum target-data control,
NOT a family of prescribed canonical current states. It refutes a same-norm
contraction inference for the isolated map, not the retained smooth finite
cycle, whose actual targets and old fields are linked and have higher-jet
and cancellation hypotheses. Actual polynomial band-frequency budgets may
absorb the epsilon factor; those budgets must be proved on the same common
tail. No new tail, derivative-dependent first band or profile is chosen here.

Next prove the actual current-target bounds in (8),(9), or identify the
complete cycle cancellation that replaces them, retaining all lower terms,
wave supports and pressure. This is the sharpened input obligation before
claiming a full pointwise residual gain.

## 5. Diagnostics and limits

Piecewise exact polynomial primitives independently integrate the radial
gap and both potential moments. Logarithms are formal log g/log b symbols
whose coefficients cancel exactly; no float evaluates a tiny tilt. A
manufactured lambda=1/4 block has rational fractional-power moments and
checks the actual two linear rows. Independent plus/minus evaluations of
the complete averaged radial operator isolate its linear response with
nonzero b and old means, including axial viscosity. Dropping that viscosity,
the gap or the physical normalization is detected. The algebraic condition
(7) and moving-pressure commutator are checked without claiming an actual
compensator choice or a full incoming state realization. Finite endpoint
fixtures are not the actual C-infinity profiles; continuum arguments are
(1)--(10), not tests or finite frequency samples.

```sh
python3 research/navier_stokes_cascade/results/core_slow_jet_loss_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_slow_jet_loss_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_slow_jet_loss_recheck/evidence.json \
  --output /tmp/core_slow_jet_loss_recheck_fresh.json
```
