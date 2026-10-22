# The next multiharmonic source and the mean-feedback gap

Parent: `6b7d106fbe9b5811c41a9e617808bd16ba6f6f05`. This checkpoint
joins the parent's generated inverse and fresh mean SIZE bounds to a
quantitative **next regular nonzero residual**. It removes the missing
single-frequency interface: the input wave has +/-1,+/-2 and the output
can contain +/-1,...,+/-4. The resulting bound has the old coarse order
B=1-kappa. The desired improved order D=3/2-2kappa does **not** follow
from those mean-size bounds. We identify its precise remaining inequality.

The general finite-harmonic identities and the half-power mean-wave loss
were already proved in `regular_source_foundation`,
`finite_induction_foundation`, and `current_residual_ledger_audit`.
We do not claim those as new results. The new work is their explicit
constant/polynomial composition with the current source-generated wave
and fresh means, and the identification of the missing bound in THIS
composition. All historical evidence and actual fields remain unchanged.

## 1. Exact current split and full pressure

Let Bbase denote the same complete background, w the old complete cut
wave, m the accumulated old mean, and f the CURRENT regular nonzero
source before the sourced solve. Write the old residual G=f+Told.
The sourced inverse retains its full normal pressure, including n.f;
its exact principal identity gives

    L_Bbase(v,pi)=-f+r.

Here r includes the full regular linear remainder on the principal
amplitude AND the principal operator on its coefficient curl. For the
retained temporal pulse cutoff psi, vpsi=curl(psi A). With all spatial,
slow and transverse localizers already in A and pi, the exact new
nonzero split after the wave and then the two fresh mean increments is

    Greg_wave=psi*r+Pperp[N(w,vpsi)+N(vpsi,w)+N(vpsi,vpsi)
                             +N(m,vpsi)+N(vpsi,m)],
    Tnew=Told+(1-psi)*f+Tpsi,
    Greg_final=Greg_wave+Pperp[N(d,w+vpsi)+N(w+vpsi,d)],
    d=d_fast+d_slow.                                    (1)

N is full cylindrical advection, including connection terms; Pperp is
the nonzero ANGULAR projection. It is not the auxiliary Haar projection.
Both old/new crosses, the new square, and BOTH mean-wave orientations
are present. Nonlinearity acts on the actual complete cut velocities,
including old compact A0/A1 velocities, not on the additive tail residuals.

The complete pressure is the SAME p_B (with Pi1 and Pi2), the current
reconstructed angular-zero p_m, and all old/new nonzero wave pressures.
The new nonzero pressure is retained in r and Tpsi. Recomputing p_m
changes the zero-mode equation; its gradient has no nonzero angular
component. That projection fact neither deletes p_m nor proves its
zero-mode residual small. `NextResidualBudget.regular_source` explicitly
calls the parent's postslow full-pressure/row reconstruction as well.

For a general cutoff the exact physical commutator is

    Tpsi=(psi_t+Bbase.grad psi)v-2(grad v)grad psi
          -(Delta psi)v+pi grad psi+L_Bbase(grad psi cross A,0). (2)

The retained Gaussian-ledger psi is temporal, so its spatial terms vanish.
Spatial/localization cutoff derivatives instead remain inside the UNCUT
coefficients and regular remainder; they cannot all be declared Gaussian
tails. The old exact tail ledger is preserved, unestimated here, and is
never sent to an inverse. Conditional Gaussian tail bounds already in
the ledger audit still require their original separation/domain inputs.

## 2. Explicit finite-harmonic constants

Use harmonic l1 factorial through-order coefficient norms; products are
submultiplicative, with finite graph/pole losses already retained by the
parent. For derivative order n put d0=n+4, dr=g*d0 and

    A=dr+d0+r0,
    H=max |m|=2, Rratio=max |m'/m|=2.                 (3)

r0 bounds the reciprocal radius; N0 bounds the full normal. All inputs
below are bounded through n+3. The full linear remainder polynomial uses

    Eslow=d0+b*dr+G*d0, Ephase=2H*Eik,
    Ebase=bR+Bz+b*r0,
    Evis=d0^2*(g^2+1)+r0*dr+r0^2,
    Emixed=4H*N0*(dr+d0)+2H*N0*(dr+d0+r0),
    Econnection=4H*N0*r0, E=sum of these six groups,
    Flin=psi_norm*[E*(T+R)+(dr+d0)*Pi
                         +(c*d0+K+4H^2*N0^2)*R].       (4)

T,R,Pi here bound the **uncut sourced increment**, respectively its
principal amplitude, coefficient curl and FULL pressure. They are
generated from the same constrained inverse, including the direct
source-normal pressure term. The coefficient curl includes cylindrical
connections. The norm psi_norm is applied after its full product rule.
Using already-cut coefficients here would not bound the uncut r in (1),
and must not be justified by division by psi.

Equation (4) follows from the full operator (1)--(2) of
`regular_source_foundation`: k<=2 epsilon^-1/2 gives the factors 2H
and 4H, while epsilon k^2 m^2<=4H^2 gives the principal curl damping.
Only the principal constrained solve cancels f. There is no homogeneous
cancellation for the coefficient curl or a truncation of its pressure.

For a complete input amplitude a_m, exact divergence gives

    div_*a_m+i k m n.a_m=0,
    i k m' (n.a_m)a_m'=-(m'/m)(div_*a_m)a_m'.          (5)

There is no division by m+m', including when it is zero. Coefficient
advection plus connections cost at most A; (5) costs Rratio*A. Thus

    Cww=(1+Rratio)*A=3A,
    Cmw=2A+2H*N0.                                    (6)

The mean phase term genuinely costs k. It has no analogue of (5) with
the mean as advector. Convolution in harmonic l1 introduces no independent
mode-count constant. Zero outputs are kept in the full mean system;
projecting them out cannot increase the l1 bound for nonzero outputs.

Let Uold,Unew be norms of COMPLETE CUT old/new velocities; Mold bounds
the old accumulated means and Md bounds d_fast+d_slow. Then explicit
positive polynomials for (1) are

    Fwave=Flin+Cww*(2Uold*Unew+Unew^2)+Cmw*Mold*Unew,
    Ffeedback=Cmw*Md*(Uold+Unew),
    Ffinal=Fwave+Ffeedback.                            (7)

This retains both old/new transport orientations. Formula (4) scales
correctly at |m|=2; neither the old m^2=1 constant nor the old ratio-one
transport constant suffices. For H={-1,1}, its linear part reduces
exactly to the earlier single-frequency recipe.

## 3. The current one-way polynomial composition

`NextResidualBudget` takes ONLY the parent's fresh-mean object. It obtains
uncut T,R,Pi from the same generated inverse, cut Unew from its complete
cut curl, Mold from the original prefix's final means, and Md from the
sum of the actual fast/slow increment majorants. In particular it does
not subtract two unrelated norm bounds or replace Md by zero. The slow
increment uses fresh postfast pressure and scalars as in the parent.

Raw inverse coefficient jets are converted to factorial through-order
norms on the same five variables. The temporal cutoff receives the same
conversion before (4); raw source jet j is bounded by j!*Ffinal_j.
Sparse polynomial hulls only increase estimates for x>=1. Every output
polynomial is evaluated at the parent's x=S+K, an estimate argument only.
No physical scale, width C^-20, band, carrier or support is changed.

For a common original finite pole p covering the requested jets, the
old wave has pole p, sourced wave 3p, old means 2p, and fresh mean
increments at most 6p. Wave crosses/squares have poles 4p/6p, old-mean
crosses 5p, and fresh mean feedback 9p. Thus 9p is a sufficient NEXT
source pole. This does not reuse an old 3p proof domain. Products retain
sqrt(zeta)*P because zeta,P<=1; fixed pole enlargement does not change
the actual width. The order-zero composed source reads original mean
inputs through order 31. That finite dependency is not an actual-input
enclosure, an all-stage constant, or a diagonal convergence theorem.

## 4. The precise gain that is still missing

Keep B=mu=1-kappa, C=B+1/2, D=C-kappa. The SOURCED VELOCITY is W_B;
its pressure and potential have order C. Confusing velocity with pressure
would incorrectly improve every interaction. The old wave is W_(1/2).
The regular source foundation and (1) give

| Group | Sufficient order |
|---|---:|
| Complete increment linear remainder | D |
| Both old/new wave crosses | D |
| New wave square | 2B-kappa >= D |
| Old mean/new wave interaction | 2B-1/2 = D |
| Fresh mean/current wave feedback, from coarse mean SIZE alone | B |

Hence (7) supplies a conditional supported regular source in W_B with
explicit polynomial coefficients and the unchanged exact tail split.
It does **not** supply W_D. The raw-source/inverse interface exposes
only this coarse B bound; it does not certify any improved iteration.

The exact missing inequality, with d=d_fast+d_slow and complete current
wave wcur, is, for every fixed n on the original accepted band,

    ||Pperp[N(d,wcur)+N(wcur,d)]||_(n, half-weight,9p)
                              <= epsilon^D Q_n(S).      (8)

One sufficient route is to prove the actual zero-Haar post-wave fast
row inputs Q1(Etheta,Ez) have order D, and the actual fresh postfast
P,Jtheta,Jz have order C (with the extra Jtheta derivative needed by
the radial slow potential). The weighted maps would then give fast
increment order D and slow increment order C, so every row in the table
meets D. These are **unproved sufficient targets**, not new assumptions
asserted as facts. Correlated cancellation in the combined feedback (8)
could also suffice without those individual size orders. The old size
polynomials and scalar norm sums alone prove neither route.

This is the precise point where the more general historical conditional
cycle uses sharper actual mean-map/target estimates. It cannot be promoted
by the newly closed coarse mean-size calculation alone.

## 5. A sharp class-level witness, not a refutation of the actual maps

On a fixed interior radial annulus choose a nonzero smooth compact
potential chi(R), f=(partial_R+1/R)chi, and a(R) overlapping f. On a
slow/temporal plateau use phase Phi=theta+Z/epsilon, n=(0,1/R,1), and

    d=epsilon^B f(R) e_z,
    w=epsilon^(1/2) a(R) (0,R,-1) cos(k Phi).           (9)

Both fields are exactly divergence free. The mean has zero angular
velocity moment and integral R*d_z=epsilon^B[R chi]=0. The wave is a
complete curl of a compact radial potential with its 1/k factor; its
coefficient divergence and normal component are both zero. Smooth
slow/transverse and temporal extensions can be made at the potential
level, with a plateau preserving (9) and the correct extra epsilon in
the mean radial velocity. Choose supports inside the overlapping allowed
interior; the weight floor there is fixed and positive.

Exactly on that plateau,

    N(w,d)=0,
    N(d,w)=-epsilon^(B+1/2) k f*a*(0,R,-1) sin(k Phi).  (10)

Thus its norm is at least c epsilon^B since k>=epsilon^-1/2.
Demanding the D bound from these size/support/divergence/moment properties
alone would require epsilon^(-(1/2-kappa))<=C Q_n(S). Along the retained
epsilon=exp(-h ell), S=ell^2, the left side is exponential in ell and
cannot be bounded by any fixed polynomial in S. The polynomial fixtures
check (10) and the moments exactly; smooth compactness and the asymptotic
non-implication are analytical arguments, not finite scans.

The mean-wave half-power phenomenon is already known; (9) makes its
specific effect on (8) explicit. These are admissible local class fields,
**not** identified outputs of the prescribed fresh fast/slow maps. This
does not refute those maps, the historical conditional cycle, or the
possibility that its needed sharper estimates hold.

## 6. Obligation ledger and reproduction

Removed: the quantitative +/-1,+/-2 source interface, explicit m^2 and
input-ratio constants, and its one-way join to generated full pressure,
complete cut velocities and the new fresh mean increments. Added
analytical assumptions: none relative to the retained original profile,
coefficient, support and temporal-tail contracts. No new gain hypothesis
is marked proved. Remaining next: establish (8) for the ACTUAL increments,
or justify its stated sufficient row/scalar orders by exact current-target
identities. Retain all zero-mode pressure/moment equations and additive
tails when joining any eventual complete residual bound.

Original actual input/recipe validity, common-chart support coverage at
every physical scale and any necessary primary/root identification remain
conditional. No full residual improvement, unconditional pulse bounds,
infinite-field realization or blow-up result is proved.

`controls.py` checks all 16 input-frequency pairs, fourth and zero outputs,
exact divergence, the input-ratio omission, an independent coefficient-norm
transport majorant, the sharp mean witness, full-pressure/uncut-source
retention and finite composed dependency. These are rational/polynomial
diagnostics, not continuum PDE certificates. Historical exact ledger and
operator evidence is bound by source hashes and replayed separately.

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
