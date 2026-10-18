# Current nonzero-angular residual: exact pressure/cutoff ledger and fixed jets

Parent: `62cce31db067535582f9b22d4ea0727a3c158f04`.
Status: `conditional-current-nonzero-residual-ledger-and-fixed-jet-tail-bound`.
On 18 October 2026 the parent's new CI step passed; the full workflow was
still running. Historical sources and frozen evidence are unchanged.

This audits the residual complement to the parent's CURRENT target/mean
cancellation. It derives an exact nonzero-angular pressure/cutoff ledger,
retains complete cut velocities in all nonlinear terms, and gives a quantitative
all-fixed-physical-jet estimate for its designated Gaussian tails. Combined
with the reviewed operator bounds it supplies a CONDITIONAL finite-cycle
residual estimate. Those operator, initial canonical seed, background/core,
common-chart and Gaussian separation inputs are listed, not assumed proved
merely by a protocol or by these finite tests. No unconditional pulse bound,
complete residual theorem, old first-band/force certification or blow-up follows.

## 1. The exact angular decomposition removes a pressure ambiguity

Let N(a,b)=(a.grad)b, including cylindrical frame connections, and
L_B v=v_t-Delta v+N(B,v)+N(v,B), viscosity one. Write the SAME current field

    u=B+m+w,     p=p_B+p_m+pi,

where B,m,p_B,p_m are angular zero, and w,pi are the complete nonzero-angular
waves and wave pressures. P0 averages cylindrical components in the angular
variable; Pperp=1-P0. This is NOT auxiliary Haar averaging. The pulse cutoffs
and coefficient fields are angular independent, so all linear operators
preserve each nonzero angular frequency, including vector viscosity.

The nonzero residual is exactly

    G=Pperp R(u,p)
     =Pperp[L_B w+grad pi+N(m,w)+N(w,m)+N(w,w)].              (1)

In particular changing the compact mean pressure, including its actual A0
term, changes the angular-zero residual but NOT G. This is a projection
identity, not permission to drop that pressure: the full mean equations,
Jz+c_rho P, all physical gradients and later force/localization still use it.
Pi1/Pi2 cancel only against the SAME complete base balance. Leading p0 is
not substituted for p_B. A pure normal sourced wave can have zero velocity
and nonzero wave pressure; its nonzero mode remains in (1).

Thus recomputing mean pressure after a wave or mean update does not secretly
move a term between the regular NONZERO source and its tail ledger. The
zero mode is reconstructed by the parent's full current covariance/mean maps.
All six covariance entries, axial velocity and both velocity moments remain.

## 2. A complete potential/pressure cutoff, with its actual sign

For an uncut wave v=curl A and pressure pi_f, put

    vpsi=curl(psi A)=psi*v+dpsi,     dpsi=grad psi cross A,
    pipsi=psi*pi_f.

Define Tpsi by the exact identity

    L_B vpsi+grad pipsi=psi*(L_B v+grad pi_f)+Tpsi,
    Tpsi=(psi_t+B.grad psi)*v-2(grad v)grad psi-(Delta psi)*v
                                      +pi_f*grad psi+L_B dpsi. (2)

Every term is retained, including viscosity, phase differentiation inside
grad v, wave pressure and the cutoff curl. This Cartesian identity transfers
to the complete cylindrical operator, not just its principal ODE. For
general viscosity nu the two displayed Laplacian product terms carry nu;
the diagnostic below uses the equivalent normalized viscosity epsilon.

In general div(psi*v)=grad psi.v is nonzero. Exactly
div dpsi=-grad psi.v, so ONLY the complete vpsi is divergence free. Discarding
dpsi because it is physically flat can invalidate the phase-transport
identity used to obtain the wave-wave order. Nonlinear sources below use
vpsi as one complete curl, never psi*v alone. Its entire covariance is
psi^2 Cov(v)+psi Bil(v,dpsi)+Cov(dpsi), with all six entries.

The pressure-only case A=0 gives Tpsi=pi_f grad psi, generally nonzero.
It cannot be removed by principal transversality or assigned the extra
harmonic gains belonging to a vector potential. The sourced inverse keeps
pi_f=i(n.Kt-nprime.t+n.f)/(km|n|^2), including BOTH normal terms.

### Distinguish the actual temporal tail cutoff from spatial localizers

Equation (2) is valid for a GENERAL cutoff. Its being an exact commutator
does NOT make it Gaussian-flat. In the reviewed pulse construction, the
cutoff used for the ADDITIVE Gaussian ledger is the temporal factor
psi=psi0(v_pulse/Ls), with frozen-label Ls constant; v_pulse is the pulse
coordinate, NOT the velocity v. The physical evaluation has grad_x v_pulse=0,
so grad_x psi=0. Thus its actual commutator simplifies to
Tpsi=psi_t*v: dpsi, spatial viscosity products and pi_f grad psi vanish
for THIS temporal factor. Our general spatial-cutoff diagnostic deliberately
does not have that property; it tests the general identity, not the actual
value or flatness of this temporal pressure term.

The transverse chi_g(xi), slow partitions and radial root localizers are
already INCLUDED in the uncut A and pi_f of (3). Their derivatives remain
inside the complete regular remainder r, including curl corrections,
pressure-amplitude gradients and mixed viscosity. They can be nonzero
where the temporal Gaussian is maximal. There is no Gaussian tail gain
there: their regular bound is the reviewed a+1/2-kappa estimate, using
all fixed coefficient jets. Do not move their full spatial commutator into
Told merely by applying (2). For a spatial chi with chi'!=0 at a pulse
maximum P=1, a pressure coefficient epsilon^a can give a term of that
same algebraic size, not O(q^N) for every N. A proposed flatness bound
would require exp((N-ha)ell)<=C for N>ha, which fails.

This corrects a possible overbroad reading of the generic cutoff formula;
it does not refute the actual temporal construction or change historical
evidence. The independent shrinking-diagonal chi(a_j q) and later global
space/time localization cutoffs are NOT Gaussian pulse tails either. Their
own cutoff product/jet estimates must be used in the single-field/force
join; no such remainder is certified by (8) below without its hypotheses.

## 3. Exact CURRENT source cancellation, not an old residual transplant

Suppose the OLD current nonzero residual has the exact split

    G=f+Told,

where f is its CURRENT regular supported source. Invert f only, with zero
entrance data, and retain the complete uncut linear identity

    L_B v+grad pi_f=-f+r.                                    (3)

Here r is the complete regular linear remainder proved in the source
foundation, not a principal-amplitude substitution. After adding vpsi
and psi*pi_f and recomputing the angular-zero mean pressure, (1)--(3) give

    Gwave=Gnew_reg+Tnew,
    Gnew_reg=psi*r+Pperp[N(m,vpsi)+N(vpsi,m)
                         +N(w,vpsi)+N(vpsi,w)+N(vpsi,vpsi)],
    Tnew=Told+(1-psi)*f+Tpsi.                                (4)

All signs in (4) follow from f-psi*f. Old tails are NOT inverted, erased,
or multiplied by a new velocity: NS is nonlinear in velocity, not residual.
The OLD actual cut velocities remain in w,m and hence in every interaction.
If different labels have different cutoffs, apply (2)--(4) mode/label-wise
and sum; the result is linear in the cutoff ledger. Actual closed enlarged
support separation, rather than a same-phase argument, handles other labels.

For a signed homogeneous wave use f=0 in (2),(3): its full regular linear
remainder and every current interaction are added to the existing regular
source, and its Tpsi is appended to the tail ledger. It does NOT cancel
the existing nonzero source merely because its two mean fluxes match a target.

For an ACTUAL angular-zero fast or slow mean increment d, with the new wave
w fixed, pressure recomputation and all pure mean terms stay angular zero.
Therefore exactly

    Gafter_mean-Gbefore_mean=Pperp[N(d,w)+N(w,d)].              (5)

The nonzero tail ledger does not change at this mean operation. Actual
compact A0/A1 terms are not deleted: angular-zero linear terms stay in the
mean ledger, and actual A1 velocities stay in both products of (5). A mean
may fill the auxiliary torus; its product still has the wave's Gaussian P,
closed slow/transverse support and entrance/exit collars. New zero outputs
from complete wave products go to the current mean system, not a denominator
depending on their output frequency.

## 4. Close the regular budget using the reviewed full operators

Retain the corrected accepted-tail domain, including ell^6>=8*Cfreeze/gamma0,
the frozen cone margin in Lreq, the actual radius, fixed positive lambda,
width C_setup^-20, canonical fixed roots/frame and common-torus charts.
For epsilon=Q^h, kappa=10^-5, use all fixed coefficient jets with half/full
zeta weights, finite edge powers and polynomial S=ell^2 losses. No stage-
dependent first band, updated positive root or carrier is introduced.

The explicit analytical INPUTS are:

* The accepted canonical initial wave and initial mean estimates, with the
  same full base, exact divergence, both Haar velocity moments and I4 hull.
* Complete regular linear remainder W_(a+1/2-kappa) for a pulse W_a,
  including pressure, slow time, radial transport and all vector viscosity.
* Complete same-label transport W_(a+a'-kappa), derived from
  div_* a+i km n.a=0, and mean-wave transport W_(a+t-1/2) for a mean of
  tangential/radial orders t/t+1. Different-label derivative supports separate.
* Sourced inverse full-path f/P and entrance jets, all normal/frame floors,
  full pressure, signed inverse half-edge calculus, and actual compact mean
  bounds. The source/edge and current-target checkpoints supply their
  functional and current-target joins relative to these reviewed estimates.
* Gaussian separation on support of 1-psi and derivatives of psi, all-fixed
  physical jet dictionary, plus the separate full base/exterior flat bounds.

In that last input psi means the temporal pulse factor just specified,
NOT an arbitrary product of every spatial/slow/diagonal cutoff.

These are estimates/identities for the SAME formula-defined fields. They are
not supplied by a flag, a manufactured polynomial source, or positivity of
just two covariance entries. Any unproved defining/core/background premise
in those inputs remains a premise of the conclusion here.
The physical-pulse foundation is used for its uniform chart/operator
dictionary, not to transplant its historical frame, amplitudes or roots.
The canonical seed uses its OWN rebuilt roots and complete initial means
from canonical_primary_foundation and canonical_initial_mean_join.

For an input S(B), B>=B0=3/2-2kappa, put C=B+1/2, mu=rho=1-kappa.
The sourced increment and subsequent signed complete wave both have order B.
The parent's target/mean recurrence supplies actual fast mean order
H=C-kappa and slow mean order C, at freshly recomputed targets. Equations
(4),(5) exhaust the nonzero change and give these sufficient orders:

| Retained group | Regular nonzero order |
|---|---:|
| Sourced / signed complete linear remainder | B+1/2-kappa |
| Sourced / signed old-wave cross, BOTH orientations | B+1/2-kappa |
| Sourced / signed actual mean-wave cross, BOTH orientations | B+mu-1/2 |
| Each new complete wave square | 2B-kappa |
| Actual fast mean-wave change | H=B+1/2-kappa |
| Actual slow mean-wave change | C=B+1/2 |

The seven linear/cross/fast rows tie Bnext=B+1/2-kappa. The square and slow
rows are stronger on the WHOLE half-line: subtracting Bnext gives affine
nonnegative slopes and nonnegative values at B0. Polynomial S factors remain
in the classes at a tie. This is not a lower-bound or saturation theorem.

Every source in (4),(5) keeps a wave factor, the original sqrt(zeta)P bound,
full smooth endpoint extension, and closed supported geometry. Source/edge
quotient estimates give its f/P jets; the fixed constrained inverse can be
applied to the NEXT regular source. No additive Gaussian tail has been
forwarded to that inverse. Combined with the parent's exact current mean
target recurrence, this proves the conditional finite update

    S(B) -> S(B+1/2-kappa),                                  (6)

with the EXACT excluded ledger (4). It does not independently prove the
listed initial/core/background/inverse estimates or one infinite field.

## 5. Quantitative Gaussian tail bounds for every fixed physical jet

The reviewed separation gives P<=exp(-c ell^2), c>0, where 1-psi is nonzero
on the source or a cutoff derivative occurs. For fixed finite stage, physical
jet order m and term in (4), the input dictionary bounds that tail by

    C q^-D ell^b exp(-c ell^2),    ell=log(1/q)>=1, D,b>=0.   (7)

Finite edge powers are first absorbed using the ORIGINAL half-exponents
t1^2/2 and 2, t1=C_setup^-20; their width-dependent constants are retained.
Physical phase/floor/coordinate derivatives only change D,b and constants.
No estimate uniform in stage, jet order, lambda^-1 or t1^-1 is asserted.

For any requested N>=0 and n=ceil(b/2), the exponential series and square
completion give the explicit whole-half-line inequality

    q^-D ell^b exp(-c ell^2)
      <= n!*(2/c)^n*exp((D+N)^2/(2c))*q^N.                  (8)

Indeed ell^b exp(-c ell^2/2)<=n!*(2/c)^n for ell>=1, while
-(c/2)ell^2+(D+N)ell <=(D+N)^2/(2c). This holds on the SAME accepted domain
for every fixed m,N; only its constant changes. Smooth physical flatness
does not mean zero residual or vanishing uncut exit data.

There is a DIFFERENT normalized-source requirement. A family
f_ell=epsilon^B ell^-a P_ell, P_ell=exp(-c ell^2), is physically flat by (8),
but claiming f_ell in W_(B+delta) would require

    exp(h*delta*ell) <= C ell^a,     delta>0,                (9)

which is false as ell tends to infinity. This scalar control is not a
counterexample to the actual construction. It prohibits absorbing a tail
into an improving normalized regular source solely because it is physically
flat. The older exact homogeneous-pulse response obstruction remains valid.

At fixed finite j,m, the regular normalization Q^(-2A-1/2), the accepted
physical derivative cost q^(-2m), Q comparable to q, and one q^-1 absorbing
the finite polynomial S yield

    |R_reg,new|_m <= C_jm q^(h Bnext-Km),
    Km=2A+3/2+2m,      A=1/2+h.                              (10)

The current mean rows/radial moment obstruction have stronger normalized
orders by the parent's fresh mean maps; their interior A0/A1 remainders fit
any fixed power. With the SEPARATELY assumed full base/exterior flat bounds,
(8),(10) give the same conditional bound for the complete finite residual.
The bare base stress divergence is canceled by the full current mean
equations, not incorrectly declared flat by itself. Constants may grow
arbitrarily with j,m: (10) is not fixed-q contraction or an arbitrary diagonal.

## 6. Checks, remaining limits and next obligation

The exact finite diagnostic independently differentiates the complete
cylindrical Fourier momentum operator, complete cut vector potentials and
pressures, and the full before/after nonlinear fields. It checks all terms
of (2),(4),(5), both cross orientations, new squares, pressure-only tails,
failure of div(psi*v)=0, mode-zero pressure independence, nonzero Pi2 axial
gradient, conjugate reality and generated harmonics. Its graph-free polynomial
fields have finite endpoint regularity and are NOT actual pulse solutions,
an evaluation of the source inverse, or C-infinity supremum certificates.
Separate time-only controls verify vanishing spatial curl/pressure commutators
for a temporal factor, and a spatial-peak control rejects Gaussian flatness
for a generic spatial localizer. These distinct controls must not be conflated.
Independent graph/mixed-viscosity and actual inverse diagnostics are replayed
from the relevant foundations. The continuum estimates are the conditional
arguments above, not a test count or an exponent scan.
The audit's recursive historical check verifies bound source bytes, not a
hidden re-execution of all past diagnostics. Fresh relevant evidence replays
are separate local/CI commands and are reported separately.

Next independently discharge the retained ACTUAL initial canonical seed,
core/background and common-domain inverse hypotheses, identifying exactly
which analytical estimates remain rather than promoting their protocol
labels. Then join this same complete residual to the distinct accepted-tail
single-field/force construction with its own controlled schedules. Historical
first-band/old-force acceptance and an unconditional theorem remain unproved.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json`.
Use fresh outputs; never overwrite historical evidence.
