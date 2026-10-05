# Repair neighborhood and conditional frequency acceptance

Status: `repair-state-bound-and-conditional-frequency-acceptance`.

This checkpoint derives the missing correction-to-state estimate on the
actual first reserved patch, then applies the finite-frequency criterion.
It depends on the earlier analytic profile, compact envelope, loop and
modulation estimates. It does not independently prove those dependencies.

The proposed **Ccorr=A^128** and coefficient radius **r=A^-128** suffice.
Together with the inherited estimates, **N=1+floor(H^32)** satisfies the
leading-profile acceptance inequalities. This is an analytic selection of a
finite integer, not a numerically resolved flow. The resulting leading cone
argument is conditional on the inherited continuum estimates and edge
arguments. The full PDE corrections, smooth forcing, actual divergence and
independent mathematical review remain unverified.

Parent: `89e00da7d41e8ee0dfb77548e7c299308e73eb55`.
The construction uses the formulas already transcribed and checked in
[conditional stress realization](../stress_realization_audit/README.md),
[compact envelope](../compact_jet_envelope/README.md) and
[loop/modulation bounds](../loop_modulation_bounds/README.md).
The supplied manuscript source hash remains
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.

## 1. Norms and unchanged patch fields

Use A=C^4096>=2^16, G=exp(A^32), H=exp(G^4), as before.
On the patch X=rhop exp(y), 0<y<5,

\[
U_0=0,\qquad E_0=e(\eta)\exp[(-1/2-\lambda)y],\qquad0<\lambda\le.01.
\]

The modulation has ended before this patch. Its field values here equal the
original field, although its cumulative moments have an incoming discrepancy.
The same is true of the fields on the gaps between bump supports.

The coefficient norm used by the existing contraction is

\[
s=\max_j\max\{\|c_j\|_\infty,\|c_j'\|_\infty\}.
\]

For field and moment products use the angular algebra norm
||f||1=||f||infinity+||f_eta||infinity. In particular ||c_j||1<=2s. Do not
confuse these norms. The earlier envelope bounds ||e||1, field mixed jets,
original moment jets, and the required positive reciprocals by A.

Let beta_j(y)=sigma'((y-yj)/.1+.5)/.1. The earlier flat-step estimates give

\[
|\beta_j|\le640,\qquad |\partial_y\beta_j|\le2\,10^6.
\]

Their five supports are disjoint and independent of eta. At each radius at
most one bump contributes, including to its angular derivative. For the
actual corrections deltaU=e sum_U c_j beta_j and deltaE=e sum_E c_j beta_j,
this proves

\[
\|\delta E\|_1,\|\delta U\|_1\le d=1280As,\quad
\|\partial_y\delta E\|_1,\|\partial_y\delta U\|_1
\le q=4\,10^6As.
\]

The factor 2 in both estimates accounts for the coefficient norm. The bump's
radial derivative contains the inverse width squared, not just inverse width.

## 2. Positive-field neighborhood and incoming memory

Take r=A^-128 and s<=r. Since E0^-1<=A on this compact region,

\[
|\delta E|/E_0\le1280A^2s<1/2,\qquad d\le1.
\]

Consequently E0+deltaE>0, its inverse is at most 2A, and its inverse Hx,
where Hx=sqrt(2X)E, is at most 2A. The field angular norms remain below 2A.

The prior modulation estimate gives each moment/Pi an angular error
m0=32A^4H/N. At the selected N>H^32,

\[
m_0<32A^4H^{-31}<1/(16A^2),
\]

because 512A^6<H^31. It follows that the incoming moment angular norms are
below 2A. The same axis pressure datum is retained. Using the original
(4.16) gives |W_in|<=A+2Am0<2A and |Ns_in|<=A+12A^2m0<2A.
The four-term numerator B_in in Qs has magnitude at most 8A.

These are comparison bounds for the already modulated state, not an
assumption that its cumulative moments still equal those of the original.

## 3. Correction-to-state estimate

The exact differences of the shears give, using original |E_y|<=A and U0=0,

\[
|\Delta a|\le4Aq+4A^3d,\qquad |\Delta b|\le4Aq.
\]

The common expression is
(16*10^6 A^2+5120 A^4)s<A^6s for A>=2^16.
No smallness of radial derivatives is inferred from field values alone.

Integrate the five original moment densities through any part of the patch.
The support length is below A, sqrt(2X)<=2A and X^-1<=A. With d<=1 and
input field norms <=A, the same product estimates as the modulation note
give each moment and Pi a repair error bounded by

\[
m=8A^3d=10240A^4s
\]

in the angular algebra norm. The U-square term in S is retained; the U/E
cross term in J vanishes because the actual supports are disjoint.

Now subtract the ORIGINAL pressure-coordinate formulas, starting from the
incoming modulated state. Let B denote the numerator in Qs. Then
|deltaB|<=4m and |B_in|<=8A. Positive denominators give

\[
|\Delta W|\le2Am,\quad
|\Delta H_x^{-1}|\le4A^3d,\quad |\Delta E^{-1}|\le2A^2d,
\]
\[
|\Delta Q_s|\le10A^2m+32A^5d,\quad
|\Delta N_s|\le12A^2m+2Ad.
\]

The pressure and angular pressure terms in Ns contribute at most 4m. The
axis datum does not change. With |Ns_in|<=2A, X/L<=A^2, this yields

\[
|\Delta p_1|\le10A^4m+32A^7d\le112A^7d,
\quad |\Delta p_2|\le24A^5m+8A^5d\le200A^8d.
\]

Thus for the maximum norm of the FOUR state coordinates,

\[
\boxed{\|\Delta(a,b,p_1,p_2)\|_\infty
\le327680A^9s<A^{11}s<A^{128}s.}
\]

This is a value estimate, not an assertion of a small angular derivative of
the pressure coordinates. Its input coefficients have a C1 norm because the
pressure formulas use angular derivatives of moments. All width, denominator
and incoming-memory costs are included. This proves the proposed Ccorr on
the specified neighborhood.

## 4. Loop gaps and a permissible tolerance

The input compact estimate supplies a>=1/A, the original third/fourth gaps
>=1/A wherever vs>=2, and protected boundary margins >=1/A. The loop uses
d0=1/(16A), delta=exp(-A^20), and excursion |t|<=T<G. As proved in the
preceding loop notes, delta<=gamma/2 and
delta<=gamma^2/[4(1+J^2)], with gamma=15/(16A).

Where the cutoff is active, the four gaps have lower bounds

\[
2/(1+T^2),\quad\delta/8,\quad3\gamma/4,\quad\gamma^2.
\]

The last follows from 2(3gamma/4)^2-delta J^2/2>=gamma^2.
Where the cutoff is zero, the original shear is retained, vs-2>=delta/4,
and the other original gaps apply. Therefore a common gap floor is
kappa=G^-4 and a common first-coordinate floor is amin=G^-3. The loop
coordinates have magnitude at most 2A, including the unchanged p coordinates.

For the original four-coordinate cone Lipschitz estimate take
R=max(2,2A+1,2G^3)<=4G^3. Then the proposed epsilon=G^-64 obeys

\[
\epsilon\le1,\quad \epsilon\le G^{-3}/2,
\quad\epsilon\le G^{-4}/[512(4G^3)^{10}].
\]

The last comparison reduces to G^30>=512*4^10. The original 256R^10
gradient bound then leaves at least half the common gap under a state
perturbation of at most epsilon. This argument concerns cone coordinates;
it does not assume the physical Navier–Stokes residual is small.

## 5. Exact repair and selected finite N

The conditional repair argument supplies the exact five-row inverse bound
beta=1000 and C1 quadratic-map bound qstar=10^5. Those bounds concern exact
unit-mass bump integrals. Their numerical quadrature matrices are diagnostics.
The normalization of the actual incoming target is bounded by D/N with
D=H^16. It retains derivatives of e and the inverse-lambda loss.

Set Cstate=D=H^16, Ccorr=A^128, r=A^-128, epsilon=G^-64, and

\[
\boxed{N=1+\lfloor H^{32}\rfloor>H^{32}.}
\]

All scales are fixed finite positive real numbers, so this defines one finite
integer. To check the acceptance criterion without representing these scales
in floating point, compare logarithms. The state-error threshold satisfies

\[
\log\frac{2(Cstate+2Ccorr\,\beta D)}{\epsilon}
\le16G^4+\log4002+128\log A+64A^{32}
<16G^4+193G<32G^4.
\]

The contraction threshold has logarithm 16G^4+log(16*10^11), and the
coefficient-radius threshold has logarithm 16G^4+log4000+128log A.
Both are below 32G^4. Also N>=4AH, and the memory condition of section 2
holds. The comparisons follow for all A>=2^16 from G>A, log A<=A and
the displayed polynomial inequalities; endpoint checks accompany them.

Hence the exact C1 small root has coefficient norm <=2beta D/N<r/2,
its contraction has Lipschitz constant less than one half, and the combined
modulation plus repair state error is below epsilon/2. All four cone gaps
remain positive on the selected compact region.

At the right end, all FIVE moment functions agree exactly with the original
profile as functions of eta. Fields there already agree, and the axis datum
is the same. Thus the original integrated pressure/stress coordinates agree
at every later radius. Both protected collars and the analytic-axis rectangle
are untouched; the remaining correction patches are still available. The
edge arguments outside this compact region are inherited, not newly proved.

## 6. Scientific scope and next work

This completes a conditional leading-profile admissibility argument. It does
not construct or sum the physical waves that cancel the residual, establish
their derivative bounds, prove smooth force extension, or show a singularity
of the complete field. Independent review of the foundation is still needed.

The next work package is the full physical residual and correction budget.
Fix the selected N before choosing later physical scales. Keep its large
radial derivative constants in every estimate. Derive the divergence-free
wave/correction fields and their full residual, then check that the claimed
improvement survives pressure, viscosity, cutoffs and all cross terms.

The regression fixtures below are explicitly manufactured. They exercise
the actual bump formulas, five normalized moment changes, pressure/angular
dependencies, positivity failure, and original source ODE identities.
They are neither measurements of the actual profile nor proof of a blow-up.
Their amplitudes make missing terms observable; they do not lie in the actual
tiny coefficient ball or supply its derivative suprema.

## Reproduce

With mpmath 1.3.0, from the repository root:

```sh
python research/navier_stokes_cascade/results/repair_state_bounds/test_audit.py -v
python research/navier_stokes_cascade/results/repair_state_bounds/audit.py \
  --verify-record research/navier_stokes_cascade/results/repair_state_bounds/evidence.json \
  --output repair_state_reproduced.json
```
