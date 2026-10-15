# Conditional initialization and complete finite residual cycle

Status: `conditional-initial-state-and-complete-finite-cycle-order-gain`.

This checkpoint derives the initial corrected state and a **conditional
complete finite-cycle residual estimate**. It includes the nonzero angular
harmonics, exact tangential mean equations, radial pressure obstruction,
all three integral defects, source supports, and all old/new nonlinear
products. One cycle improves the declared normalized orders by 1/10.

The result is conditional on the inherited continuum leading/background
estimates, pulse/inverse estimates and auxiliary support geometry listed
below. It is not independent verification of those foundations, a numerical
evaluation of the actual global field, a convergent infinite construction,
or a blow-up proof. No unconditional complete-residual or forcing status is
promoted. The gain is in an asymptotic order, not a measured contraction of
the error at one fixed scale.

Parent: `73116848b1e53cad36e3c0e314aaaac452518a2e`.
The reference under audit is the supplied
[proposed manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(6.8)--(6.32), (7.39)--(7.42), and (9.1)--(9.10), Propositions 9.3--9.6.
Its conclusion is not assumed. PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
All earlier audit sources, evidence and GOAL.md remain unchanged.

## 1. Precisely what is assumed and what is proved here

Keep the same fixed N, h in (0,1/100), positive patch exponent lambda,
base, primary waves, discrete labels and inverse maps. Write epsilon=Q^h,
S=ell^2, k=ceil(epsilon^-1/2), and kappa=10^-5. We use the inherited classes
W_alpha, M_alpha and S_alpha, with every fixed coefficient derivative:
waves have weight sqrt(zeta)*delta^-d*P, means have zeta*delta^-d, and moments
have no radial weight. Polynomial powers of S, constants and inverse-distance
powers may depend on the finite stage and derivative order. They are uniform
in band, label, copy and point. For 0<epsilon<1, a larger exponent is a
stronger estimate.

The premises, with their provenance, are:

1. The [summed background](../background_summation/README.md) has the exact
   retained tangential stress, a flat augmented error, all fixed normalized
   derivatives, b=O(epsilon), and V,G=O(1) on the wave and mean supports.
   Its reserved mean patch has G=0 and the unchanged nonzero power-law V.
2. The actual leading profile has the inherited positivity, shear and cone
   gaps and flat edge estimates. These remain conditional, as in the
   [actual-band checkpoint](../actual_band_pulse_bounds/README.md).
3. The prescribed common-torus geometry has separated CLOSED enlarged
   supports for distinct interacting labels, smooth extensions, bounded
   covering differences and the regular slow mesh S^-3. We do not replace
   that hypothesis by a sampled rectangle plot. Its product rule is that
   derivatives of waves of different labels have zero products.
4. The [pulse checkpoint](../actual_band_pulse_bounds/README.md) supplies
   the fixed sourced inverse W_alpha to W_alpha, the homogeneous pulses,
   pressure order alpha+1/2, full linear remainder alpha+1/2-3kappa,
   Gaussian endpoint tails and the signed inverse at the FIXED positive
   primary amplitudes. Sources for a global linear solve must extend
   smoothly by zero at both pulse endpoints.
5. The exact curl construction, real-primary covariance estimate and signed
   actual-velocity covariance identity hold as derived in the
   [physical interface](../physical_wave_interface/README.md) and pulse
   checkpoint. Later complex corrections use the conservative signed
   covariance remainder, not the special real-primary cancellation.
6. The [compact mean maps](../compact_mean_corrections/README.md) have their
   proved conditional estimates, exact pressure/flux identities and fixed
   five-equation inverse, with all stated input regularity. The two velocity
   moments are preserved. Their slow update gains 0.9-2kappa on the three
   defects when the current total means have the declared bounds below.

The new work is the nonlinear interaction derivation, closure of the source
hypotheses under the full sequence, initialization, and the inequalities for
EVERY finite cycle. No cycle-gain hypothesis is on this list.

On a fixed common chart the graph rules are
D_r:C_alpha->C_(alpha-kappa), D_z=epsilon partial_Z:C_alpha->C_(alpha+1),
slow time -epsilon partial_T:C_alpha->C_(alpha+1), and fast time preserves
the order. Coefficient derivatives preserve the order. The explicit number
of graph derivatives in an operator determines its loss; requesting more
coefficient derivatives of an already estimated expression does not repeat
that graph loss.

## 2. The fast quadratic term uses the COMPLETE curl amplitude

Let w=a exp(ikm Phi) be one complex curl harmonic, m!=0, and n its normalized
spatial phase gradient. Its cylindrical coefficient divergence is

\[
\operatorname{div}_*a=D_ra_r+a_r/R+D_za_z.
\]

Exact incompressibility, before any approximation, gives

\[
\operatorname{div}_*a+ikm(n\cdot a)=0,
\qquad n\cdot a=\frac{i}{km}\operatorname{div}_*a.
\tag{1}
\]

Consequently a in W_alpha has n.a in W_(alpha+1/2-kappa). This statement
applies to a=t+r, including its curl correction; it is not the assertion
n.a=0 for the complete amplitude.

For two harmonics of the SAME label, the coefficient of transport of b by a
is exactly

\[
a_rD_rb+a_zD_zb+\frac1R(-a_\theta b_\theta,a_\theta b_r,0)
 +ikm'(a\cdot n)b.\tag{2}
\]

The last term seems to lose 1/2, but (1) supplies it back. The remaining
amplitude derivatives cost at most kappa and the frame terms cost none.
Thus the symmetric nonlinear interaction has order alpha+alpha'-kappa.
For a zero output harmonic it is in the mean class; otherwise it is a wave
with the same label and harmonic m+m'. Leibniz' rule proves the estimate
at each fixed coefficient derivative order. Products have weight zeta P^2;
P^2<=P and zeta is stronger than sqrt(zeta), so the output retains its wave
envelope when needed. The same conclusion holds for a wave's self-product.

For a mean m with tangential order mu and radial order mu+1, its transport
of a wave phase costs at most k=O(epsilon^-1/2). All other contributions
are better: wave transport of a mean loses at most kappa; radial mean
transport has the extra radial power; axial transport has D_z; cylindrical
connections cost none. Hence

\[
(w\cdot\nabla_*)m+(m\cdot\nabla_*)w
 \in W_{\alpha+\mu-1/2}.\tag{3}
\]

This half-power loss really occurs, for example when a constant axial mean
advects cos(kz)e_x. It must not be replaced by an amplitude derivative.
Also, (1) cannot cancel a.n' for a DIFFERENT phase. Two overlapping
divergence-free Cartesian plane waves cos(kx)e_y and cos(ky)e_x have cross
transport with RMS k/sqrt(2). Their same-label cancellation is unavailable.
This explains the specific role of the inherited closed-support separation.

The local checks independently differentiate complex Cartesian curl fields
and compare with (1)--(2), including negative and zero output harmonics.
They also check the complete same-field identity

\[
\begin{aligned}
R(u+d,p+\pi)-R(u,p)={}&\partial_td-\Delta d+\nabla\pi\\
 &+(u\cdot\nabla)d+(d\cdot\nabla)u+(d\cdot\nabla)d.
\end{aligned}\tag{4}
\]

The fixture uses full compact old/new curl waves and full compact old/new
mean curls. Pressure, viscosity, old/new cross products and the new square
are all nonzero. It tests the identity, not decay of the actual global flow.

## 3. Source closure survives a mean that fills the auxiliary torus

Here is an induction on the expressions in a finite construction.

* Actual waves are curls of potentials cut off in the fixed slow,
  transverse and temporal supports. Differentiation preserves closed
  support and smooth zero extension, although a pre-cutoff inverse solution
  need not vanish at the pulse exit.
* A mean operation creates no angular harmonic. It can aggregate sources
  across the radial shell or fill the auxiliary torus, while preserving its
  M bounds and smooth shell extension. A subsequent nonzero harmonic
  contains a wave factor, so multiplication and differentiation retain that
  wave's slow box, rectangle, temporal support and envelope P.
* Same-label wave products add the integer harmonic multipliers and obey
  (2). Different-label derivative products vanish by premise 3. In
  particular a zero harmonic is a genuine coefficient mean, with no
  unresolved difference of two distinct fast phases.
* Radial stress inversion uses the auxiliary average after subtracting its
  weighted moment bump. It creates an auxiliary-independent signed target.
  The fixed primary inverse divides by the original positive amplitude,
  with the inherited half-edge-weight estimate; it takes no square root of
  an updated stress. Applying its slow cutoff and exact curl restores the
  same wave supports and P envelope.
* Forward pulse propagation fixes the slow and transverse point along the
  prescribed common-torus path. A zero source on an entire parameter path
  has zero solution by uniqueness. Multiplication of its potential and
  pressure by the final temporal cutoff restores smooth endpoint zero
  extension. No divisibility by an old cutoff is used.

For clarity, a simple exact example shows why a mean cannot itself be
assumed rectangle-supported. If chi is a smooth supported torus bump,
N_t chi has zero Haar mean, but the zero-mean inverse is chi-<chi>, which
is nonzero outside the bump. Its product with a supported wave still has
the wave's support. The diagnostic exercises this example rather than
silently assigning the inverse its source's auxiliary support.

Each regular slow coordinate interval has O(S^3) grid cells in a fixed
enlarged chart, so the three-dimensional product mesh has O(S^9) boxes.
Two signs, boundedly many overlapping bands and bounded covering-index
differences add fixed factors. Radial mean aggregation therefore adds a
polynomial S cost, not an epsilon loss. It creates no new band because
q=q(z,t) is constant along the radial path. Actual coloring and separation
remain inherited premises; this count does not independently certify them.

Initially waves have m=+/-1 and their residual has |m|<=2. If a residual
has |m|<=M_j, a sourced update has those same harmonics, a signed update
adds only +/-1, and means add zero. A new quadratic residual has
M_(j+1)<=2 max(M_j,1). Thus every finite stage has a finite angular set
independent of band. Its size and harmonic constants may grow with stage.

Gaussian cutoff errors are flat in every fixed physical jet. Mean radial
cutoff errors are flat by repeated integration by parts with finitely more
input derivatives for any requested output order. Keep these errors in
the exact residual or its flat summand, and keep their actual velocities
in EVERY later product. A zero-mean flat coefficient can acquire a nonzero
mean when multiplied by another coefficient: eta cos(y)*cos(y) has mean
eta/2. Flat is not zero. Flat terms are reclassified or removed from the
flat ledger when an exact mean operation cancels them; they are not counted
twice. Nonzero Gaussian cutoff tails are excluded from subsequent forward
sources. The remaining sources have smooth temporal zero extension by the
support induction above, supplying the extra hypothesis of the saved
global linear remainder estimate.

## 4. Construct the initial state with the exact primary covariance

Let w0^tan be the assembled real transverse primary field, w its exact
curl velocity, and initially beta=v=gamma=0. Reconstruct pressure from the
FULL W_ab=<w_a w_b>_theta by the previous compact mean map.

The saved bounds give w in W_1/2, w-w0^tan in W_(1-kappa), full raw means
and defects of order H0=1-kappa, and linear/nonlinear nonzero residual
orders 1-3kappa and 1-kappa, respectively.

There is a better auxiliary-averaged tangential estimate. For these REAL
primary carriers the exact curl covariance has no mixed sine/cosine term.
Its normalized tensor error is M_(2-2kappa). The squared slow partition
recovers the leading stress exactly, including after differentiation; both
signs of one box count once. The summed base stress differs from that
leading stress by M_3 (epsilon times the O(epsilon^2) background correction).
Therefore, before any mean velocity update,

\[
\begin{aligned}
\bar E_\theta&=D_2(\bar W_{r\theta}-\Sigma_\theta)
                  +\epsilon\partial_Z\bar W_{\theta z},\\
\bar E_z&=D_1(\bar W_{rz}-\Sigma_z)
                  +\epsilon\partial_Z\overline{W_{zz}+p_m}
 \quad\in M_{2-3\kappa}.
\end{aligned}\tag{5}
\]

The exponent in (5) bounds BOTH rows. The axial pressure is included:
p_m has order H0, so its axial derivative has order 2-kappa. The row
identity is a difference of entire assembled tensors; derivatives of
partition functions are not dropped separately.

Apply the temporal mean update at H0, realize its axial part by the compact
potential, and recompute pressure and all defects. The increments have
tangential order H0 and radial order H0+1. After fast-time cancellation,
slow time, radial/axial fluxes, viscosity and changed pressure give remainder
order H0+1-2kappa=2-3kappa. Equation (3) gives wave error order H0.
Defects remain of order H0.

Apply the slow five-equation map to those CURRENT defects at H0 and
recompute pressure again. Its two zero-moment rows preserve the velocity
moments exactly, and its other rows remove the linear defect changes.
The exact nonlinear remainder from the preceding checkpoint gives defects
of order H0+0.9-2kappa=1.9-3kappa. Tangential and wave errors have the
same lower orders just established. Thus the constructed initial state
satisfies

\[
G_{wave}^{[0]}\in W_{1-3\kappa},\quad
(E_\theta,E_z)^{[0]}\in M_{2-3\kappa},\quad
(P,J_\theta,J_z)^{[0]}\in S_{1.9-3\kappa}.\tag{6}
\]

These are stronger than the convenient stage-zero bookkeeping B0=0.7,
C0=1.2. They are conditional bounds for the constructed functions, not
values obtained by evaluating the actual global leading profile.

## 5. A complete finite correction state and cycle

At stage j put B=0.7+j/10 and C=B+1/2. The state hypotheses are

\[
G_{wave}\in W_B,\quad (E_\theta,E_z)\in M_C,\quad
(P,J_\theta,J_z)\in S_C,
\tag{7}
\]
\[
w\in W_{1/2},\quad w-w0^{tan}\in W_{0.68},\quad
v,\gamma,p_m\in M_{0.9},\quad\beta\in M_{1.9}.
\tag{8}
\]

The velocity is the sum of the fixed base, actual curl waves, a mean
azimuthal-potential curl and a direct azimuthal mean. It is divergence-free
and real, and int R^2 bar(v)=int R bar(gamma)=0. Pressure is the linear
compact reconstruction from the current full radial source. In addition
to the nonzero harmonics and exact E rows, the radial residual is -rho P
plus its flat cutoff term. The base flat error is retained separately.
This specifies every component of the full momentum residual.

### Step 1: sourced wave curls

Solve the supported G_wave source with the fixed principal inverse, then
cut off the POTENTIAL and pressure and take the exact curl. The increment
is W_B, with pressure W_(B+1/2). The remaining nonzero orders, from the
full linear estimate and (2)--(4), are

| Contribution | Order |
|---|---:|
| full linear remainder, including viscosity | B+1/2-3kappa |
| cross with the complete old wave | B+1/2-kappa |
| new wave square | 2B-kappa |
| cross with the complete old mean | B+0.4 |

The covariance changes have order at least B+1/2=C, since 2B>=C.
Recomputing the mean equations, radial source and pressure gives all full
means and defects order C-kappa. The exact integrated identities give the
STRONGER weighted residual moments

\[
\int R^2\bar E_\theta=\epsilon\partial_ZJ_\theta,\qquad
\int R\bar E_z=\epsilon\partial_Z(J_z+c_\rho P)
 \quad\in S_{C+1-\kappa}.\tag{9}
\]

The derivative acts on the entire moving-bump product c_rho P. The moment
gain in (9) is essential for the next step.

### Step 2: signed averaged-stress update

Subtract normalized interior bumps with the moments (9), take the compact
weighted radial primitives of the auxiliary averages, and use their
negative as the covariance target Sigma. It is M_(C-kappa) and independent
of the auxiliary torus. The signed primary-amplitude update is W_(B-kappa).
Take its exact curl and associated pressure and recompute from the new
wave. Its nonzero errors have orders

| Contribution | Order |
|---|---:|
| linear remainder | B+1/2-4kappa |
| cross with total means | B+0.4-kappa |
| cross with old exact waves | B+1/2-2kappa |
| signed new wave square | 2B-3kappa |

For the auxiliary average, split the old wave into FIXED w0^tan plus its
W_0.68 remainder. The exact covariance change has its prescribed cross
entries Sigma, and the following remainders AFTER radial divergence:

| Contribution | Order |
|---|---:|
| old primary remainder times new transverse increment | C+0.18-2kappa |
| complete old wave times new curl correction | C+1/2-3kappa |
| square of signed new wave | 2B-3kappa |
| retained weighted-moment bumps | C+1-kappa |
| axial flux and changed radial pressure | C+1-2kappa |

In particular the square of the signed increment is retained. Since B>=0.7,
2B-3kappa=C+(B-1/2)-3kappa >= C+0.2-3kappa. Every row exceeds C+0.17
at kappa=10^-5. Meanwhile the full, unaveraged means and defects have
order H=C-2kappa, including the changes of radial source and pressure.
The conclusion after pressure reconstruction is

\[
\bar E_\theta,\bar E_z\in M_{C+0.17},\qquad
E_\theta,E_z\in M_H,\qquad(P,J_\theta,J_z)\in S_H.
\tag{10}
\]

### Step 3: temporal mean inverse

Use the CURRENT zero-auxiliary-mean E rows in the temporal inverse. Its
actual tangential/radial increments have orders H and H+1; retain the flat
axial realization remainder in the actual velocity. Fast time cancels the
chosen source. In the remaining tangential equations slow time gains one,
radial flux gains at least 1-kappa, axial flux gains one, and viscosity gains
1-2kappa. The pressure change has order H and its axial effect gains one.
Thus the new E rows have order

\[
\min(C+0.17,H+1-2\kappa)=C+0.17.\tag{11}
\]

Wave interactions with this update have order H by (3). The three defects
still have the conservative bound S_H. Exact equations use the realized
Delta gamma, not just its desired value before the cutoff correction.

### Step 4: current five-equation mean update

Apply the fixed slow inverse to the CURRENT defects at H and reconstruct
pressure again. Tangential errors are bounded by the same H+1-2kappa
ledger; no fast-time source is introduced. Wave changes again have order H.
All hypotheses of the saved nonlinear defect lemma hold: H>=1.19998>0.9,
the current total means obey (8), and the potential support remains in the
reserved patch. It gives

\[
(P,J_\theta,J_z)_{new}\in S_{H+0.9-2\kappa}
 =S_{C+0.9-4\kappa}.\tag{12}
\]

Combining ALL four steps, every wave row is at least B+0.1, every
tangential mean row at least C+0.1 and (12) exceeds C+0.1. The dominant
mean bound C+0.17 has 0.07 of slack above the requested gain. The auxiliary
covariance margin needed to establish that 0.17 is at least 0.00998.
`cycle_bounds.py` represents each row as an affine function of B and proves
these comparisons on the entire half-line B>=0.7 by the sign of its slope
and its value at 0.7. Sampling a finite list of stages is not that proof.

Each new wave correction has order at least B-kappa>=0.69999>0.68. Each
new mean/pressure increment has order at least H>=1.19998>0.9, with radial
order H+1>1.9. Finite sums therefore preserve (8). Reality follows from
conjugate Fourier pairs, real signed targets and uniqueness of the real
inverse equations. The two velocity moments remain exact. Section 3
provides the sources for the next cycle. This completes the conditional
finite induction, including the radial residual and flat part.

## 6. Common domain and fixed physical derivative losses

All inverses in the construction use the fixed background, phase normals,
positive PRIMARY amplitudes, mean bump matrices and supports. They are
linear for arbitrary finite source size. No updated cone or nonlinear
small-target solve is required. Higher harmonics add forward damping and
change constants, not the pulse domain. A single inherited qstar therefore
works for every finite stage and every fixed coefficient derivative; the
constants need not be uniform in stage or derivative order.

For use in a later full-field sum, a normalized order must survive fixed
Cartesian derivatives with a loss independent of stage. On the active
annulus R is uniformly positive and q is comparable to Q. Coefficient
graph derivatives and coordinate changes have fixed powers of Q and
polynomial S bounds. The phase formula is

\[
\Phi=p\theta+p_zZ/\epsilon+x_0R-v(pF+p_zG).
\tag{13}
\]

Every fixed derivative of its coefficient part is bounded by epsilon^-1
times a polynomial in S; all label choices are held fixed. A physical time
operator is -Q^-1 partial_T+Q^(-1-h)c_i N_i. Including k, a deliberately
coarse phase-derivative cost is at most Q^(-1-5h/2), which is better than
Q^-2 for h<1/100. Spatial derivatives, inverse-radius frame factors and
graph derivatives are also bounded by Q^-2 per fixed physical derivative.
Repeated chain/product rules give Q^(-2m) for a jet through order m.
Finite harmonic numbers enter constants only. The flat edge weight absorbs
every fixed inverse-distance power.

The residual normalization is Q^(-2A-1/2). Absorb a fixed-stage polynomial
in S by one further q^-1. A sufficient stage-independent jet loss is

\[
K_m=2A+3/2+2m,\qquad
|R(u^{[j]},p^{[j]})|_m\le C_{j,m}q^{hB_j-K_m},
\quad B_j=0.7+j/10,
\tag{14}
\]

with the flat part absorbed at any required finite order. Outside the
correction annulus only the inherited flat base error remains. At each
fixed derivative order, (14)'s exponent increases by h/10 per cycle.
Its initial exponent can still be NEGATIVE: for diagnostic h=1/200 and
m=0 the conservative exponent is -5013/2000. Thus the newly proved order
gain does not assert that a low-stage physical residual tends to zero,
or that its norm contracts at a fixed q. Constants C_j,m can grow without
any uniform control in j.
The quantifiers fix j before q approaches zero. Substituting a growing
j(q) into the constants in (14) would require a separate argument.

## 7. What this checkpoint verifies and what remains

The new analytic conclusion is conditional initialization, source closure
and full finite-cycle order gain, with a fixed physical jet loss. The
tests independently exercise the longitudinal identity, Cartesian
transport, exact residual increments, actual signed covariance remainder
orders, support failures, nonlocal mean inverses and all-stage rational
inequalities. They are not a numeric solve of the global constructed flow,
interval enclosures of its constants or independent foundation review.

```sh
python3 research/navier_stokes_cascade/results/residual_cycle_closure/test_audit.py -v
python3 research/navier_stokes_cascade/results/residual_cycle_closure/audit.py \
  --verify-record research/navier_stokes_cascade/results/residual_cycle_closure/evidence.json \
  --output /tmp/residual_cycle_closure_recheck.json
```

Next: sum the FULL wave and mean VECTOR POTENTIALS with shrinking q cutoffs,
including every curl-cutoff commutator, and compare with one fixed finite
stage using (14). The raw infinite correction series is not known to
converge. A full realization must preserve the exact divergence, support,
core velocity lower bound and finite energy; then localization and
initialization in physical space-time must produce admissible smooth
compact forcing through the singular time. These are separate open
obligations, alongside independent review of the inherited leading profile
and global auxiliary geometry. No blow-up theorem is claimed here.
