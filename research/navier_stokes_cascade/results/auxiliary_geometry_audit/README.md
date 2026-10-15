# Constructive global auxiliary supports and common-chart compatibility

This batch replaces the support-separation black box used by the pulse
and residual checkpoints with an explicit construction. For all bands
ell>=16, with any further fixed lower threshold required by the profile,
it supplies:

* A proper **1,250-color rule** for the doubled slow-box interaction graph.
* Band-torus centers c_n=(0,n/2003), 1<=n<=1250, and any fixed radius
  \[
  0<r_0\le r_{\rm geo}=1/6057072.
  \]
  Closed enlarged absolute preimages of interacting labels are disjoint.
* A fixed common-covering bound Delta<=3, at most 2744 common-torus copies,
  preserved Haar averages and compatible physical derivative operators.
* Slow-partition coverage, pointwise bounded overlap, polynomial label
  counts under radial integration, and compatible lifted pulse paths.

These are direct geometric arguments with exact integer/rational
certificates and separate finite diagnostics. They do not construct or
independently verify the leading profile, actual pulse amplitudes or full
nonlinear correction sequence. The downstream exclusion implication remains
conditional on those analytic inputs. No blow-up theorem is asserted.

Parent: `464244c7334eb0a73579431e5155f1699565f61a`. The source under audit
is the proposed [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
Section 6, especially (6.1)--(6.22), Lemmas 6.1--6.3. PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Earlier sources, evidence and GOAL.md are frozen and unchanged.

## 1. Verify the source formula before auditing its consequences

The displayed equation (6.2), checked visually on PDF page 63, is

\[
J=\begin{pmatrix}3&1\\1&5\end{pmatrix},\quad
\Lambda=4-\sqrt2,\quad T_g=4+\sqrt2,\quad b=\sqrt2-1,
\]
\[
v_r=(1,-b),\quad v_t=(b,1),\quad
\rho=\frac{\log\Lambda}{\log T_g},\quad
\kappa=10^{-5},\qquad d_r=2\big((1+h)\rho-h\kappa\big).
\tag{1}
\]

An initial text-extraction reading lost the OUTER parentheses and suggested
a factor-of-two mismatch. Visual inspection rejects that suspicion: the
outer 2 multiplies BOTH terms. The frozen `actual_band_pulse_bounds/charts.py`
uses the correct displayed formula. **No source/repository normalization
defect was found, and no repair of that code is made.** The misparsed
variant is retained only as a deliberate failure control in this batch.
The visual check is a recorded source inspection, not a PDF render carried
out by the CI tests.

Fix the inherited 0<h<1/100. Let Q=2^-ell, epsilon=Q^h, S=ell^2 and

\[
i(\ell)=\left\lfloor
 \frac{(1+h)\ell\log2-2\log\ell}{\log T_g}\right\rfloor,
\quad c_i=T_g^iQ^{1+h},\quad M_i=\Lambda^iQ^{d_r/2}.
\tag{2}
\]

If theta is the fractional part of the expression in the floor, then

\[
c_i=T_g^{-\theta}S^{-1},\qquad
M_i=\Lambda^{-\theta}\epsilon^{-\kappa}S^{-\rho},\qquad 0\le\theta<1.
\tag{3}
\]

Thus T_g^-1/S<c_i<=1/S and
Lambda^-1 epsilon^-kappa S^-rho<M_i<=epsilon^-kappa S^-rho.
In particular M_i^-1<=Lambda epsilon^kappa S^rho, exactly the gain used
by the radial Fourier integration-by-parts argument. The proof is algebraic
on the ENTIRE band half-line, not inferred from samples. Log-domain
diagnostics include ell=10^12 without materializing Q, carrier integers
or the huge matrix J^i.

The eigenvector identities J v_r=Lambda v_r, J v_t=T_g v_t are immediate
from b^2+2b-1=0. The exponent d_r is positive: for example Lambda>2,
T_g<6, log2>1/2 and log6<2 give rho>1/4, which dominates h kappa.

## 2. Uniform covering-index difference

The real expression defining i has derivative
((1+h)log2-2/ell)/log T_g. For ell,ell'>=16 with |ell-ell'|<=4,
the mean value theorem and the cost of flooring give

\[
|i(\ell)-i(\ell')|
\le1+\frac{4(1+h)\log2+8/16}{\log T_g}.
\tag{4}
\]

The elementary bounds log2<7/10 and log T_g>5/3 yield the rational
upper bound 1873/625<3. `exponential_upper` proves these logarithm bounds
using rational Taylor sums with a geometric tail, rather than rounding a
floating logarithm. We choose the convenient, non-sharp Delta_max=3.
This bound remains valid if the lower band threshold is increased later.
The same elementary bounds show i(ell)>=0 on ell>=16.

## 3. Slow partitions, closed boxes and an explicit coloring

Let beta(s)=exp(-1/(1-s^2)) for |s|<1 and zero otherwise. Normalize its
integer translates by the square root of their squared sum. Every real
s has a nearest integer with |s-a|<=1/2, so the denominator is at least
exp(-4/3). The normalized translates are smooth, flat at their support
boundaries and have squared sum exactly one. Their fixed derivatives
are uniformly bounded by periodicity and compactness.

Use this construction on log2(q) for the dyadic partition chi_ell, and
as a tensor product on (R,Z,T)/delta_ell for the slow mesh, where
delta_ell=S^-3=ell^-6. Mesh box B_(ell,a) is centered at delta_ell a,
a in Z3, with closed support halfwidth delta_ell. Its doubled box has
halfwidth 2 delta_ell. On the active shell the selected indices are those
whose closed boxes meet the active chart set, including tau=0 BEFORE
representatives and other discrete data are fixed. This set is finite
at every fixed band: with s=q/Q in [1/2,2],

\[
\sqrt{X_a}\le R\le2\sqrt{X_b},\qquad
|Z|\le2^D,\qquad 0\le T\le2.\tag{5}
\]

If needed increase ell0 so delta_ell<=sqrt(Xa)/6. Then the selected
doubled boxes stay a positive distance from R=0. The local extension at
tau=0 is the inherited smooth background extension; no new boundary
labels are chosen as t approaches one. All discrete band, mesh and
frequency choices are frozen during differentiation.

For eta_(ell,a,sigma)=chi_ell chi_(ell,a), sigma in {0,1}, the two signs
duplicate the same slow cutoff. The squared partition is counted ONCE
per box, not twice per sign:

\[
\sum_{\ell,a}\eta_{\ell,a,0}^2=1,\qquad
\eta_{\ell,a,0}=\eta_{\ell,a,1}.\tag{6}
\]

Join labels whose doubled slow boxes intersect in physical coordinates
and whose band difference is at most four. This graph contains every
pair with overlapping original slow supports: their dyadic supports
already imply |ell-ell'|<=2. Define the color

\[
\nu(\ell,a,\sigma)=1+\sigma+
2\{[a_1]_5+5([a_2]_5+5([a_3]_5+5[\ell]_5))\},
\tag{7}
\]

where residues are in {0,1,2,3,4}. There are exactly 5*2*5^3=1250 colors.
For bands differing by one through four, their band residues differ,
regardless of the spatial indices. At the SAME band, intersecting doubled
boxes have |a_j-a'_j|<=4 in every coordinate. Equal residues modulo five
then force a=a'; opposite signs still have different colors. Thus (7)
is a proper coloring for ALL indices, including negative mesh indices.
No countable greedy selection or unspecified generic tuple is needed.

A period-four spatial rule would fail for the chosen CLOSED doubled-box
graph: centers four mesh steps apart have a common boundary point.
Period five is used for this particular graph, not claimed necessary
for every possible localization method.

## 4. Explicit rational centers and closed support separation

Choose the prime P=2003, verified by exact trial division through 44,
and c_n=(0,n/P), 1<=n<=1250. The relevant matrix powers are

\[
J^2=\begin{pmatrix}10&8\\8&26\end{pmatrix},\qquad
J^3=\begin{pmatrix}38&50\\50&138\end{pmatrix}.\tag{8}
\]

For Delta=1,2,3 the first coordinate of J^Delta c_n is respectively
n/P, 8n/P, 50n/P modulo one. Each is nonzero because P is prime,
P>1250 and P>50. Consequently for ALL colors m,n, even m=n,
c_m != J^Delta c_n mod Z2. At Delta=0 distinct colors have distinct
centers. Every relevant forbidden center difference therefore has
torus sup-distance at least 1/P from the integer lattice.

The closed enlarged band rectangle is c_n+xi v_r+eta v_t with
|xi|,|eta|<=2r0. Its displacement e from c_n has ||e||infinity<=4r0,
since ||v_r||infinity=||v_t||infinity=1. Also

\[
\max_{0\le\Delta\le3}\|J^\Delta\|_\infty=188.\tag{9}
\]

If two enlarged ABSOLUTE preimages at levels i and i+Delta intersect,
their centers must satisfy
c_m-J^Delta c_n=J^Delta e_n-e_m mod Z2. For r0<=1/6057072,

\[
\|J^\Delta e_n-e_m\|_\infty
 \le4(188+1)r_0\le\frac1{4P}<\frac1P.\tag{10}
\]

This contradicts the center-congruence gap. Interacting same-level labels
have distinct colors by (7); different levels satisfy the stronger
constraint for every color pair. Thus all required CLOSED enlarged
absolute preimages are disjoint. The congruence test has remaining margin
3/(4P). This is NOT a uniform Euclidean separation distance between
absolute-torus components: their sizes and distances may shrink as i grows.

The same radius makes each band rectangle injectively parametrized: a
difference of two enlarged displacements has sup norm <=8r0<1, so it
cannot be a nonzero integer vector; the eigenvector basis is invertible.
Any further fixed decrease of r0 preserves the proof. Pulse and background
small-q thresholds may then be chosen with these fixed constants BEFORE
the correction iteration, as in the earlier local estimates.

The earlier r0=.04 in the manufactured band fixture is NOT this global
support radius. No actual global radius was numerically fixed by that
fixture. The next pulse/source review must use the explicit geometry above
and its associated fixed constants, allowing the common q threshold to
shrink before iteration. No radial frequency N or h is silently changed.

Pairwise distinct band centers alone would not suffice: c and Jc may be
distinct centers while their level-i and level-(i+1) absolute preimages
share a point. This failed shortcut is retained in the diagnostics.

## 5. Common tori, lift counts and physical derivative operators

Fix (z,t) with q>0 and a sufficiently small neighborhood U whose possible
dyadic indices span at most four. Set i0=min i(ell) over this BAND set,
not merely labels active at a particular radial point. With H=J^i0 Y,

\[
Y_i=J^\Delta H,\quad 0\le\Delta=i-i_0\le3,\qquad
R_\gamma^H=(J^\Delta)^{-1}(R_\gamma).\tag{11}
\]

The covering has degree det(J)^Delta=14^Delta, at most 2744. The exact
integer kernel enumeration in `kernel_points` constructs its points as
adj(J^Delta)k/det(J^Delta) modulo Z2 and checks all images and counts.
Small injective rectangles lift to that many disjoint components. The
change-of-representation matrices have norm at most 188. Haar averages
are preserved: a nonzero Fourier character n pulls back to nonzero
(J^Delta)^T n, so only the zero character contributes to either average.

The absolute physical evaluation is

\[
Y(r,t)=v_r r^{d_r}+v_t t\pmod{\mathbb Z^2}.
\tag{12}
\]

It is smooth for r>0. Auxiliary-dependent corrections vanish near r=0
at each t<1; their smooth zero extensions, not smoothness of a fractional
power at the axis, supply the Cartesian field there. Chain rule and (1)
give exactly

\[
D_r=\partial_R+M_i d_rR^{d_r-1}L_i,\quad
D_z=\epsilon\partial_Z,\quad
t_*= -\epsilon\partial_T+c_i N_i,
\tag{13}
\]

with L_i=v_r.grad_(Y_i), N_i=v_t.grad_(Y_i). In the common torus,
M_i0=Lambda^-Delta M_i and c_i0=T_g^-Delta c_i. The pullback derivatives
contribute the inverse factors Lambda^Delta,T_g^Delta, so the physical
operators are unchanged. Slow charts and physical normalizations are
also converted when Q changes; their factors are bounded for a band
difference at most four. Comparisons of representatives use one common
physical normalization, not an equality of differently normalized arrays.

## 6. Polynomial counts and globally compatible operations

At a point at most three closed dyadic supports meet, at most three
one-dimensional mesh supports meet per coordinate, and there are two
signs: at most 162 slow labels. After auxiliary localization, the
separation proof gives at most one actual wave label at a point of the
extended domain or its physical evaluation. Derivatives of smooth zero
extensions stay in their closed supports, so all cross-label derivative
products vanish as required by the finite-cycle argument.

Per band, a compact chart box with side lengths L_j meets at most
2 product_j(ceil(L_j S^3)+3)=O(S^9) labels. A bounded radial line at
fixed Z,T meets at most 18(ceil(L_R S^3)+3)=O(S^3). Fixed-factor
enlargements and at most 2744 common-torus lifts change constants only.
In particular radial integration may collect many mesh boxes, but its
loss is polynomial in S, not a new power of q. Since q(z,t) is independent
of r, it NEVER introduces a new band or a larger covering-index gap.

The underlying graph also has uniformly bounded degree: for band changes
at most four and ell>=16, chart scaling is at most 16 and the mesh ratio
at most (5/4)^6. A doubled neighboring interval has center-index length
at most 4(1+16(5/4)^6)<249, so 256 candidate indices per coordinate
suffice. This gives the coarse degree-plus-one bound 18*256^3; coloring
(7) is much smaller and does not need that bound.

Compatibility is proved in the absolute variable, not imposed separately
on each common chart. A primary coefficient localized in an injective
band rectangle has a unique local coordinate formula, cut off smoothly
before its boundary. It defines one periodic band function and hence one
absolute pullback. On overlapping U,V, pull their representatives back
to the same absolute Y. Both equal the same sum; any band absent from
one neighborhood has identically zero slow cutoff on the overlap.
Smooth zero extension identifies all derivative jets at support changes.

If a band is absent near (z,t), its contribution is zero along the WHOLE
radial integral and auxiliary pulse path: both keep (z,t), hence q, fixed.
For the pulse inverse, zero initial data then give zero contribution by
uniqueness. Thus changing the common index cannot leave hidden coarse-band
memory that would fail to descend to the new common torus.

Products, sums and derivatives of these physical pullbacks remain
compatible. So do Haar averages and torus translations. Directional
Fourier inverses preserve the pullback frequency lattice; for example
N_abs^-1 Pi_i=T_g^-i Pi_i N_i^-1 on zero-mean functions. Radial primitives
use the SAME physical operator

\[
I g(r,Y)=\int g(r',Y+((r')^{d_r}-r^{d_r})v_r)\,dr'.\tag{14}
\]

Deck translations commute with this shift, and q/band membership is fixed
throughout the integral. Compact pressure and mean-velocity repairs must
use the previously specified physical q-based bumps, such as
chi_m(r^2/(2q)); independently chosen fixed-R chart bumps would not glue.
These facts prove preservation of descent and overlap identities for
the prescribed operations, not for arbitrary unrelated local formulas.

## 7. Lifted pulse paths, including mixed-band sources

Let lambda_t(w)=v_t.w/(1+b^2). On one common-torus lift set
eta_g(H)=lambda_t(J^Delta H-c_gamma-k_copy),
v(H)=(eta_g(H)+r0)/c_i, and L_s=2r0/c_i. The path to pulse time w is

\[
H_w=H+T_g^{-\Delta}c_i(w-v(H))v_t,\qquad 0\le w\le L_s.
\tag{15}
\]

Multiplication by J^Delta shows its band point has the SAME transverse
coordinate and eta_g=c_i w-r0. It remains in that lifted rectangle.
Moreover

\[
D_Hv=T_g^\Delta\lambda_t/c_i=O(S),\quad
D_H H_w=I-v_t\lambda_t,\quad D_H^k H_w=0\ (k\ge2),
\tag{16}
\]

and D_r v=D_z v=0, t_*v=1. In absolute coordinates the same path is
Y'=Y+T_g^-i(eta'_g-eta_g)v_t, independent of the chosen common chart.
A deck translation preserving H translates the path by the same amount
and merely reindexes the copy. Uniqueness of the linear zero-initial-data
ODE therefore preserves the pullback identity under chart changes.

A source built from several bands may not descend to Y_i. It must be
evaluated on EACH common-torus lift along (15), not averaged over other
preimages. The argument proves compatibility of that liftwise inverse;
its pulse growth and derivative estimates remain the earlier analytic
estimates to be reviewed with the inherited background. Transverse
cutoff and the time cutoff strictly inside the rectangle then give the
required globally smooth wave extensions.

## 8. Verification scope and next research step

Exact certificates cover the residue coloring, all covering-center
constraints, prime, radius margin, determinants and kernel copies.
High-precision finite diagnostics check squared partitions, floor scales,
actual phase-map chain rules and every lifted pulse-path identity.
They do not replace the all-band algebraic proofs above. Controls reject
mere center distinctness, an excessive radius, a period-four closed-box
coloring, omitted graph derivatives and the misparsed PDF formula.

```sh
python3 research/navier_stokes_cascade/results/auxiliary_geometry_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/auxiliary_geometry_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/auxiliary_geometry_audit/evidence.json \
  --output /tmp/auxiliary_geometry_audit_recheck.json
```

The pure support/covering geometry is now explicit. The next batch should
trace the primary pulse phase/frame and source estimates with these fixed
geometric constants, then critically audit the inherited leading profile,
common background collar, mean inverses and parameter inequalities. None
of those analytic inputs is independently re-proved by this geometry
certificate. The conditional exclusion chain is not promoted to a verified
blow-up construction, and no obstruction to the construction is claimed
from the rejected text-extraction false alarm.
