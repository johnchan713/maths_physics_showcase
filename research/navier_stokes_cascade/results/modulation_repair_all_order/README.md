# Fixed-frequency modulation and the SAME exact repair at every fixed order

Parent: `9f4ae8973fdd1384530085505010597033e206cf`.
This checkpoint proves an all-order MAP estimate for compact modulation,
the moving phase inverse and the accepted five-row repair. The selected
N=1+floor(H^32), all bump supports, lambda>0 and pressure datum are unchanged.
It avoids an unjustified simultaneous contraction in infinitely many norms.
Historical sources, evidence and tolerances are preserved.

The input is a smooth PRE-MODULATION profile on the WHOLE compact rectangle
from modulation entry through repair, with the positive floors and finite-order
acceptance estimates already stated in
[loop bounds](../loop_modulation_bounds/README.md),
[repair bounds](../repair_state_bounds/README.md) and
[higher jets](../higher_jet_budget/README.md). The parent independently supplies
all orders through the inner continuation/shape part, NOT the entire exterior
prefix. This theorem does not silently supply the missing annular/exterior
input. It does not prove unconditional actual-pulse bounds, complete residual
improvement or blow-up. Tests are manufactured equation checks, not a
continuum proof or evaluation of the enormous actual field.

Defining source: (C.5), (C.9)--(C.16) of the retained
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The exact five-row map is transcribed in
[stress realization](../stress_realization_audit/README.md), section 5.
Its final theorem is not a premise.

## 1. Smooth variance root INCLUDING zero variance and zero pressure

Write q(z,theta)=exp(z sin theta)/I0(z), z=mu p,
R(z,theta)=int_0^1 q_z(tz,theta)dt and R2(z)=<R(z,theta)^2>.
The REMOVABLE formula is

\[
 t=t_s+d_0\mu R(\mu p,\theta),\quad
 V=d_0^2\mu^2R2(\mu p),\quad
 F(\mu,p)=d_0\mu\sqrt{R2(\mu p)}.                 \tag{1}
\]

I0 is positive on the real line. Thus R2 is smooth and positive: for
z!=0 it is <(q(z)-1)^2>/z^2>0, while R2(0)=<sin^2 theta>=1/2.
On the frozen compact |mu p|<=G every fixed derivative of R2 and its
positive square root has a finite bound. F is the SIGNED root, not sqrt(V)
on both sides of zero. F_mu(0,p)=d0/sqrt(2); F_p(0,p)=0.
The tilted-variance argument gives F_mu>0 away from zero. On the compact
root family it therefore has a positive minimum m_F. The previous bound
F_mu>=(d0/32)exp(-4|mu p|) gives an explicit valid smaller floor.

Set s=zeta(vs)*sqrt((2+delta_L/2-vs)/a), extending it by zero outside
the flat cutoff support. Where zeta is nonzero the radicand numerator is
at least delta_L/4, and a has a positive floor. Hence s is smooth across
both cutoff boundaries. The equation F(mu,p)=s has the SAME unique root
as the accepted nonnegative variance equation. Near s=0 use its signed
extension to apply the ordinary implicit function theorem; no division by
V_mu, p, mu or s is allowed there.

Here is an explicit all-order recursion. At any base point write normalized
Taylor jets with h=(Delta y,Delta eta), mu=mu0+sum mu_alpha h^alpha.
After coefficients of degree <n are known, let mu_<n be that polynomial.
For every |alpha|=n,

\[
 mu_\alpha=-{[h^\alpha]\{F(mu_{<n},p(h))-s(h)\}\over F_\mu(mu0,p0)}.
                                                               \tag{2}
\]

Only the highest coefficient is missing, so its multiplier is exactly
F_mu at the base point. Multivariate Taylor composition gives finitely
many lower-order terms, each bounded using input jets and m_F^-1.
This proves finite mixed jets of every fixed order by induction. Compactness
and the same positive floor make those bounds uniform over the rectangle.
The recurrence concerns jets of the smooth solution furnished by IFT,
not convergence of a radial Taylor series across flat cutoffs.

## 2. A moving phase, not a frozen auxiliary angle

Use the exact mean-preserving lift

\[
 w={a(1+t^2)\over2\pi v},\quad
 \phi(x,\theta)=\int_0^\theta w(x,r)dr,\quad
 \phi(x,\Theta(x,\varphi))=\varphi.                 \tag{3}
\]

The variance equation gives int_0^(2pi)w=1 exactly; w has the retained
floor G^-3. Apply the same coefficient recurrence (2), replacing F by
phi and the denominator by w, now with h=(Delta y,Delta eta,Delta varphi).
This supplies ALL mixed derivatives of Theta, including terms where
derivatives hit w. In particular Theta_x=-phi_x/w and Theta_varphi=1/w.
The composition t(x,Theta(x,varphi)) needs both pieces
t_x+t_theta Theta_x; fixing Theta misses real terms at higher orders too.

The lift itself is not bounded over all real phases. Work on one period;
Theta(x,varphi+1)=Theta(x,varphi)+2pi. Positive-order slow derivatives and
pulled-back periodic functions have uniform bounds over every phase.
Consequently the exact loop and its zero-mean periodic primitives Acal,Bcal
have every fixed joint (y,eta,varphi) derivative bounded, with constants
depending on the frozen smooth input, cap, floors and delta_L. Periodic
integration commutes with slow differentiation; phase derivatives recover
the loop sources. Flat radial collars make the zero extensions smooth.

This proves smoothness on a REAL neighborhood of the closed eta interval,
including +-1, by compactness and the positive floors. It does NOT assert
global complex-eta analyticity: zeta(vs(eta)) can be a nonanalytic flat step.
The original analytic core collar is unchanged by these later supports.

## 3. Sharp fixed-order frequency powers, with full pressure memory

Hold the input, loop and primitives fixed before varying N>=N0, where N0
is the already accepted low-order threshold. For any smooth periodic Acal,

\[
 \partial_y^r\partial_\eta^k
 {Acal(y,eta,Ny)\over N}
 =\sum_{j=0}^r {r\choose j}N^{j-1}
 (\partial_y^{r-j}\partial_\eta^k\partial_\varphi^jAcal)(y,eta,Ny).
                                                               \tag{4}
\]

Thus the sum of its joint input seminorms bounds the right side by
2^r C_rk N^(r-1), also for r=0. Phase is independent of eta.
For E_N=E exp(Acal/N), U_N=U+Bcal/N, the exponential chain rule has
at least one small factor in every difference term. A product with j>=1
factors costs no more than N^(r-j). Use the retained value bound on
exp(Acal/N), not exp of a high-order norm. This proves

\[
 |\partial_y^r\partial_\eta^k\Delta(E,U)|\le C_{rk}N^{r-1}.
                                                               \tag{5}
\]

For r=0 this is O(1/N); r=1 can be order one and r>=2 can grow.
It is not uniform over r or over all other frozen parameters.
Each original density is a fixed smooth weight times U, E, UE,
U^2-E^2/2, or E^2. Product differences satisfy the same power (5).
Keep all incoming moments; use actual E_N,U_N in the forward densities.
For every one of the five moments and Pi=Pi_axis+Cp, integration gives

\[
 |\partial_y^r\partial_\eta^k\Delta(M,I,J,S,Cp,Pi)|
       \le C_{rk}N^{\max(r-2,-1)}.                    \tag{6}
\]

At r=0 integrate a value error over fixed support length. At r>=1,
differentiate the density r-1 times, giving N^(r-2). The first derivative
therefore remains O(1/N). For r=0,1 the exponent is -1. This is an
extra radial gain for MOMENTS, not for the fields. It uses no oscillatory
cancellation or pressure reset. Ordinary X derivatives retain every term
of X^-r product_{j=0}^{r-1}(D_y-j), and X has a positive compact floor.

The original pressure-coordinate map uses one MORE eta derivative of
moments and fields, and positive reciprocals X,E,sqrt(2X)E,L. All fixed
angular errors of p_s are O_k(1/N). Their radial errors can instead have
the field power N^(r-1), because the formula also contains E^-1 and U.
Shear-versus-loop errors likewise obey O_k(1/N) in angular norms; they
are NOT errors of the shear relative to the original unmodulated shear.

## 4. Same accepted root: invert its pointwise Jacobian, not every Ck ball

The first reserved patch has U0=0,E0=e(eta)exp((-.5-lambda)y),
lambda>0 FIXED, and disjoint U/E bump supports. Scaling both field changes
by this SAME original e makes the transformed moment map

\[
 Bc+Q(c,c)=z_N(eta)                                    \tag{7}
\]

independent of eta in B,Q. Keep normalization derivatives of e and the
lambda^-1 transformed target cost; neither can be dropped. Q retains the
U-square and E-square energy terms and E-square Cp term. Disjoint supports
make the J cross term zero, not either square. Use the exact unit-mass bump
integrals, not their quadrature approximations. The retained estimates give
beta=||B^-1||<=1000, q=||Q||<=10^5 and the accepted C1 root
||c||C1<=2beta D/N, D=H^16, with 4beta^2 q D/N<1/2.

Therefore for the SAME root the pointwise Jacobian obeys

\[
 J_c=B+2Q(c,\cdot),\quad
 \|B^{-1}(J_c-B)\|\le2\beta q\|c\|_\infty<1/2,
 \qquad \|J_c^{-1}\|\le2\beta.                       \tag{8}
\]

Pointwise IFT and uniqueness in the original small ball give a smooth root
on a real neighborhood of [-1,1], IDENTICAL to the accepted C1 root.
No contraction in Ck is needed, and N is not reselected with k.
For k>=1 let d_k=sup|c^(k)|/k!, z_k=sup|z_N^(k)|/k!. Differentiating
(7) gives the explicit normalized recurrence

\[
 d_k\le2\beta\left[z_k+q\sum_{j=1}^{k-1}d_jd_{k-j}\right].     \tag{9}
\]

The two highest terms Q(c,c^(k)),Q(c^(k),c) belong to J_c, NOT the
lower sum. Equation (6), the fixed normalizers and lambda>0 give
z_k<=Z_k/N for every fixed k. Set D0=2beta Z0 and inductively
Dk=2beta[Zk+q sum_{j=1}^{k-1}Dj D(k-j)]. Then d_k<=Dk/N for
all N>=N0 (use N^-2<=N^-1), with constants independent of N, but not
uniform in k or lambda. The elementary scalar model c+c^2=z has a
square-root branch point at z=-1/4: small C1 size does not imply a common
smallness bound for all high derivative norms.

Fixed bump widths now give all mixed repair field jets O_rk(1/N), without
losing inverse width powers. With the incoming modulation memory retained,
the FINAL repaired fields and moments retain (5)--(6). At the right endpoint
ALL five moments match exactly as functions of eta. All their eta
derivatives therefore match; fields already agree. The axis datum is the
same, so full pressure and its integrated coordinates agree at EVERY later
radius, not just at a test angle. Restoring only four rows leaves a defect
c(eta) in pressure and 4Aphys eta c-(1-eta^2)c_eta in Ns.

## 5. What joins to the cones, and what remains open

Only the original finite-order acceptance criterion is needed for the
state/cone comparison. Its Cstate=D=H^16, Ccorr=A^128, radius A^-128,
epsilon=G^-64 and N=1+floor(H^32) are unchanged. The all-order arguments
above identify the SAME loop, fields and repair as those low-order estimates;
they impose no extra infinitely many smallness conditions. Thus, relative
to the inherited compact input floors/gaps and boundary comparisons, their
joined leading-cone conclusion persists. This is not a new independent
derivation of the inner/outer continuum margins.

To turn this MAP theorem into a full independent finalized-profile
certificate, the next obligation is to audit the remaining annular and
exterior pre-modulation all-order inputs and exact global pressure/moment
identities, including supports and weighted edges. The parent supplies only
the inner/shape portion. Subsequent uniform physical pulses, means and the
complete residual budget cannot be promoted unconditionally before that
join. There is no unconditional construction or blow-up claim here.

Checks use exact rational jets through orders six/eight: variance-root
coefficients are recovered independently by Lagrange inversion; a moving
phase inverse is recovered by both homogeneous-degree inversion and Newton
jets; five quadratic repair jets are compared with independent fixed-point
series; phase substitution and the extra moment radial gain are checked.
They are manufactured controls, not actual giant fields or continuum bounds.
The original finite-order acceptance evidence is replayed separately.
Generated records snapshot cached diagnostics before exposing them; a mutation
test checks that passed gates cannot hide edits to their canonical reviews.

```sh
python3 research/navier_stokes_cascade/results/modulation_repair_all_order/test_audit.py -v
python3 research/navier_stokes_cascade/results/modulation_repair_all_order/audit.py \
  --verify-record research/navier_stokes_cascade/results/modulation_repair_all_order/evidence.json \
  --output /tmp/modulation_repair_all_order_fresh.json
```
