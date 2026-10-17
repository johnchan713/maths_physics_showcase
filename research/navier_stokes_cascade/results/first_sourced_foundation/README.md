# First sourced-pulse foundation: full complex covariance and current pressure

Parent: `9dbaa3cf768f3009f4f17d9b5c5b47e4df613a6b`.
The [regular-source foundation](../regular_source_foundation/README.md)
supplies the initial supported source at B=1-kappa, kappa=10^-5. This
checkpoint carries that stronger estimate through the SAME existing
signed/mean prefix and first forward solve, without changing its source,
fields, operation order, parameters or cutoffs. It then independently joins
the complete multiharmonic covariance to the recomputed mean pressure and
three defects. Historical sources and frozen evidence remain unchanged.

Result, relative to the retained exact-core and finite background/pulse
existence and gap lemmas: the actual first sourced correction has all fixed
weighted coefficient and physical jets on the unchanged accepted-band
domain. Its complete covariance change has order C=3/2-kappa. The new
regular nonzero residual and FULL tangential mean rows have sufficient
order D=3/2-2kappa, whereas their Haar averages and compatibility defects
have order C. They must be recomputed; their stronger pre-wave orders do
not persist. This is ONE finite-step foundation, not independent closure
of every later complex state or an unconditional complete residual theorem.

Source under audit: [supplied manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(7.5)--(7.16), (7.38)--(7.42), (8.12)--(8.27), (9.1)--(9.10), SHA256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its final conclusion is not an assumption. Analytical derivations below
are separate from the exact manufactured controls and exponent bookkeeping.

## 1. Preserve the actual pre-solve prefix, rather than invent a new state

Keep Md64, fixed h,N, positive lambda, width C^-20, unit viscosity,
primary roots, radius, covers, representatives and carriers. Write
epsilon=Q^h, S=ell^2. All classes include EVERY fixed ordinary coefficient
derivative, polynomial S costs and fixed inverse edge-distance powers.
W has sqrt(zeta)P, M has zeta; P<=1. Dr costs kappa, Dz and slow time
gain one. Angular average and auxiliary Haar average are different.

Do not skip the historical real signed correction and its following mean
relaxation. Put alpha=2-3kappa, H=alpha-kappa, mu=1-kappa. From the
actual initial mean state, its current signed target is M_alpha, real
and auxiliary-independent. Exact weighted moment subtraction makes its
radial primitive compact. The ORIGINAL primary inverse is

\[
d=H_{pulse}^{-1}(\Xi/\epsilon),\qquad
\delta a_\sigma=d_\sigma/(2\sqrt{y_\sigma}). \tag{1}
\]

The parent primary foundation supplies the fixed-radius inverse and
positive original-root bounds. Each differentiated reciprocal-root term
has edge exponent (-n-1/2)+n+1=1/2 after multiplication by d. Thus the
same signed pulse is W_(alpha-1/2), with curl remainder W_(alpha-kappa),
all fixed jets and original smooth supports. No updated root, sign or
small-target condition is used. The squared slow partition recovers the
physical target; no partition is divided out.

The initial REAL covariance identity retains the signed square, both
complete sine curls and all six entries. Its Haar error order is
min(alpha+1-2kappa,2alpha-1)=3-6kappa. Reconstruct full pressure and all
defects, apply the zero-Haar fast map at H to the CURRENT rows, then
the slow five-map at the CURRENT S_alpha defects, as already defined.
The mean-map foundation supplies their all-fixed-order bounds and exact
moments on the convex I4 potential hull, including its nonzero radial
background and both exceptional pressure tails. Fast correlations with
the old fast mean are not zero; the subsequent slow map uses updated
defects. The accumulated slow Haar means remain nonzero of orders (2,1,1).
The conservative pre-source full mean/defect orders remain 3-6kappa and
3-3kappa; total tangential/radial means are M_mu/M_(mu+1).

Using the independently reviewed regular linear and complete transport
operators, the nonzero source orders across this SAME prefix are:

| Retained group | Sufficient order |
|---|---:|
| original primary regular source | B=1-kappa |
| signed linear remainder, primary cross, actual mean cross | alpha-kappa |
| signed square | 2alpha-1-kappa |
| post-signed fast mean--wave increment | H |
| post-signed slow mean--wave increment | alpha |

Every row is at least B. The old complete wave U is W_(1/2), and
U-W0_as is W_B: its primary curl remainder is W_B and its first signed
increment is stronger. This is a sharper estimate for the unchanged
field, not a change of construction. Every source retains its wave factor,
P and smooth full endpoint extension; closed separation removes distinct
label derivative products. The actual input multipliers are +/-1,+/-2.
The exact Gaussian linear tails are excluded from f and remain additive.

## 2. First forward solve: source-normal pressure and every fixed jet

For each actual nonzero f_m in W_B along its entire common lifted path,
solve with zero entrance data

\[
t'+Kt+m^2dt+i km n\pi=-f_m,\quad n\cdot t=0,
\quad d=\epsilon k^2|n|^2,\quad
\pi={i(n\cdot Kt-n'\cdot t+n\cdot f_m)\over km|n|^2}. \tag{2}
\]

Both n'.t and n.f_m stay. Projection gives
t'=A_n t-m^2dt-Proj_(n-perp)f_m; differentiating n.t then gives zero
exactly. Positive/negative conjugate sources yield conjugate t,pi and
potentials, since all background coefficients are real and i/m transforms
correctly. A finite nonzero m is not an inverse angular-component floor.

The parent's propagator bound is ||V_m(v,w)||<=C P(v)/P(w), with extra
nonnegative scalar damping for |m|>1. Hence zero-data Duhamel gives the
W_B bound, retaining sqrt(zeta), because slow/transverse parameters are
fixed along the path. For input polynomial degrees a_j, multiplication
by the frame/projection gives g_j=max_(0<=i<=j)(i+a_(j-i)); sufficient
solution degrees obey

\[
b_0=g_0+1,\qquad
b_j=\max\{g_j+1,\max_{1\le i\le j}(i+1+b_{j-i})\}. \tag{3}
\]

Every requested order uses finitely many input jets. Pulse derivatives
follow the ODE, and ordinary common-torus derivatives cost polynomial S,
not new epsilon powers. Pressure (2) and the fixed potential
i(n cross t)/(km|n|^2) have order C=B+1/2. This now applies to the
ACTUAL bounded/supported f, rather than an arbitrary unverified source.
No giant parameter value, starting band or actual ODE profile is sampled.

A source zero on a full slow/transverse parameter path gives zero response
there by uniqueness. Deck translations commute with the path; no distinct
preimages are averaged. The uncut exit need not vanish. Apply the SAME
pulse cutoff to potential AND pressure before curl and zero extension.
The resulting actual delta w is W_B, with complete curl remainder W_(C-kappa),
real, divergence-free and smoothly supported on the original closed labels.
All scales use the same accepted domain and margin 4Q<qstar. In particular,
an exit value merely bounded by a tiny positive Gaussian is not declared zero.

## 3. Complete complex covariance, including the second harmonic square

Keep actual mean velocities fixed during this wave operation. For each
label combine equal angular frequencies BEFORE taking the angular mean.
Writing old and new cosine/sine coefficients a_m,b_m and c_m,d_m,

\[
\Delta W={1\over2}\sum_{m>0}
(a_m c_m^T+c_m a_m^T+b_m d_m^T+d_m b_m^T
                         +c_m c_m^T+d_m d_m^T). \tag{4}
\]

Old coefficients absent at a frequency are zero. Old U has |m|=1; the
sourced increment has |m|=1,2. Its |m|=2 square contributes even though
its cross product with the old |m|=1 averages to zero. Different labels
instead have zero derivative products by closed separation, not a presumed
frequency cancellation. Every cylindrical tensor entry is retained.

The cross term is M_(B+1/2)=M_C, while the new square is M_(2B) and
2B>=C. Thus Delta W is M_C at every fixed coefficient derivative.
No initial real-primary cosine/sine shortcut applies: the sourced pulse
is generally complex in its factored coefficients. Old sine-curl times
new sourced sine-quadrature terms can survive. Delta W, as a difference
of positive covariance tensors, need not be positive semidefinite;
only the new square is. For example removing e_z cos(theta), a locally
divergence-free axial wave, gives Delta Wzz=-1/2. A radial cutoff independent
of z preserves that divergence identity. No positivity of updated amplitudes is required.

## 4. Reconstruct the actual pressure and ALL current defects

Use the same complete background pressure p_B, including Pi1=K1-mu0^2/(4X)
and Pi2=K2 on I4. It cancels only against that SAME full base balance,
not leading p0. A coefficient radially constant in X can have a nonzero
physical axial gradient. The new wave pressure is angular-mean zero,
not pointwise zero. Both facts leave the following mean pressure change:

\[
\Delta g_r=-D_1\Delta W_{rr}-D_z\Delta W_{rz}+\Delta W_{\theta\theta}/R,
\quad\Delta P=\int\overline{\Delta g_r}\,dR,
\quad\Delta p_m=T_0(\Delta g_r-\rho\Delta P). \tag{5}
\]

This is subtraction of two FULL recomputed states with the SAME nonzero
old means, not a reset to zero mean. Pointwise Delta g_r,Delta p_m are
M_(C-kappa); Haar translation invariance removes the displayed graph
derivative, so their averages are M_C. The actual A0 pressure error has
zero Haar mean and remains a flat pointwise residual. Exactly

\[
\int R\overline{\Delta p_m}
 =-\tfrac12\int R^2\overline{\Delta g_r}+c_\rho\Delta P,
\quad c_\rho=\tfrac12\int R^2\rho,
\]
\[
\Delta J_\theta=\int R^2\overline{\Delta W_{\theta z}},\qquad
\Delta J_z=\int R\overline{\Delta W_{zz}}
             -\tfrac12\int R^2\overline{\Delta g_r}. \tag{6}
\]

Keep the WHOLE partial_Z(c_rho Delta P). The physical unit-integral bump
rho=s^-1/2 rhohat(R/sqrt(s)) moves with s=q/Q. Its derivative times
Delta P is generally nonzero. Current full tangential row changes are
D2 Delta W_rtheta+Dz Delta W_theta z and
D1 Delta W_rz+Dz(Delta W_zz+Delta p_m), respectively. Their orders are
M_(C-kappa), and their Haar averages M_C. All three current defects are
S_C, not their stronger pre-wave values. Both old mean velocity moments
and old Haar mean velocities are unchanged EXACTLY by the nonzero-angular
wave addition. Total p_m stays M_mu since C-kappa>mu. The radial residual
is the CURRENT -rho P plus the actual flat A0 term, not the old P.

## 5. One finite regular step, physical jets, and the next gate

The exact linear cancellation keeps
psi times the uncut remainder PLUS (1-psi)f+T_psi, with the complete
potential/pressure commutator T_psi from the parent. Only the first
summand is a subsequent regular source. Both tails remain exact additive
residuals: a flat early source can have a nonflat forward response.
Actual flat cutoff/A0/A1 velocities remain in EVERY nonlinear product;
the base flat residual remains a separate summand without assumed shell
support or weight. Reclassification by mean reconstruction uses exact
subtraction, not double counting or deletion.

The new regular NONZERO groups have sufficient orders

| Retained group | Order |
|---|---:|
| full linear remainder | B+1/2-kappa=D |
| complete old/new wave interaction | B+1/2-kappa=D |
| new wave square | 2B-kappa |
| actual total mean/new wave interaction | B+mu-1/2=D |

The regular nonzero rows are therefore W_D, with |m|<=4, full source
endpoint extension and the original Gaussian/edge weights. Full mean
rows are M_D and defects S_C as above. These are ONE finite normalized
step, gain D-B=1/2-kappa, relative to the retained analytic dependencies;
they are NOT an unconditional global complete-residual improvement.
The complete old-wave error is now only W_B, not W_D. Do not assign the
improved residual order to the velocity error used in the next inverse.

For every fixed spatial/time jet count a,b, sx=1/2+h/2, st=1+3h/2,
the SAME physical transfer gives powers (up to fixed weighted S factors)

\[
\delta u: -A+hB-a s_x-b s_t,\quad
\delta\pi: -2A+hC-a s_x-b s_t,\quad
\delta\mathcal A: 1/2-A+hC-a s_x-b s_t. \tag{7}
\]

These are uniform over every accepted band/label/copy, not in derivative
order, h, N or the actual huge construction constants. The zero-jet
regular residual power -2A-1/2+hD is still negative (also at the
manufactured h=1/200): normalized order gain is not physical absolute
decay, fixed-scale contraction or numerical evidence of blow-up.

Next use the CURRENT full mean rows and defects to form the generic
complex signed target at order C, with old error rho=B, not an initial
real shortcut or the residual order D. Then reapply fast/slow mean maps
with pressure/defects recomputed after each operation. Every later target,
source family and full finite residual induction still needs its independent
join before any all-family or global residual claim. Infinite realization,
forcing/energy admissibility and unconditional blow-up are not established.

## 6. Controls and reproduction

New independent rational controls check forced complex moving-normal jets
and conjugate harmonics; complete curls at harmonics 1 and 2 via Fourier
zero-mode convolution; and full before/after mean pressure, axial moments,
viscosity and nonzero old means. The forced-jet fixture checks compatibility,
not an actual supported giant pulse. Polynomial potentials are only C5
after zero extension (their curls C4), not the analytical C-infinity fields.
All new fixtures are diagnostics, not continuum suprema or certificates.

```sh
python3 research/navier_stokes_cascade/results/first_sourced_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/first_sourced_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/first_sourced_foundation/evidence.json \
  --output /tmp/first_sourced_foundation_replay.json
```

Use fresh output paths; never overwrite historical evidence.
