# Generic complex old-wave signed covariance correction

Status: `current-complex-wave-signed-covariance-substep-audited`.
Parent: `5e7f1ee0103be6f9ba9ad6b1f9e4fec64efffb31`.

Using the explicitly retained continuum lemmas, this checkpoint realizes
the next signed covariance correction about the CURRENT complex wave.
Its auxiliary-averaged tangential residual has order 2-6kappa. The generic
old-error cross term and signed square both limit this estimate. Full
unaveraged means and all three defects are recomputed, NOT relaxed here.
The nonzero-angular regular residual stays at 3/2-6kappa, plus exact
designated tails. This is one signed substep, not an all-cycle residual
theorem, an infinite construction, or blow-up. Diagnostics are not
continuum proofs or evaluations of the actual giant first band.

All previous source files, frozen evidence, GOAL.md and fixed setup remain
unchanged. The proposed manuscript is the [supplied source under audit](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(7.39)--(7.42), (8.12)--(8.27), (9.3)--(9.6); its conclusion is not assumed.
PDF SHA256: `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The exact formulas below use the frozen [wave interface](../physical_wave_interface/README.md),
[actual-band inverse](../actual_band_pulse_bounds/README.md),
[real signed update](../real_signed_covariance_update/README.md) and
[current sourced state](../first_sourced_wave_update/README.md).

## 1. Correct the observed replay failure before advancing

On 2026-10-16 the parent's [CI run 38929673498](https://github.com/johnchan713/maths_physics_showcase/actions/runs/38929673498)
FAILED in `test_frozen_evidence_reproduces`, at
`full_viscous_potential_pressure_cutoff.maximum_relative_commutator_error`.
The job log shows that all other tests and the commutator acceptance gate
passed. The inherited comparison has relative tolerance 1e-7 but absolute
tolerance 1e-70. Comparing roundoff-sized errors near 6e-15 that way demands
agreement to roughly 6e-22; different NumPy/BLAS arithmetic need not supply it.
This was a validation portability error, not a failed continuum estimate.
The previous run's local pass is not a claim that its remote CI passed.

`replay.py` adds an explicit compatibility adapter for ONLY three fields of
that one manufactured Cartesian fixture: maximum relative commutator error,
maximum curl-product error and maximum divergence error. At exactly those
paths both old AND replay values must be finite, nonnegative and <=2e-12.
This ceiling is stricter than the unchanged original gates (2e-10, 2e-10,
2e-9). It is not a new tolerance for the actual PDE or an interval bound.
All other scalars, checks, source hashes, keys, types and historical
comparisons remain unchanged. Missing gates, NaN, infinity, negative error,
or excessive roundoff still fail. Frozen files are not regenerated to hide
the mismatch. The adapter runs the original 35 tests and original audit;
CI's parent step now invokes it. New audits recursively replay the same
parent through this adapter. The correction is documented HERE rather
than modifying the parent's historical record and breaking its hashes.

## 2. Current inputs and the fixed signed target

Keep the same Md64, N, fixed 0<h<1/100, positive lambda, width C^-20,
labels, closed enlarged supports, common torus, background, qstar and
positive ORIGINAL primary amplitudes. Below alpha and rho denote orders,
not the setup's large constant or the radial moment bump. Use
epsilon=Q^h, S=ell^2 and kappa=1/100000 with all fixed coefficient jets
and the frozen W/M/S edge weights, graph rules and Gaussian envelope P.

The current state has

\[
B=1-3\kappa,\quad \alpha=B+\tfrac12,\quad D=\tfrac32-6\kappa,
\qquad U\in W_{1/2},\quad U-W_0^{as}\in W_\rho,\quad \rho=B.
\tag{1}
\]

Tangential means have order mu=1-kappa, radial mean mu+1. Full tangential
mean residuals have order alpha-kappa, their auxiliary Haar averages order
alpha, and three defects S_alpha. The old slow Haar mean is nonzero. Both
background pressure tails, full axial wave covariance, mean pressure and
actual A0/A1 flat velocities remain. All designated Gaussian linear tails
are exact separate residuals, never sources for a forward pulse solve.

From the CURRENT Haar mean equations and moments form

\[
\Xi_\theta=-T_2[\bar E_\theta-\sigma_\theta\epsilon\partial_ZJ_\theta],
\quad
\Xi_z=-T_1[\bar E_z-\sigma_z\epsilon\partial_Z(J_z+c_\rho P)]. \tag{2}
\]

The exact integrated identities make both bracketed weighted integrals
zero. Compact radial inversion is therefore exact; it costs no graph
derivative after averaging. Moment-bump sources have order alpha+1, so
Xi is M_alpha, real and auxiliary-independent, with the original smooth
shell/slow extensions. The WHOLE c_rho P product is differentiated.

Use the original audited H and positive roots a_sigma=sqrt(y_sigma):
d=H^-1(Xi/epsilon), delta a_sigma=d_sigma/(2a_sigma). The inverse is
linear, signed, and does not require updated amplitudes to be positive.
No square root of y+d or division by a slow partition is taken. The
reciprocal-root edge calculus retains sqrt(zeta) for every fixed derivative:
each term has exponent (-n-1/2)+n+1=1/2. The fixed radius-dependent costs
of H^-1 stay in the constants, not a claim uniform as radius tends to zero.
The squared slow partition restores the physical target exactly across
overlapping bands. These are the SAME analytical inverse/pulse premises
as the frozen real update, not a new numerical inverse assumption.

Cut potential AND wave pressure before curl. The principal signed field
L has order b=alpha-1/2=B; its complete curl velocity V=L+R has that order,
with R in W_(alpha-kappa) and pressure/potential order alpha. Its target is
auxiliary-independent, so it is a homogeneous principal pulse, not a new
forward sourced solve. All cutoff commutators are retained exactly.

## 3. The generic complete-wave identity, not the initial shortcut

For FULL tensors set C(u)=<u tensor u>_theta and
B(u,v)=<u tensor v+v tensor u>_theta. Let E=U-W0^as. Exactly,

\[
C(U+V)-C(U)=\Lambda+B(E,L)+B(U,R)+C(V),\qquad
\Lambda=B(W_0^{as},L). \tag{3}
\]

The selected two radial/tangential entries of the auxiliary average of
Lambda equal Xi; its other entries need not vanish. Equation (3) retains
ALL six entries and all quadratic terms. General complex old amplitudes
have both cosine and sine quadratures. Their old cosine error times L
does NOT vanish on angular averaging. Real physical velocity is not the
same as the phase-locked real coefficient structure of the INITIAL step.

The generic error orders after Haar averaging are

| Exact group | Order |
|---|---:|
| B(E,L) | rho+alpha-1/2 = 2-6kappa |
| B(U,R) | alpha+1/2-kappa = 2-4kappa |
| C(V), including the full signed square | 2alpha-1 = 2-6kappa |

Thus nu=2-6kappa, and the gain over alpha is 1/2-3kappa. In general the
estimate gives a strict gain precisely when

\[
\min(\rho-\tfrac12,\tfrac12-\kappa,\alpha-1)>0. \tag{4}
\]

At rho=1/2 or alpha=1 this bound gives no positive gain; tests are not a
way around that inequality. The bound is conservative, not a lower bound
or a claim that every individual pulse saturates it. Extra angular
orthogonality can only remove terms: the current old |m|=2 components do
not average against the fresh |m|=1 signed field. Different labels have
zero derivative products by CLOSED support separation, not by treating
distinct overlapping phases as the same. The complete tensor CHANGE has
order alpha because Lambda does; it is not all M_nu.

## 4. Reconstruct pressure, defects and every generated nonzero term

Hold actual mean VELOCITIES fixed, recompute all six W entries, then the
full radial source and compact pressure. With Delta denoting this substep,

\[
\begin{split}
\Delta\bar g_r&=-(\partial_R+1/R)\Delta\bar W_{rr}
 -\epsilon\partial_Z\Delta\bar W_{rz}+\Delta\bar W_{\theta\theta}/R,\\
\Delta P&=\int\Delta\bar g_r,\quad
\partial_R\Delta\bar p_m=\Delta\bar g_r-\rho_{bump}\Delta P,\\
\int R\Delta\bar p_m&=-\tfrac12\int R^2\Delta\bar g_r+c_{\rho_{bump}}\Delta P,\\
\Delta J_\theta&=\int R^2\Delta\bar W_{\theta z},\quad
\Delta J_z=\int R\Delta\bar W_{zz}-\tfrac12\int R^2\Delta\bar g_r.
\end{split} \tag{5}
\]

Use unaveraged sources for unaveraged pressure. Full pressure/source changes
are M_(alpha-kappa); their Haar averages are M_alpha. After (2), the averaged
tangential rows consist of radial divergence of the error in (3), old
moment bumps of order alpha+1, new axial covariance flux and changed axial
pressure of order alpha+1. Haar averaging removes graph derivatives, so
BOTH averaged tangential rows reach M_nu. Never omit W_zz, W_theta z,
Delta p_m, or the moving-bump derivative. All three defects are recomputed
and are conservatively S_alpha; full tangential means stay M_(alpha-kappa).
Neither is set to zero or promoted to nu. The radial residual is the new
moment obstruction plus the actual flat reconstruction term.

Nonzero waves preserve both mean-velocity moments, mu/mu+1 velocity bounds,
and the nonzero old slow mean. Total p_m stays M_mu. Both exceptional p_B
tails cancel only inside the SAME full base equation. Actual flat A0/A1
velocities remain in EVERY later nonlinear product. Any flat ledger moved
by pressure reconstruction is updated by exact subtraction, not counted
twice. The new homogeneous potential/pressure cutoff commutator is added
to the exact excluded tail ledger, not erased or forwarded to an inverse.

Complete same-label curl incompressibility gives n.a=i div_*a/(km), so
wave-wave products lose kappa, whereas mean-wave transport can lose 1/2.
This substep's generated nonzero regular orders are

| Full term | Order |
|---|---:|
| Full linear remainder, viscosity and pressure included | alpha-3kappa = D |
| Complete old/new wave transport | alpha-kappa |
| Actual mean/new wave transport | b+mu-1/2 = alpha-kappa |
| New signed wave square | 2alpha-1-kappa |

All retain a wave factor, sqrt(zeta), P (using P^2<=P), the same closed
slow/transverse/temporal supports and smooth endpoint zero extension.
Finite harmonic multipliers stay <=4 in the residual. The old regular
nonzero residual was W_D: its sum with these new terms is still W_D,
NOT improved. After this wave change U_new-W0^as remains W_B.

For every fixed spatial/time jet count a,b the physical signed increments
have powers -A+h(alpha-1/2)-a sx-b st for velocity,
-2A+h alpha-a sx-b st for pressure and
1/2-A+h alpha-a sx-b st for potential, A=1/2+h,
sx=1/2+h/2, st=1+3h/2. Bounds are uniform across every accepted late band,
with the same enlarged-chart requirement 4Q<qstar and fixed polynomial S
losses. At h=1/200 the zero-jet powers are -0.50000015, -1.00250015,
+0.00249985. The unchanged residual physical exponent is -1.5025003;
normalized order gain is NOT absolute physical decay or fixed-scale
contraction. Stage/derivative constants are not uniform in infinite iteration.

## 5. Checks and next obligation

The fixture rotates an exact complete old curl into both quadratures, adds
a new nonconstant signed potential, and verifies (3) by independent exact
quarter-phase averages. It detects all three omitted error groups, checks
all six tensor entries, divergence, compact pressure, both integrated flux
identities, nonzero axial pressure and moving-bump terms. Fixtures have
finite endpoint regularity and are NOT actual global fields. The fixed
inverse/factor2 and squared-partition checks are replayed from the earlier
signed checkpoint. The parent is replayed with the documented adapter;
655 inherited source hashes are checked without changing historical bytes.

Next apply current fast mean inversion, recompute pressure and defects,
then the slow five-map and recompute AGAIN. Old slow Haar means and fast
correlations cannot be replaced by their initial values. Later generated
sources, every finite cycle, infinite potentials, forcing/energy, full
residual closure and blow-up remain unverified.

```sh
python research/navier_stokes_cascade/results/complex_signed_covariance_update/test_audit.py -v
python research/navier_stokes_cascade/results/complex_signed_covariance_update/audit.py \
  --verify-record research/navier_stokes_cascade/results/complex_signed_covariance_update/evidence.json \
  --output complex_signed_replay.json
python research/navier_stokes_cascade/results/complex_signed_covariance_update/replay.py --tests
```

Choose fresh output paths. No frozen record is overwritten.
