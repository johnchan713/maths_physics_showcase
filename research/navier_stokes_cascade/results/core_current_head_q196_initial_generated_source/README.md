# Q196: corrected initial generated-source stage and quotient bound

Parent: Q195, `fb74dd3975157277afb5d41f4e8b2033ff18d13e`.
The live branch, Q195, its active CI, the regular-source and source-envelope
foundations, Q184's logarithmic-loss correction, the first-sourced provider
and their frozen evidence were read first. Historical evidence, the selected
profile, pressure gauge, roots, supports, accepted scales, frequencies and the
`C^-20` activation width remain unchanged. Q195 CI had no failed step and was
still in progress when this checkpoint began.

Q195 correctly closes the complete real-primary covariance and initial mean
update. Its stated next inequality is not the immediate source obligation: it
uses the post-first-sourced order `3/2-2*kappa`, and repeats Q183's exact-power
supremum without the polynomial `S` loss already corrected by Q184. This
checkpoint repairs the stage order and derives the complete source estimate
which follows directly from the Q195 primary-plus-initial-mean state.

## 1. The immediate source order is `1-kappa`

Use the normalized chart operators

\[
 D_r=\partial_R+Md_rR^{d_r-1}L,\qquad
 D_z=\epsilon\partial_Z,\qquad
 t_*=-\epsilon\partial_T+c_iN_t,
\]

with `epsilon=Q^h`, `S=ell^2`, angular carrier
`omega=km`, and the Q194/Q195 complete real primary order

\[
 \alpha={1\over2}.
\]

For a coefficient amplitude `a` and wave pressure `pi`, the complete
cylindrical linear remainder retains slow time, radial and axial transport,
background derivatives, radial and axial pressure gradients, full vector
viscosity, mixed phase/amplitude viscosity, phase divergence and angular
vector-Laplacian connections. The limiting gain over `alpha` is

\[
 \min\left\{1-\kappa,{1\over2},1,
 {1\over2}-\kappa,1-2\kappa\right\}
 ={1\over2}-\kappa.                                  \tag{1}
\]

The source-normal pressure term `n dot f` and the moving-normal terms remain
in the constrained pressure formula. The complete curl remainder is included;
its fast operator does not acquire a second graph derivative merely because
`N_t` is requested as a coefficient jet.

The full primary-plus-initial-mean nonzero source has three groups. Their
orders are

\[
 \begin{array}{c|c}
 \text{complete primary linear remainder}
   &\alpha+{1\over2}-\kappa=1-\kappa\\
 \text{complete primary square}
   &2\alpha-\kappa=1-\kappa\\
 \text{actual initial mean--primary interaction}
   &\alpha+(1-\kappa)-{1\over2}=1-\kappa.
 \end{array}                                         \tag{2}
\]

The square calculation uses the exact complete-curl longitudinal identity,
not principal transversality. The mean interaction retains its genuine
half-frequency cost. Hence the immediate regular source order is

\[
 B=1-\kappa={99999\over100000}.                      \tag{3}
\]

The order

\[
 D={3\over2}-2\kappa={74999\over50000}               \tag{4}
\]

is reached only after solving the first sourced correction and recomputing
its complete covariance and means. Assigning (4) to the source directly after
Q195 skips a construction stage.

## 2. Exact support grammar and pressure memory

Every forwarded linear source contains the smooth temporal cutoff. Every
forwarded quadratic or mean-wave term contains a complete cut velocity,
including radial and axial components, cutoff curls and compact `A0/A1`
velocities. Distinct labels have zero derivative products on their separated
closed supports. The nonzero harmonics at this stage are contained in
`{+/-1,+/-2}`.

The source retains:

* source-normal and moving-normal pressure;
* both `Pi1/Pi2` tails in the same base and mean equations;
* vector and mixed viscosity;
* every cylindrical connection and differentiated graph multiplier;
* both interaction orientations and the complete square; and
* radial, angular and axial velocity.

After the exact principal cancellation, the regular forwarded part is the
cutoff times the uncut remainder. The terms `(1-psi)f+T_psi` are additive
Gaussian/cutoff tails. They are not forwarded through the inverse. Their
cut velocities still participate in every later nonlinear product. A positive
Gaussian value at a pulse endpoint is never treated as support.

The actual smooth cutoffs make every mixed endpoint jet of the forwarded
source vanish. Since `P` is positive at each finite pulse endpoint, the mixed
endpoint jets of `f/P` vanish as well. This is support and quotient algebra,
not a claim that the uncut sourced response has zero exit data.

## 3. Quotient jets and the missing polynomial degree

For every fixed mixed jet `I`, the current all-fixed coefficient bounds give
finite integers `d,r,s` and a constant such that

\[
 \left|D^I\left({f_{initial,m}\over P}\right)\right|
 \le C_I\epsilon^{1-\kappa}S^d\sqrt\zeta\,
       y_a^{-r}y_b^{-s}.                              \tag{5}
\]

The Bell recurrence for derivatives of `P^-1` uses finitely many actual
`D^J log P` bounds. It introduces only another finite polynomial in `S`; no
cutoff, partition or oscillating amplitude is divided out.

When `d>0`, (5) does not imply an all-band supremum at exactly the exponent
`1-kappa`, because

\[
 \epsilon^{-(1-\kappa)}\epsilon^{1-\kappa}S^d
 =\ell^{2d}\longrightarrow\infty.                   \tag{6}
\]

This is the second correction to Q195's displayed frontier. Q184 had already
made the same correction for the later `B_j` stages.

For any `eta>0`, with `D=ceil(d)` and `epsilon=2^{-h ell}`, the exponential
series and `log 2>=1/2` give

\[
 \ell^{2d}\epsilon^\eta
 \le(2D)!\left({2\over h\eta}\right)^{2D}.           \tag{7}
\]

Choosing `eta=kappa` converts (5) into the uniform pure-power estimate

\[
 \left|D^I\left({f_{initial,m}\over P}\right)\right|
 \le C'_I\epsilon^{1-2\kappa}\sqrt\zeta\,
       y_a^{-r}y_b^{-s}.                              \tag{8}
\]

For the exact diagnostic `d=1`, `h=1/200`, `kappa=1/100000`, the explicit
constant in (7) is `3200000000000000`. Its size is irrelevant to asymptotic
existence but emphasizes that this is not a numerical evaluation of the first
accepted band.

## 4. Scope and next obligation

Equations (1)--(5) are an analytical composition of the independently
assembled Q194 field, Q195 complete initial mean state, and the explicit
regular-source operator. Equation (7) is an exact elementary conversion.
They establish the immediate initial generated-source quotient conditionally
on the named profile/frame inputs. The executable tests verify the term and
exponent ledger; they do not establish continuum suprema by finite sampling.

This checkpoint does not independently solve the first sourced wave. The next
obligation is to solve the zero-entry forced ODE from the exact source (5),
retaining `n dot f` in pressure, and then independently rederive on the same
domain:

1. the complete real sourced curl with harmonics `+/-1,+/-2`;
2. all six complex covariance increments, including every new-harmonic
   square;
3. the compact pressure and whole moving axial pressure moment; and
4. the recomputed current defects and the post-source order (4).

No unconditional pulse hierarchy, complete residual improvement, admissible
force, comparison exclusion or Navier--Stokes blow-up is claimed.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q196_initial_generated_source/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q196_initial_generated_source/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q196_initial_generated_source/evidence.json \
  --output q196-fresh.json
```
