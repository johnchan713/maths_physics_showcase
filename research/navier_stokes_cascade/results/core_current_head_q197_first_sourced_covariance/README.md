# Q197: first sourced inverse, complex covariance and current means

Parent: Q196, `6803cff8a2f77fcd0ab721e02768f938cda48e53`.
The live branch, Q196, its active CI, the explicit sourced inverse, the frozen
first-sourced covariance/pressure foundation and Q184's logarithmic-loss
correction were read first. Historical evidence, the selected profile,
pressure gauge, roots, supports, accepted scales, frequencies and unchanged
`C^-20` activation width are preserved. Q196 CI had no failed step and was
still in progress when this checkpoint began.

Q196 supplies the complete initial generated source at log-renormalized order
`B=1-kappa`. This checkpoint carries that source through the actual existing
prefix, solves the first forced pulse, and independently checks the complete
complex covariance, compact pressure and current mean orders. The conclusion
is conditional on the named profile/frame inputs and is one finite step.

## 1. Preserve the signed/mean prefix before the sourced solve

The sourced solve is not applied directly to the bare Q195 state. The
construction first keeps the existing real signed correction with

\[
 \alpha=2-3\kappa,
 \qquad H=\alpha-\kappa=2-4\kappa,
\]

and its recomputed fast/slow mean relaxation. The complete nonzero source
orders through that prefix are

\[
\begin{array}{c|c}
\text{original primary-plus-initial-mean source}&B=1-\kappa\\
\text{signed linear, primary cross, actual mean cross}&\alpha-\kappa\\
\text{signed square}&2\alpha-1-\kappa\\
\text{post-signed fast mean--wave}&H\\
\text{post-signed slow mean--wave}&\alpha.
\end{array}                                             \tag{1}
\]

Every added row is stronger than `B`, so the Q196 source remains the limiting
input. This ordering retains the original positive roots, the signed square,
nonzero old fast correlations and the full pressure recomputation. Skipping
the prefix would describe a different state.

## 2. Exact forced constraint and zero-entry solve

For each nonzero input harmonic `m in {+/-1,+/-2}`, solve on the unchanged
lifted pulse path

\[
 t'+Kt+m^2dt+ikm n\pi=-f_m,\qquad n\cdot t=0,
 \qquad t(0)=0.                                       \tag{2}
\]

Differentiating the constraint and dotting (2) with `n` gives

\[
 \pi={i\,[n\cdot Kt-n'\cdot t+n\cdot f_m]
          \over km|n|^2}.                            \tag{3}
\]

Both the moving-normal terms and the source-normal term remain. Substitution
of (3) makes `(n dot t)'=-m^2d(n dot t)`, so zero entrance data preserve the
constraint exactly. A source-normal pressure cannot be assigned the stronger
potential damping: for a purely normal source one can have `t=0` and nonzero
`pi`.

In moving-frame coordinates, after dividing by the positive reference
envelope `P`, the propagator satisfies

\[
 \|V_m(v,w)\|\le e^H{P(v)\over P(w)}
 e^{-(m^2-1)d_*(v-w)}.                               \tag{4}
\]

Thus the residence factor is `L` for `|m|=1` and at most
`[(m^2-1)d_*]^-1` for `|m|>=2`. Differentiating the ODE gives a finite
polynomial `S` degree at every fixed order; it does not change the epsilon
order or accepted domain. The Q196 mixed entrance compatibility makes all
fixed pulse/parameter jets legal.

With

\[
 B=1-\kappa,\qquad C=B+{1\over2}={3\over2}-\kappa,
\]

the uncut tangential response and the complete cut sourced velocity have
order `B`; pressure and vector potential have order `C`; the complete curl
remainder has order `C-kappa`. The uncut exit need not vanish. The same pulse
cutoff is applied to potential and pressure before the complete curl and zero
extension.

## 3. Complete equal-frequency covariance convolution

Hold the current mean velocities fixed during this wave operation. For each
label combine equal positive angular frequencies before averaging. If
`(a_m,b_m)` are old cosine/sine coefficient vectors and `(c_m,d_m)` are the
new ones, then

\[
 \Delta W={1\over2}\sum_{m>0}
 \bigl(a_mc_m^T+c_ma_m^T+b_md_m^T+d_mb_m^T
       +c_mc_m^T+d_md_m^T\bigr).                    \tag{5}
\]

Old coefficients absent at frequency `m` are zero. The old wave has harmonic
one, while the sourced correction has harmonics one and two. Therefore the
harmonic-two square in (5) contributes even though its cross with the old
harmonic-one wave averages to zero. Different labels vanish by separated
closed support, rather than by a frequency assumption.

Equation (5) retains all six entries

\[
 \Delta W_{rr},\Delta W_{r\theta},\Delta W_{rz},
 \Delta W_{\theta\theta},\Delta W_{\theta z},\Delta W_{zz}.
\]

The old/new cross has order `B+1/2=C`; the new square has order `2B`, and

\[
 2B-C={1\over2}-\kappa>0.                            \tag{6}
\]

Hence `Delta W` is in `M_C`. The increment need not be positive semidefinite;
only the new square is. No real-primary cosine/sine shortcut is reused for
this generally complex sourced field.

## 4. Recompute compact pressure, axial moment and defects

Subtract the two full current states with identical nonzero old means. The
complete radial source and compact pressure changes are

\[
 \Delta g_r=-D_1\Delta W_{rr}-D_z\Delta W_{rz}
              +{\Delta W_{\theta\theta}\over R},
\]
\[
 \Delta P=\int\overline{\Delta g_r}\,dR,
 \qquad\Delta p_m=T_0(\Delta g_r-\rho\Delta P).       \tag{7}
\]

The same background pressure with both `Pi1/Pi2` tails is subtracted; the
new wave pressure is zero only in angular mean. Pointwise (7) has order
`C-kappa`, while Haar translation invariance gives averaged order `C`.
Exactly,

\[
 \Delta J_\theta=\int R^2\overline{\Delta W_{\theta z}}\,dR,
\]
\[
 \Delta J_z=\int R\overline{\Delta W_{zz}}\,dR
             -{1\over2}\int R^2\overline{\Delta g_r}\,dR,       \tag{8}
\]

and the axial pressure identity differentiates the whole product
`c_rho Delta P`, including `(partial_Z c_rho) Delta P`. All three recomputed
defects belong to `S_C`; the full tangential current rows belong to
`M_(C-kappa)`, and their Haar averages belong to `M_C`. Old mean velocity
moments are unchanged by adding a nonzero angular wave.

## 5. Post-source regular rows and pure-power ledger

The next complete nonzero regular groups have orders

\[
\begin{array}{c|c}
\text{full linear remainder}&B+1/2-\kappa\\
\text{complete old/new interaction}&B+1/2-\kappa\\
\text{new wave square}&2B-\kappa\\
\text{actual total mean/new wave}&B+(1-\kappa)-1/2.
\end{array}
\]

Their minimum is

\[
 D={3\over2}-2\kappa.                                \tag{9}
\]

This is the first point at which the `D` order named prematurely in Q195 is
available. Fixed polynomial `S` losses remain in the log-renormalized
classes. Applying Q184's one-`kappa` absorption only after collecting each
fixed request gives the parallel pure-power ledger

\[
 \widetilde B=1-2\kappa,qquad
 \widetilde C={3\over2}-2\kappa,qquad
 \widetilde D={3\over2}-3\kappa.                    \tag{10}
\]

The physical zero-jet regular residual still has exponent

\[
 -2A-{1\over2}+hD=-{15025001\over10000000}<0.        \tag{11}
\]

Thus the normalized gain is not absolute physical decay, a numerical PDE
verification or a complete residual theorem.

## 6. Remaining boundary

The next state must form the generic complex signed target from the current
full rows and `S_C` defects. Its old wave error has order `B`, not the improved
residual order `D`. The next checkpoint must derive the weighted original-root
update and then recompute fast/slow means and full pressure after each
operation. Reusing the initial real covariance shortcut or preserving the
pre-wave defects would be invalid.

The exact tests check (1)--(11), forced constraint algebra, finite degree
recurrences, equal-frequency covariance and immutable provenance. They are
not a continuum proof. No unconditional finite-cycle residual, admissible
force, comparison exclusion or Navier--Stokes blow-up is claimed.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q197_first_sourced_covariance/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q197_first_sourced_covariance/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q197_first_sourced_covariance/evidence.json \
  --output q197-fresh.json
```
