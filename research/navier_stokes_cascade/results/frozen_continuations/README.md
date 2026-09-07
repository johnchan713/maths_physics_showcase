# Frozen-field continuations

This is a completed, bounded `N=32/64` continuation campaign. All three
initial fields and viscosity were held fixed. Every grid pair eventually
failed an empirical spatial-convergence gate while critical-norm growth
remained positive. That is a resolution limit, not evidence for singularity
and not proof that these initial fields cannot develop one.

[`manifest.json`](manifest.json) records all 20 stage commands, outcomes,
source/executable hashes and final checkpoint identities. Controls were
`nu=0.02`, `E(0)=10`, initial bandwidth `K=3`, evolution cutoffs `10/21`,
maximum timestep `0.000125`, CFL target `0.4`, and observations every `0.01`.
Primary and denser physical samples use grids `64` and `128` respectively.
Stages were checked at `0.08, 0.10, 0.12, 0.16`; each pair stopped advancing
after a failed stage. Cutoff energy was monitored at every accepted step.

## Results of the grid-pair gates

| Frozen seed | Last passing stage | Fine-grid H1/2 gain there | First failed stage | Failure |
|---|---:|---:|---:|---|
| Previous packet checkpoint | 0.08 | 7.91% | 0.10 | Vorticity disagreement |
| Fresh packet, start 0 | 0.12 | 11.04% | 0.16 | Vorticity and stretching disagreement |
| Tube, start 0 | 0.12 | 5.76% | 0.16 | Vorticity and stretching disagreement |

These are passing **continuation-stage** gates, not complete new independent
trajectory or timestep validations for every field/horizon. The earlier
three-field independent checks were at `T=0.06`. Longer states require their
own confirmation before promotion. Seed names record origins; Fourier
optimization does not preserve a tube or packet as a geometric constraint.

| Measurement at the failed stage | Previous packet, T=.10 | Fresh packet, T=.16 | Tube, T=.16 |
|---|---:|---:|---:|
| Fine-grid H1/2 / initial | 1.099118173 | 1.139234184 | 1.085819500 |
| Fine-grid L3 / initial, 128 samples | 0.976061783 | 0.940664020 | 0.943611420 |
| Fine-grid vorticity / initial, 128 samples | 1.930019689 | 2.029743321 | 2.248627037 |
| Maximum vorticity-ratio discrepancy along paired observations | 16.23% | 20.64% | 19.83% |
| Maximum late stretching discrepancy | 8.14% | 21.65% | 20.94% |
| Fine-grid peak cutoff energy fraction | 0.0000133013 | 0.000175836 | 0.000125190 |

The discrepancies divide by the larger absolute value. The vorticity gate
is 10%; late stretching is also 10%. H1/2 and L3 ratio gates are 2%, and
enstrophy and characteristic-wavenumber ratio gates are 5%. Late H1/2 must
increase by at least 0.1% over the final `0.02` window, with no sampled
reversal; stretching must exceed viscous destruction in that window.
These are triage thresholds, not certified error bounds.

Small cutoff **energy** is insufficient to establish pointwise vorticity
convergence. For the leading field, denser sampling changes the vorticity
measurements by at most 1.224% through `.10`, whereas the `32/64` discrepancy
on the common 128 sampling grid reaches 16.23%. This identifies evolution
resolution as the larger uncertainty. A separate `64/128` check addresses
that uncertainty; it is not supplied by zero padding.

The independent radix-2/FFTW full trajectory replay of the leading field at
`N=64, T=.10` passed 880 shared safe steps, with peak relative state
disagreement `1.77007940722e-15` and peak scaled diagnostic disagreement
`2.47194776478e-15`. See [`independent64`](independent64). The comparison
uses its own safe timestep schedule; agreement between two truncated solvers
does not estimate the missing PDE modes.

## Reproduction and evolved checkpoints

Run from the repository root, using a new output directory:

```sh
python3 research/navier_stokes_cascade/scripts/continue_candidates.py \
  --solver ./build/research/navier_stokes_cascade/navier_stokes_continue \
  --output-dir continuation-replay
```

Each `grid-*` directory retains stage CSVs/logs and `current.chk.gz`, which
contains its final **evolved** Fourier field and restart settings. Initial
coefficient CSVs remain in the preceding campaign. To inspect/continue a saved
state, decompress it to a fresh local path and use `--restart`, an absolute
later `--final-time` on its `.01` clock, and new evidence output. A pair that
failed spatial comparison is not promoted merely because its individual
checkpoint can still evolve.

The initial campaign used the original four-derivative FFTW RK4 storage.
[`solver_source.tar.gz`](solver_source.tar.gz) contains its exact source and
all header dependencies; the hashes match the manifest. The current compact
storage implementation reproduced the entire 1,284-step, `N=32, T=.16`
fresh-packet checkpoint byte for byte; see [`compact_identity.json`](compact_identity.json).
The frozen clock, baselines and BKM telemetry are part of that comparison.
[`checkpoint_payload_audit.json`](checkpoint_payload_audit.json) accounts for
the numerical schema and public initial-field provenance of all six files.

The `N=64, T=.08` replay reproduces the preceding H1/2 result to about
`1e-14` despite its different safe-step count (`645` here versus `667`
previously). The scientific comparison uses fixed physical times. Step
counts alone are not evidence of timestep accuracy.

## A PDE budget for the next search objective

For `A=(-Delta)^(1/2)` and the real Fourier inner product, let

\[
Q=\langle u,Au\rangle=\|u\|_{\dot H^{1/2}}^2,\qquad
F(u)=\mathbb P(u\times\omega)+\nu\Delta u.
\]

The truncated Navier–Stokes equation gives

\[
\dot Q=\underbrace{2\operatorname{Re}\sum_k|k|\,
\overline{\hat u_k}\cdot\widehat{\mathbb P(u\times\omega)}_k}_{P_H}
-\underbrace{2\nu\sum_k |k|^3|\hat u_k|^2}_{D_H},\qquad
\gamma_H=\frac{d\log\|u\|_{\dot H^{1/2}}}{dt}=\frac{P_H-D_H}{2Q}.
\]

[`analyze_critical_budget.py`](analyze_critical_budget.py) evaluates this
identity with NumPy FFTs, separately from the C++ FFTW trajectory. It passes
an exact viscously decaying shear check and agrees with six shared FFTW
diagnostics for all three final `N=64` snapshots to at most `1.38e-14`
relative error. Since `3K<N`, the native-grid cubic mean has no Fourier
alias into the constant mode. It remains floating-point snapshot analysis,
not a third independent trajectory or a continuum error bound.

| Frozen seed | Initial P_H/D_H | Final P_H/D_H | Final gamma_H |
|---|---:|---:|---:|
| Previous packet, T=.10 | 5.69464 | 4.14776 | 0.86706 |
| Fresh packet, T=.16 | 4.85935 | 2.18804 | 0.54099 |
| Tube, T=.16 | 0.87091 | 2.26588 | 0.58954 |

The tube starts with negative critical growth and develops positive growth
later. This argues for a finite-horizon search with late-window objectives,
rather than selecting only a positive instantaneous derivative at time zero.
The final budgets above are still subject to the endpoint resolution limits.

An integral of `gamma_H` alone adds no new objective: it is exactly the log
endpoint H1/2 ratio already used by the optimizer. A genuinely different
trial objective is a smooth minimum of `gamma_H` over a fixed late window,
combined with the existing cutoff, field diversity and holdout checks. It
would penalize late stalls while allowing delayed amplification.

\[
\Phi=-\tau\log\left(\frac1m\sum_{j=1}^{m}e^{-\gamma_H(u(t_j))/\tau}\right),
\qquad w_j=\frac{e^{-\gamma_H(u(t_j))/\tau}}{\sum_i e^{-\gamma_H(u(t_i))/\tau}}.
\]

The local gradient is available from the existing adjoint RHS:

\[
\nabla\gamma_H(u)=\frac{AF(u)+F'(u)^*Au}{Q}
-\frac{2\langle Au,F(u)\rangle}{Q^2}Au.
\]

[`check_critical_rate_gradient.cpp`](check_critical_rate_gradient.cpp) checks
this snapshot-level formula against central differences on `N=16/32`.
Its values are recorded in `critical_rate_gradient.csv`. Integrating these
weighted gradient sources into the reverse trajectory remains the next
implementation step: add `w_j * gradient(gamma_H)` at each fixed observation,
then project the initial gradient onto the existing bandwidth, solenoidal,
reality and fixed-energy constraints. Use a frozen safe timestep schedule
when differentiating and independent adaptive held-out runs afterward;
the snapshot formula does not differentiate adaptive timestep decisions.
The six direction/grid checks attain best relative slope errors below
`1.53e-10` over the three recorded difference spacings. This validates the
local gradient source, not the proposed full trajectory objective yet.
Sustained positive
critical growth is a search heuristic, not a necessary pattern for every
possible singularity mechanism, and none of this constitutes a proof.
