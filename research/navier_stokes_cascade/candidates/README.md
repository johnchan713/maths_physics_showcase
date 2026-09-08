# Screening candidates

These full-precision Fourier coefficient files preserve reproducible search
checkpoints. They are numerical screening artifacts, not certified solutions,
proof data, or candidates that passed every promotion gate.

`wave_k3_late_rate_t004_screening.csv` is the leading frozen initial field
from the matched late-growth pilot. Two fixed-energy ascent steps optimize
the weakest sampled late critical-growth rate through `T=.04`. It reaches
about 8.96% H1/2 growth in the separate `N=32/64, T=.06` replay, passing
timestep, denser-sampling and independent FFTW checks; L3 remains
below its initial value. The original `16/32` held-out screen failed its
cutoff gate, and no singularity is established. See the
[complete pilot and validation status](../results/late_growth_pilot/README.md).
The SHA-256 digest of this coefficient file is
`893ef61450d670936dc851d28ab3cf7c1aa315c6ff23934646e74137c40ac6ef`.

The unchanged field subsequently reaches 14.55% H1/2 growth through
`T=.10` on `N=64`; L3 decreases by 2.34% and relative H1/2 growth slows.
The `32/64` pointwise and stretching gates pass narrowly. See the
[frozen continuation and validation record](../results/late_growth_continuation/README.md)
for the numerical limitations, checks and evolved restart states. The
[completed 128-grid replay](../results/late_growth_resolution128/README.md)
reproduces `14.55%` H1/2 growth through `.10` and passes the retained comparison
and endpoint sampling gates. Vorticity reaches `1.982` times its initial value
on `256` physical samples, but differs by `4.47%` from the 64-grid evolution
on those same samples. L3 and the relative critical growth rate still decline.

The [unchanged-field extension to .12](../results/late_growth_t012/README.md)
also passes the preliminary gates. H1/2 grows `16.65%` and densely sampled
vorticity reaches `2.290` times its initial value, while L3 falls `3.77%` and
relative H1/2 growth slows further to `.792`. The archived comparison records
both the worst pointwise gap at `.11` and the closer endpoint agreement at
`.12`. No new initial coefficient candidate is introduced.

The [extension to .14](../results/late_growth_t014/README.md) passes the
declared checks at `.13` but fails the `64/128` vorticity agreement gate at
`.14`: the discrepancy is `12.39%`, and remains `12.65%` with denser physical
sampling, above the fixed `10%` limit. Fine H1/2 growth reaches `18.24%`, but
its relative rate falls to `.557` and L3 decreases `5.23%`. The initial field
remains unchanged. Further time promotion of this grid pair is paused pending
finer evolution checks at the same horizon; no converged local maximum or
singularity claim follows from the larger fine-grid vorticity.

The [subsequent snapshot diagnosis and cost probe](../results/late_growth_peak_audit/README.md)
retain that failure and introduce no new coefficients. The maximum sampled
vorticity-vector difference is `21.34%` of the fine peak, larger than the
scalar-peak gap. Modes absent from the coarse grid have a disproportionate
effect on vorticity despite their small energy, and shared modes also differ.
The short 128/256 capacity test passes; it is not a full `.14` replay. Its
measured cost makes a validated performance improvement the next practical
step before higher-resolution continuation.

For this field, replay with `--search-track amplification --growth-objective
late-rate --grid 32 --fine-grid 64 --seed-bandwidth 3 --iterations 0
--final-time .06 --dt .000125 --fine-max-dt .000125`, and distinct fresh
output paths. Raising the writing grid changes the CSV's `source_grid`
and `simulation_cutoff` metadata, not the initial Fourier coefficients.

`wave_k3_t008_screening.csv` is the best paired state from the first
684-variable, fixed-energy adjoint pilot. It passed timestep, cutoff, and
coarse/fine agreement checks, but its rescaled-profile drift remained about
6.3 against the required threshold of 1. The main research README records the
complete interpretation and measurements. Its SHA-256 digest is
`c94bdc0d38f1bd52fcca2072e9dcbe9103a4421d50fa9f10ea84541a4c99a128`.

`wave_k2_profile_t004_screening.csv` is the lowest-drift state retained from
the first deterministic multi-start, smooth-profile adjoint pilot and its
profile-dominant continuation. At `T=0.04`, the coarse/fine H1/2 ratios are
`1.013871` and `1.013876`, while the strict fixed-threshold profile drifts are
`5.795320` and `5.838034`. Its fine cutoff-shell fraction is about
`2.02e-8`, and the coarse/fine pair agrees, but the drift remains almost six
times the promotion threshold and L3 decreases. Its SHA-256 digest is
`dbcdbe925512e7e4896dc3b972a491e27e5b9454e38a4f86f515718f49f4d9df`.

`wave_k3_path_t008_screening.csv` is the best state from the first
path-dependent, four-snapshot adjoint search. At `T=0.08`, its coarse/fine
H1/2 ratios are `1.054411` and `1.055435`, and its strict profile drifts are
`4.629944` and `4.605072`. The coarse cutoff fraction is `0.00997863`, only
just inside the one-percent gate, while the fine value is `0.000122778`.
Coarse/fine vorticity amplification differs by about 14%, so the candidate
fails the 10% agreement gate and is not eligible for refinement. An exact
independent FFTW trajectory replay agreed with the internal state to
`8.27e-16`. Its SHA-256 digest is
`cd73fe95e4a8c69f193dff3ebe5abd07a8716bf3739ae412a30b95b008b59491`.

`wave_k3_smoothmax_t008_screening.csv` is the best paired local control from
the eight-snapshot smooth-maximum experiment around the preceding checkpoint.
At `T=0.08`, its coarse/fine H1/2 ratios are `1.054363` and `1.055388`; its
strict drifts are `4.606862` and `4.574379`. The worst drift and vorticity
agreement improve slightly, and coarse cutoff loading falls to `0.00997425`,
but fine-grid drift worsens and vorticity disagreement remains `14.0576%`.
It is retained for exact follow-up, not promoted. An independent FFTW replay
agreed to `7.85e-16`. Its SHA-256 digest is
`75818222d22aaaeb0b92f1a417a0309400e810335712b1df9aef729884fcc7b1`.

Replay the `K_seed=3` checkpoint without changing the coefficients:

```bash
./build/research/navier_stokes_cascade/navier_stokes_state_optimize \
  --initial-family wave-packets --seed-bandwidth 3 \
  --state-input research/navier_stokes_cascade/candidates/wave_k3_t008_screening.csv \
  --final-time 0.08 --iterations 0 \
  --output replay.csv --state-output replayed-state.csv
```

Replay the profile-oriented checkpoint with:

```bash
./build/research/navier_stokes_cascade/navier_stokes_state_optimize \
  --initial-family wave-packets --seed-bandwidth 2 \
  --profile-shape-weight 1 \
  --state-input research/navier_stokes_cascade/candidates/wave_k2_profile_t004_screening.csv \
  --final-time 0.04 --iterations 0 \
  --output profile-replay.csv --state-output profile-replayed-state.csv
```

Replay the path-dependent checkpoint through the optimizer with:

```bash
./build/research/navier_stokes_cascade/navier_stokes_state_optimize \
  --initial-family wave-packets --seed-bandwidth 3 \
  --profile-shape-weight 0.25 --profile-path-weight 1 \
  --profile-path-aggregation mean --profile-path-samples 4 \
  --state-input research/navier_stokes_cascade/candidates/wave_k3_path_t008_screening.csv \
  --final-time 0.08 --iterations 0 \
  --output path-replay.csv --state-output path-replayed-state.csv
```

Replay the smooth-maximum local checkpoint with:

```bash
./build/research/navier_stokes_cascade/navier_stokes_state_optimize \
  --initial-family wave-packets --seed-bandwidth 3 \
  --profile-shape-weight 0.25 --profile-path-weight 1 \
  --profile-path-samples 8 --profile-path-temperature 0.01 \
  --state-input research/navier_stokes_cascade/candidates/wave_k3_smoothmax_t008_screening.csv \
  --final-time 0.08 --iterations 0 --diagnostic-every 1 \
  --output smoothmax-replay.csv --state-output smoothmax-replayed-state.csv
```

Replay the exact same coefficients through both independent evolution paths:

```bash
./build/research/navier_stokes_cascade/navier_stokes_fftw_compare \
  --grid 32 --viscosity 0.02 --energy 10 --dt 0.005 \
  --final-time 0.08 --diagnostic-every 25 \
  --state-input research/navier_stokes_cascade/candidates/wave_k3_smoothmax_t008_screening.csv \
  --output smoothmax-fftw-comparison.csv
```
