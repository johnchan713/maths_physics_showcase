# Screening candidates

These full-precision Fourier coefficient files preserve reproducible search
checkpoints. They are numerical screening artifacts, not certified solutions,
proof data, or candidates that passed every promotion gate.

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
