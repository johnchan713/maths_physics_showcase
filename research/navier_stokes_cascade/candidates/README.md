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

Replay it without changing the coefficients:

```bash
./build/research/navier_stokes_cascade/navier_stokes_state_optimize \
  --initial-family wave-packets --seed-bandwidth 3 \
  --state-input research/navier_stokes_cascade/candidates/wave_k3_t008_screening.csv \
  --final-time 0.08 --iterations 0 \
  --output replay.csv --state-output replayed-state.csv
```
