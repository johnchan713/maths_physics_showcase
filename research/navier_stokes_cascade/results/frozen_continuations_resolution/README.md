# Leading-field resolution check through T=0.10

The frozen leading packet passes the predeclared **preliminary finite-
amplification** comparison after raising the evolution grids from `32/64`
to `64/128`. The vorticity discrepancy falls from 16.23% to 4.63%. This is
stronger numerical evidence for finite amplification; the remaining
pointwise discrepancy is material and no continuum singularity is established.

The field is unchanged from the earlier search (`E(0)=10`, initial `K=3`,
SHA-256 `5f3dc4449abd58f5ae5e1aca3443b2a0fa7ad06029f9dadae939ff2b9a3b51d8`).
Viscosity stays `0.02`. Both evolutions use maximum timestep `0.000125`,
CFL target `0.4`, and observations every `.01`. The finer trajectory starts
from the same initial coefficients at **time zero**, not by padding an
already evolved coarse field.

## Measurements

[`assessment.json`](assessment.json) compares every common observation.
[`fine128_evidence.csv`](fine128_evidence.csv) is the complete fine trajectory;
the coarse data are in the preceding frozen-field campaign.

| At T=.10 | N=64, cutoff 21 | N=128, cutoff 42 |
|---|---:|---:|
| H1/2 / initial | 1.099118173 | 1.099119723 |
| L3 / initial, common 128 samples | 0.976061783 | 0.976075429 |
| Vorticity / initial, common 128 samples | 1.930019689 | 2.023656715 |
| Enstrophy / initial | 1.682256114 | 1.682498648 |
| Characteristic wavenumber / initial | 1.313185871 | 1.313280739 |
| Stretching / viscous enstrophy destruction | 6.707380621 | 6.652871114 |
| Peak cutoff energy fraction, every accepted step | 1.33013009e-5 | 7.29780316e-9 |
| Accepted steps | 834 | 1433 |

Maximum relative discrepancies over the saved trajectory are:

- H1/2 ratio: `0.000141%`;
- L3 ratio: `0.001398%`;
- enstrophy ratio: `0.014415%`;
- characteristic-wavenumber ratio: `0.007224%`;
- sampled vorticity ratio: `4.62712%`;
- late stretching: `0.455441%`.

H1/2 increases another `1.8563%` over the last `.02` window on the fine
grid, with no sampled late reversal. L3 still decreases. The 10% vorticity
gate is an empirical screening threshold, not a 10% error bound for the PDE.
The small H1/2 disagreement does not imply equally small vorticity error.

## Independent and perturbed checks

The full independent radix-2/FFTW replay at `N=64, T=.10` passes 880 shared
safe steps. Its largest whole-state discrepancy is `1.77007940722e-15`.
The `oracle64` files here are the same run also archived with the preceding
campaign; they do not represent an additional independent experiment.

The separate `N=64` timestep rerun halves the maximum timestep and raises
the accepted step count from 834 to 1,604 (`1.92326x`). Maximum relative
changes along the matched observation times are `7.26e-14` for the H1/2
ratio and `6.46e-12` for the sampled-vorticity ratio. See
[`timestep64/assessment.json`](timestep64/assessment.json). This timestep
perturbation and the independent full trajectory are at `N=64`; no separate
full `N=128` timestep study or radix-2 trajectory is claimed.

[`sample_checkpoint_maximum.py`](sample_checkpoint_maximum.py) uses NumPy
real FFTs to sample the **same** saved field on nested `128/256` grids.
At `.09`, the final maximum increases by `0.07355%`. At `.10`, the two
sampled final maxima coincide to roundoff (`141.4498486`), while the denser
initial maximum is slightly larger. Consequently the final amplification
ratio changes from `2.023657` to `2.022529` (`0.05574%`). The native NumPy
maximum matches FFTW to approximately `2e-16` relative error. Equal maxima
on nested grids do not certify the continuous maximum or control all times.

## Restart and provenance

The first `N=128` segment was deliberately interrupted after the `.03`
checkpoint when the independent `N=64` replay had finished. Compact RK4
storage reproduced that complete `.03` checkpoint byte for byte. The final
segment resumed it through `.10`, with the same field, viscosity, cutoff,
timestep controls and observation clock. Both CSV segments agree exactly
at the shared boundary. [`provenance.json`](provenance.json) identifies the
interruption, implementation versions and successful final completion.
The compact implementation also reproduced a full 1,284-step `N=32`
research checkpoint exactly in the earlier campaign.

The final 128-grid checkpoint contains `100,663,792` uncompressed bytes.
Its gzip stream is stored in ordered parts with individual and combined
hashes in [`fine128_checkpoint.json`](fine128_checkpoint.json). Restore it
from a repository checkout to a new local filename:

```sh
python3 research/navier_stokes_cascade/results/frozen_continuations_resolution/restore_checkpoint.py \
  --output fine128-restored.chk
```

The restored bytes have been compared with the actual solver output, and
existing output files are protected. The uncompressed SHA-256 is
`1b3f01637ce4686df5d10c02119ee947cb4442a10064fbf9cab399fec11d2a9b`.
The checkpoint preserves the evolving state and its original normalization;
it can resume without repeating the expensive start. The single-file
`timestep64/final.chk.gz` retains the timestep-refined state as well.

All archive contents are generated numerical metadata and Fourier doubles
derived from the public initial field and solver source. Their schema and
identities are checked in [`checkpoint_audit.json`](checkpoint_audit.json).
[`manifest.json`](manifest.json) records the archived file hashes.

To reproduce the fine trajectory directly with the current implementation,
use the same frozen CSV and the options above with `--grid 128`,
`--sampling-grid 128`, `--dense-sampling-grid 128`, `--final-time .10`, and
fresh output/checkpoint paths. `commands.json` and the resume/timestep
records retain the exact commands actually executed in this session.

## Next experiment

The next candidate-search change is to insert the checked critical-growth
gradient into the existing reverse trajectory and optimize a smooth minimum
over a late time window. The [budget derivation and local gradient checks](../frozen_continuations/README.md#a-pde-budget-for-the-next-search-objective)
are ready; a full optimizer for that objective has not yet been implemented.
Keep independent held-out horizons and resolution gates, preserve delayed-
growth seeds, and recheck pointwise convergence for any stronger finalist.
More finite amplification, by itself, still cannot supply the unresolved
Fourier-tail bounds or limiting argument needed for a proof.
