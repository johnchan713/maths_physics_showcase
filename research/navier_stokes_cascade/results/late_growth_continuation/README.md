# Frozen late-growth candidate through T=0.10

The unchanged candidate continues amplifying through both held-out horizons,
`.08` and `.10`. The `32/64` pair passes the declared preliminary spatial
gates, narrowly at `.10`. Fine-grid H1/2 growth reaches **14.55%**, while
its relative growth rate keeps declining and L3 falls. This is evidence
for sustained finite amplification, not a demonstrated finite-time singularity.

## Frozen experiment

[`protocol.json`](protocol.json) was written before the runs. The initial
coefficient file is
[`wave_k3_late_rate_t004_screening.csv`](../../candidates/wave_k3_late_rate_t004_screening.csv),
SHA-256 `893ef61450d670936dc851d28ab3cf7c1aa315c6ff23934646e74137c40ac6ef`.
No coefficients, energy, initial bandwidth, or viscosity were optimized
on these holdouts. Controls are `E(0)=10`, initial bandwidth `3`, `nu=.02`,
evolution cutoffs `10/21`, maximum timestep `.000125`, CFL target `.4`,
and observations every `.01`. Physical samples use grids `64` and `128`.

Both evolution grids start from the same initial coefficients at time zero.
The solver preserves its evolved state, original baselines and observation
clock at every checkpoint. Segments meet exactly at their shared CSV row.
The `.06` results reproduce the earlier optimizer-based validation despite
using a different safe timestep clock.

## Measurements

These are the `N=64` reference values. L3 and vorticity use the common
`128` physical sampling grid, which is **not** a `128` evolution.

| Time | H1/2 gain | L3 change | Sampled vorticity gain | Instantaneous d(log H1/2)/dt |
|---|---:|---:|---:|---:|
| .06 | 8.9597% | -0.1361% | 34.4423% | 1.458609 |
| .08 | 11.9654% | -1.0630% | 62.0566% | 1.255145 |
| .10 | 14.5519% | -2.3379% | 89.4253% | 1.026694 |

The norm rises, but each unit of norm grows more slowly near the endpoint.
All saved critical rates from `.04` to `.10` decline on both grids. The
fine-grid ratio of nonlinear critical-norm production to viscous destruction
falls from `9.53574` at `.04` to `3.68195` at `.10`; it stays above one.
This budget concerns H1/2 squared. It is distinct from the enstrophy
stretching/destruction ratio, which is `5.78215` at `.10`.

The rates come from independent NumPy FFT analysis of actual evolved
checkpoints using `gamma_H = (P_H-D_H)/(2*||u||_Hhalf^2)`. The analysis
passes exact viscous-shear decay and compares its energy, enstrophy, H1/2,
stretching and dissipation against each saved C++ observation. No finite-time
singularity law is fitted. Monotonically increasing gamma is a search
preference, not a necessary condition for every possible blow-up mechanism.

## Resolution and numerical limits

[`assessment.json`](assessment.json) compares every common saved time,
not only endpoints. The percentages below are relative discrepancies,
divided by the larger absolute value.

| Maximum discrepancy through horizon | T=.08 | T=.10 | Gate |
|---|---:|---:|---:|
| H1/2 ratio | .00742% | .02404% | 2% |
| L3 ratio | .04741% | .05899% | 2% |
| Enstrophy ratio | .21480% | .89277% | 5% |
| Characteristic wavenumber ratio | .10780% | .44984% | 5% |
| Sampled vorticity ratio | 2.48044% | 9.08869% | 10% |
| Stretching in the final .02 window | 3.23934% | 9.42021% | 10% |

At `.10`, the peak per-step cutoff energy fraction is `.00212943` on
`N=32` and `.0000321403` on `N=64`: **.212943%** and **.00321403%**,
respectively, against the `.8%` gate. The largest change in raw sampled
vorticity from `64` to `128` spatial sampling is `1.34523%`, below the
2% sampling gate. The dimensionless critical-rate gap
`T*max|gamma_32-gamma_64|` is `.00087748` at `.10`, below `.002`.

The passing margins for vorticity and stretching are small. They are
screening limits, not 10% error bounds for the PDE. This campaign does not
include `128` dynamics; that is the next required comparison before trusting
pointwise growth at a longer horizon. Small retained-edge energy does not
bound the absent Fourier tail.

## Independent and perturbed checks

The separate `refined64` run starts again from the original initial field,
halves the maximum timestep to `.0000625`, and halves the observation
interval to `.005`. Comparing common physical times tests both temporal
integration and sensitivity to extra observation boundaries; the added
observations also check for a missed late norm reversal. This is a combined
perturbation, not two separately isolated error estimates. Exact critical
budget samples remain on the recorded `.01` clock.

The reference/refined step counts rise from `662` to `1,284` at `.08`
and from `880` to `1,608` at `.10`, confirming actual refinement. Maximum
relative changes over common observations are `1.51e-13` for the H1/2
ratio and `5.09e-12` for the sampled vorticity ratio. No sampled late H1/2
reversal appears. These tiny temporal changes do not reduce the much
larger spatial-resolution uncertainty reported above.

`independent64` compares the complete radix-2 and FFTW trajectories through
`.10`, starting from the same frozen initial data. Its manifest, CSV and log
record the actual outcome. Agreement concerns two implementations of the
same truncated equation, not a rigorous continuum error estimate.

The independent replay passes **945 shared safe steps**, with peak
whole-state relative disagreement `2.45875287219e-15` and peak scaled
diagnostic disagreement `2.23401474344e-15`. Its own shared safe clock
differs from the continuation clock, hence 945 versus 880 reference steps.
The snapshot budget crosschecks have maximum relative diagnostic
disagreement `9.12e-14` over all 21 saved reference/refined snapshots.

Final completeness and validation are checked by [`audit.py`](audit.py).
Only a successfully produced `audit.json` with status
`passed-declared-finite-amplification-gates` establishes that all required
legs have completed and their archived identities and gates passed.
The completed archive passes this audit, including all six restored
checkpoint identities and reproduction of the previous `.06` H1/2 ratios.

## Reproduction and saved states

From the repository root, run one grid into a new output directory:

```sh
python3 research/navier_stokes_cascade/results/late_growth_continuation/run.py \
  --solver ./build/research/navier_stokes_cascade/navier_stokes_continue \
  --grid 64 --output-dir new-continuation/reference64
```

`--solver` selects the built continuation executable. `--grid` selects
the evolution resolution; repeat with `32` for the other primary leg.
`--output-dir` must be new so evidence cannot be overwritten. The frozen
field, stages, viscosity and gates are recorded by the driver. For the
perturbed leg, add `--max-dt .0000625 --observation-interval .005` and
use a new `refined64` directory. `assess.py --root new-continuation`
reports whichever legs have completed and never counts missing legs as passes.

Each trajectory directory preserves `t080.chk.gz` and `t100.chk.gz`.
These contain evolved Fourier states, not optimizer initial coefficients.
Decompress a chosen file to a fresh local path, then pass it through
`--restart` with a later absolute `--final-time` on its frozen observation
clock and fresh evidence output. All scientific settings and normalization
are restored from the checkpoint. The archive audit verifies compressed and
restored SHA-256 identities and matches each time and grid to its manifest.

## Research decision

Keep this field as a finite-amplification lead. The held-out growth is
stronger than at `.06`, but the rate declines and L3 decreases. There is
no demonstrated accelerating divergence, no rigorous tail estimate and no
proof or disproof of regularity. The next bounded experiment should evolve
these same initial coefficients on `N=128` from time zero through `.10`,
compare `64/128` on a common physical grid, and only then consider `.12`.
Padding the evolved `N=64` endpoint cannot supply the missing finer evolution.
