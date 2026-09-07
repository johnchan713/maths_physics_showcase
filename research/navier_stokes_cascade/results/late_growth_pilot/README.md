# Matched late-growth / endpoint search pilot

The new objective gives small improvements in all four matched discovery
starts. Neither low-grid search arm passes its full held-out screen, but
the frozen leading late-mode field passes a **separate** `32/64` validation
through `.06`, including timestep, denser sampling and independent FFTW
checks. Its H1/2 gain is about **8.96%**. This is finite amplification,
not a proof of blow-up or regularity.

## Protocol frozen before the pilot

This is a bounded comparison of candidate-search objectives, not a blow-up
test or a certified bound for the PDE. The new objective replaces endpoint
critical-norm gain with `T * softmin(gamma_H)`; characteristic-scale reward
(`0.15`) and terminal cutoff cost (`0.04`, threshold `.01`) are unchanged.
Both searches use the amplification track, without shape penalties.

- Four matched starts: the built-in wave packet, vortex pair, orthogonal
  bundle, and the frozen previously optimized packet
  `candidates/wave_k3_amplification_t004_screening.csv` (SHA-256
  `5f3dc4449abd58f5ae5e1aca3443b2a0fa7ad06029f9dadae939ff2b9a3b51d8`).
- Initial Fourier bandwidth `3`, initial energy `10`, viscosity `.02`.
- Discovery on evolution grids `16/32`, fixed differentiated coarse step
  `.0005`, adaptive fine maximum `.0005`, horizon `.04`, two ascent steps
  maximum per start. Every step requires the geodesic finite-difference gate.
- Late window `[.5T,T]`, five equally spaced observations, temperature `.1`
  in inverse-time units. This clock is a subset of the historical eight
  profile/evidence times, so matched discovery runs have the same step
  boundaries. No post-result window or temperature tuning.
- At most two nondominated, numerically usable finalists per method, with
  the existing family-diversity rule. The late method gives first shortlist
  priority to its conservative minimum rate; the endpoint method retains
  its historical ordering. Thus holdout shortlists are not necessarily matched.
- Frozen selected fields are replayed through held-out time `.06`; then
  with materially finer timesteps and doubled spatial/time sampling. The
  late window is `[.03,.06]` there; no coefficients are optimized on it.
- All existing robust-search cutoff (`.008`), common-grid convergence,
  stretching, late-growth, gain/spread, and independent FFTW gates stay on.
  The late method additionally requires every sampled late rate positive,
  `T * max_j |gamma_coarse-gamma_fine| <= .002`, and changes of its minimum
  under timestep/sampling perturbations of at most `.001/T`.
- A positive normalized soft minimum alone is insufficient: it is greater
  than or equal to the actual minimum. Denser checks can reveal a missed
  reversal. None of these finite samples controls every point or time.

The two `manifest.json` files retain commands, hashes, complete traces,
rejected runs, finalist decisions, and raw physical evidence. Equal numbers
of ascent iterations are a matched iteration budget, not equal CPU cost:
the new objective needs extra RHS and adjoint-source evaluations.

## Mathematical objective

With the real Fourier inner product, set `A=(-Delta)^(1/2)`,
`Q=<u,Au>=||u||_(H1/2)^2`, and let `F(u)` be the retained Navier--Stokes RHS
including viscosity. Then

```text
gamma_H(u) = d log(||u||_(H1/2))/dt = <Au,F(u)> / Q,
grad gamma_H = (A F(u) + F'(u)^* A u) / Q - 2 gamma_H A u / Q,
Phi = -tau log((1/m) sum_j exp(-gamma_H(u(t_j))/tau)),
J = T Phi + .15 log(k_rms(T)/k_rms(0)) - .04 cutoff_cost(T).
```

The adjoint source at each observation is `T*w_j*grad gamma_H`, where
`w_j=exp(-gamma_j/tau)/sum exp(-gamma/tau)`. Insert it before reversing the
step that ends at that observation, including the terminal source once.
The initial gradient is projected onto the real, solenoidal, band-limited,
fixed-energy constraint. Timesteps are frozen while differentiating; trial
states violating the conservative timestep bound are rejected. Independent
adaptive replays are validation, not derivatives of adaptive decisions.

This score can distinguish equal endpoint gains with different late growth.
Simply integrating `gamma_H` would reproduce the old endpoint log ratio.
The objective is our experimental choice, not a theorem or a claim from
the literature. The general use of adjoint searches for extreme finite
Navier--Stokes amplification has primary precedent in
[Kang, Yun and Protas](https://arxiv.org/abs/1909.00041), whose tested flows
also showed finite amplification without establishing a singularity.

## Outcome

Both four-start discovery arms completed, with all 16 gradient checks below
`1.255e-5` relative discrepancy (gate `.01`). Each accepted two ascent
iterations per start. Replaying the endpoint-selected fields under the new
score, with **zero** optimization steps, preserves their coefficient CSVs
and every exported physical observation byte for byte. See
[`comparison.json`](comparison.json) and the four `endpoint_late_replays`.

At discovery time `.04`, conservative values take the smaller result of
the two evolution grids:

| Matched initial seed | Endpoint H1/2 gain | Late-mode H1/2 gain | Endpoint weakest late rate | Late-mode weakest late rate | Rate improvement |
|---|---:|---:|---:|---:|---:|
| Fresh packet | 4.4153% | 4.5085% | 1.128686 | 1.157010 | 2.509% |
| Vortex pair | 1.5904% | 1.6340% | .395837 | .407011 | 2.823% |
| Orthogonal bundle | 2.9712% | 3.0286% | .752726 | .769258 | 2.196% |
| Previous optimized packet | 5.5889% | 5.6578% | 1.415650 | 1.440033 | 1.722% |

Rates are in inverse time. All four late-mode states also improve the
critical-norm ratio over the last half of discovery, but improvements are
small. Four starts and two iterations do not establish statistical or global
superiority of the new optimizer. No window/temperature tuning was done.

Both methods shortlisted the previous packet and the bundle. **Neither
arm passes all low-grid held-out gates.** At `T=.06`, the new packet exceeds
the coarse hard cutoff gate; the endpoint packet reaches a coarse cutoff
fraction `.0089406`, beyond the stricter `.008` screening margin. The bundle
reference step is too large under the fixed-step bound. Halving it completes
the bundle replay, but the late-mode bundle also fails the new cross-rate
agreement gate. Those failures are retained, not replaced by successful runs.

The strongest new field is frozen at
[`late/checkpoint-discovery/state.csv`](late/checkpoint-discovery/state.csv).
`validate_leader.py` performs a **separate** `32/64` resolution follow-up at
the same held-out `.06`, with zero ascent iterations and a maximum step
`.000125`, followed by actual timestep refinement, doubled sampling, and
an independent `64`-grid FFTW trajectory if the prior gates pass. Its
[`leader_32_64/manifest.json`](leader_32_64/manifest.json) now records
`validated-finite-amplification-shortlist`, with no failed follow-up gates.
That label does not retroactively pass the original low-grid experiment.

### Separate higher-resolution result

The completed `32/64` reference replay at `.06` passes its preliminary
spatial, cutoff and sampled late-growth gates:

| Quantity at T=.06 | N=32, cutoff 10 | N=64, cutoff 21 |
|---|---:|---:|
| H1/2 / initial | 1.089587921 | 1.089596688 |
| L3 / initial, common 64 samples | .998517288 | .998639098 |
| Vorticity / initial, common 64 samples | 1.369796764 | 1.344960506 |
| Characteristic wavenumber / initial | 1.237828697 | 1.237956734 |
| H1/2 gain during the last half | 4.724873% | 4.725712% |
| Minimum sampled late critical rate | 1.457503198 | 1.458609239 |
| Peak per-step cutoff energy fraction | 7.721435e-5 | 1.742562e-8 |
| Accepted reference steps | 480 | 493 |
| Accepted refined-timestep steps | 960 | 992 |

The smaller-step replay also passes. Its final H1/2 ratios change by at
most `1.364e-13` relatively; sampled final vorticity ratios change by at
most `2.689e-12`. This is an actual timestep refinement, not merely a
smaller unused adaptive maximum. It does not establish continuum accuracy.

The critical logarithmic growth rate stays positive at the saved late
times, but falls from about `1.570` at `.0375` to `1.459` at `.06` on the
fine grid. Thus the critical norm keeps growing while its relative growth
rate already slows. This remains a finite-amplification lead, not evidence
of a divergent growth law.

The denser check doubles the spatial sampling grid to `128` and late-rate
intervals from four to eight. The fine vorticity amplification ratio changes
by `0.03996%` to `1.344423140`; the coarse change is `0.17095%`. The actual
minimum late rate remains at the final time and remains positive. No
sampled late reversal appears. These are samples of `32/64` evolutions,
not new `128`-grid dynamics or certified continuous maxima.

The independent internal-FFT/FFTW comparison passes **491 shared safe
steps**, with peak whole-state disagreement `1.41234478668e-15` and peak
scaled diagnostic disagreement `1.42675977754e-15`. See
[`fftw.log`](leader_32_64/fftw.log) and
[`fftw.csv`](leader_32_64/fftw.csv). Its shared clock differs slightly from
the fixed observation clock, explaining 491 versus 493 reference steps.
Agreement concerns the same truncated equation, not the uncomputed tail.

Across both evolution grids, timestep refinement and denser sampling,
the smallest H1/2 gain is `.08958792145`, while the observed spread of
final ratios is `8.766322e-6`. The gain/spread gate passes. The largest
reference vorticity-ratio disagreement between grids is `1.8131%` at the
endpoint; this is empirical convergence evidence, not a PDE error bound.

## Reproduction and validation

Use `scripts/robust_search.py` with the protocol above and fresh output
directories; select `--growth-objective endpoint` and `late-rate` for the
two arms. `compare.py` scores the endpoint arm under the late objective
without modifying it. `record_validation.py` captures actual local test
logs and exit codes. The new CLI test checks evidence completeness in
addition to exit status; any saved failure remains evidence of a failed
check until a diagnosed correction and successful rerun are recorded.

The derivative suite covers stable soft-minimum weights, negative true
minima behind positive smooth scores, exact viscous shear decay, bad
options and data, source timing (including the final observation), missing
trajectory data, and full fixed-energy geodesic finite differences in
three directions on both `N=8` and `N=16`. The six best relative errors
are at most `2.148e-8`. Fixed differentiated schedules are not gradients
of adaptive-step branch decisions.

The initial validation directory preserves an intermittent empty-evidence
observation in this runner. The final code explicitly closes and checks
each output before reporting success, and the CLI test requires every
evidence row before proceeding. The final local twelve-check pass is saved
separately in `validation_final`; the precise cause of the earlier empty
file was not established. Clean Ubuntu CI independently exercises the same
completeness checks. [The core-code workflow](https://github.com/johnchan713/maths_physics_showcase/actions/runs/34096901426)
passed all 27 steps, including ASan/UBSan and the new full-gradient and
CLI suites. No failed run is silently counted as successful.

The pilot and first `32/64` replay used the pre-I/O-hardening executable.
`pre_io_source` retains the two subsequently changed source files, with
hashes matching the original discovery manifests; the other mathematical
sources are unchanged. The final code also avoids storing a reverse
trajectory when `--iterations 0` makes it unnecessary. Frozen-state and
physical-evidence bytes match in the short old/new executable replay.
The higher-grid driver originally mistook changed CSV grid metadata for
changed coefficients. Its correction compares all 342 physical coefficient
rows and constraints, excluding only the two writing-grid metadata columns;
the first manifest records the postprocessing failure and recovered status.

[`audit.json`](audit.json) inventories the final archived files and verifies
recorded hashes, source provenance, completed regressions and the final
follow-up status. The canonical initial-field copy is
[`wave_k3_late_rate_t004_screening.csv`](../../candidates/wave_k3_late_rate_t004_screening.csv),
SHA-256 `893ef61450d670936dc851d28ab3cf7c1aa315c6ff23934646e74137c40ac6ef`.
These optimizer `state.csv` files contain **initial coefficients**, not
saved evolved states.

Finally, a fixed Fourier cutoff cannot literally produce the desired PDE
blow-up while respecting the unforced energy inequality. On the retained
cube, `|k| <= sqrt(3)*K`, so
`||u||_(H1/2)^2 <= 2*sqrt(3)*K*E(0)` for the exact Galerkin solution.
The search can identify finite amplification and shapes worth resolving,
but a proof would still require a limiting argument and rigorous control
of the omitted Fourier modes. A better score does not remove that gap.

## Next experiment

Freeze this initial field and use the independent FFTW checkpointed
continuation path for held-out `.08` and `.10` horizons, raising evolution
resolution only as the common-grid pointwise and cutoff tests require.
Do not optimize on those holdouts. Track whether the positive but declining
critical rate recovers or approaches a plateau. Additional starts and
temperatures would be a new, separately declared search campaign, not
post-hoc tuning of this four-start comparison.
