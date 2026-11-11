# Q182: fixed-gauge comparison and protected-core replay

Parent: Q181, `c82212fb774e7a83b460c9ee90fb5b3ce761c9e9`.
The current remote head, Q181, Q163 and the corrected Q169/Q170 comparison and
core audits, their frozen evidence and current CI were read first. Historical
evidence, the inherited physical field, compact force, pressure gauge,
localization plateau, supports and `C^-20` width remain unchanged.

Q181 closes the physical diagonal, endpoint and localization/force chain
relative to its named continuum providers and leaves fixed-gauge comparison
and the protected-core limit next. Q163 contains this implication, while
Q169/Q170 supply the later source-level corrections used here. The conclusion
is conditional on those continuum estimates. Tests and hashes verify exact
arithmetic, interfaces and provenance; they are not a continuum proof.

## 1. Same field, force and pressure gauge

Q181 and Q163 contain the byte-identical Q162 evidence. Q181 and Q169 contain
the byte-identical Q168 evidence, and Q170 consumes the frozen Q169 record.
Thus the comparison uses Q181's exact localized field, compact smooth force,
fixed-gauge pressure, annular support hole and unit localization plateau. No
parallel field or changed pressure normalization is inserted.

Fix `T<1`. Let `(u,p)` be the Q181 field and pressure, and let `(v,P)` be a
locally smooth competitor on `[0,T]` with the same force, zero datum and
uniform `L2` bound. Put

\[
 w=v-u,\qquad \pi=P-p,
 \qquad g=w\otimes w+w\otimes u+u\otimes w.          \tag{1}
\]

The individual pressure equations retain the force divergence:

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f,
 \qquad
 -\Delta P=\partial_i\partial_j(v_i v_j)-\nabla\cdot f. \tag{2}
\]

Subtracting (2) cancels the common force exactly once:

\[
 w_t+\operatorname{div}g=\Delta w-\nabla\pi,
 \qquad -\Delta\pi=\partial_i\partial_jg_{ij}.       \tag{3}
\]

Neither individual pressure is replaced by an unforced formula, and the
force is not subtracted twice.

## 2. Pressure recovery needs momentum and finite energy

Set

\[
 \pi_*=\sum_{i,j}R_iR_jg_{ij}.                       \tag{4}
\]

The `L2` hypothesis makes `g` uniformly `L1`, so `pi_*` lies in `H^-2`.
Poisson's equation alone only says that `pi-pi_*` is harmonic. For each
compact time test, the momentum equation places its gradient in `H^-3`. A
harmonic `H^-3` distribution has a weighted-`L2` Fourier transform supported
at the single point zero and therefore vanishes. Hence

\[
 \nabla\pi=\nabla\pi_*                              \tag{5}
\]

as a space-time distribution. A time-dependent spatial pressure constant may
remain, but its cutoff flux vanishes by `div(w)=0`. No competitor pressure
decay is assumed.

The finite-energy hypothesis is essential. Without it,
`v=t e_1, P=-x_1` and `u=p=0` have the same zero datum and force but a nonzero
harmonic pressure gradient.

## 3. Expanding-ball estimate

Split `g=w tensor w+h`, where `h=w tensor u+u tensor w`. The Q181 velocity has
fixed compact support, so its mixed pressure has pointwise annular decay
`R^-3` and flux decay `R^-5/2`.

For `chi_R=phi(x/R)^8`, the quadratic pressure commutator at exponent `4/3`
uses

\[
 R^{-4/3}\int_0^Rr^{-2/3}\,dr
 +\int_R^\infty r^{-2}\,dr={4\over R}.              \tag{6}
\]

Thus the commutator norm is `O(R^-3/4)`. The eighth cutoff power is the first
one satisfying the required weight comparison; the seventh fails. After
local Sobolev and Young inequalities, the two nonlinear rows have decay
`R^-4` and `R^-14/5`, both with `R^-5/2` Sobolev offsets. The coarse
cutoff-Laplacian term is limiting:

\[
 {1\over2}E_R'+{1\over2}D_R^2
 \le\|\nabla u(t)\|_\infty E_R+C_TR^{-2},
 \qquad E_R(0)=0.                                   \tag{7}
\]

No competitor derivative, dissipation or boundary-flux estimate enters (7).
For fixed `T<1`, Gronwall gives

\[
 E_R(t)\le2C_TT e^{2TL_T}R^{-2}.                    \tag{8}
\]

Sending `R` to infinity at this fixed horizon yields `v=u` on `[0,T]`.
Neither `C_T` nor `L_T` needs to remain bounded as `T` approaches one.

## 4. Protected physical core

Keep

\[
 h={1\over200},\qquad A={101\over200},\qquad
 X_{in}={1\over101}<X_a={4\over101}.                 \tag{9}
\]

For `tau=1-t`, choose

\[
 x_\tau=(\sqrt{2X_{in}\tau},0,0).                   \tag{10}
\]

On this path `q=tau`, `eta=0` and `X=X_in`. Every annular diagonal correction
has closed support beginning at `X_a`, so it vanishes identically near the
path. The shared positive-order Borel background does not vanish; it satisfies

\[
 \tau^A u_\theta(x_\tau,1-\tau)=e_0+R_{bg}(\tau),
 \qquad |R_{bg}(\tau)|\le C_{bg}\tau^{2h},
 \qquad e_0>0.                                      \tag{11}
\]

Let

\[
 \rho=\min\{1,e_0/(2C_{bg})\},\qquad
 q_{core}=\min\{q_0,\rho^{100}\}.                  \tag{12}
\]

The exponent `100=1/(2h)` is exact. After also imposing the time and spatial
localization plateau restrictions, one obtains a strictly positive `tau_*`
such that

\[
 |u_\theta(x_\tau,1-\tau)|
 \ge {e_0\over2}\tau^{-101/200},
 \qquad 0<\tau\le\tau_*.                            \tag{13}
\]

This uses a same-profile physical background estimate. It does not identify
the canonical and physical primary families and does not change the pressure
gauge.

## 5. Legal limit order and scope

For each fixed protected `tau`, set

\[
 t_\tau=1-\tau,\qquad T_\tau=1-{\tau\over2},
 \qquad 0<t_\tau<T_\tau<1.                          \tag{14}
\]

First fix `T_tau` and recover the pressure gradient. Next complete the
`R -> infinity` comparison limit at that fixed horizon, obtaining equality at
`(x_tau,t_tau)`. Only then let `tau` decrease to zero. The path tends to
`(0,1)` while (13) diverges. A globally smooth competitor is bounded on a
compact neighborhood of `(0,1)`, giving a contradiction.

Therefore, relative to every named continuum provider, no globally smooth
same-force zero-datum competitor in the fixed-interval finite-energy class
exists. This closes P7 as a conditional implication on Q181's physical
object. It does not construct a future field, independently reprove all named
providers, establish an unconditional admissible force, or prove
Navier--Stokes blow-up.

Next: perform a current-head provider-discharge audit of Q172--Q182, avoiding
Q171's stale frontier, and identify the first named continuum estimate still
only routed rather than independently proved.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q182_comparison_core_exclusion_replay/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q182_comparison_core_exclusion_replay/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q182_comparison_core_exclusion_replay/evidence.json \
  --output q182_replay.json
```
