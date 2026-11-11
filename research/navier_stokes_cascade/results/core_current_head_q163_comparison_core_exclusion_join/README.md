# Q163: same-object comparison and protected-core exclusion

## Scope

Q162 closes the physical dependency frontier P1--P6 relative to the named
continuum providers and leaves P7: fixed-gauge pressure recovery, expanding-
ball comparison and the protected-core limit on that exact localized field.
Q156 and Q157 already contain the required comparison and core estimates, but
they were frozen before Q160--Q162 supplied the physical target, source, tail,
diagonal, endpoint and force premises together. This checkpoint replays P7 on
the Q162 object. Historical evidence is unchanged.

The conclusion remains conditional on the named continuum construction
providers. Exact arithmetic, tests and dependency hashes verify interfaces and
provenance; they do not prove those provider theorems.

## 1. Complete fixed-gauge pressure comparison

Let `u,p,f` be the Q162 localized field, fixed-gauge pressure and compact
smooth force on `t<1`, and let `v,P` be a locally smooth competitor with the
same force and zero datum and with a finite `L2` bound on every fixed interval
`[0,T]`, `T<1`. Put

```text
w = v-u,
pi = P-p,
g = w tensor w + w tensor u + u tensor w.
```

The common force divergence cancels exactly once, giving

```text
w_t + div g = Delta w - grad pi,
-Delta pi = partial_i partial_j g_ij.
```

The complete constructed pressure is not replaced by an unforced formula.
The pressure difference satisfies

```text
grad pi = grad sum_ij R_i R_j g_ij
```

in space-time distributions. This recovery uses the momentum equation to
remove the harmonic-gradient ambiguity; no decay assumption on the competitor
pressure is made. A time-dependent spatial pressure constant has zero cutoff
flux.

Split `g=w tensor w+h`, where `h=w tensor u+u tensor w`. Since the localized
`u` has fixed compact spatial support, so does `h`. Its pressure flux on the
outer annulus is `O(R^-5/2)`. The quadratic pressure commutator, transport and
cutoff terms give

```text
0.5 E_R' + 0.5 D_R^2
  <= ||grad u||_infinity E_R + C_T R^-2.             (1)
```

The exact ledger uses the eighth cutoff power. The two nonconstant Young
pairs are `(a,b)=(3/2,1)` and `(3/4,7/4)`, producing decays `R^-4` and
`R^-14/5`; both Sobolev offsets and the compact mixed pressure decay as
`R^-5/2`. The coarse cutoff-Laplacian term is the limiting `R^-2` remainder.
No competitor derivative, dissipation or boundary-flux bound enters (1).

For fixed `T<1`, Gronwall yields

```text
E_R(t) <= 2 C_T T exp(2 T L_T) R^-2,
L_T = sup_[0,T] ||grad u||_infinity < infinity.       (2)
```

Sending `R` to infinity at this fixed `T` gives `v=u` on `[0,T]`.

## 2. Protected core on the same localized field

The Q157 path uses

```text
A = 101/200,      X_in = 1/101 < X_a = 4/101,
x_tau = (sqrt(2 X_in tau),0,0),    t=1-tau.
```

At `z=0`, the implicit chart has `q=tau`, `eta=0`, and fixed `X=X_in`.
Every annular diagonal correction vanishes in a neighborhood of this path
because its closed `X` support begins at `X_a`. The positive-order Borel
background is small rather than zero:

```text
tau^A u_theta(x_tau,1-tau) = e0 + R_bg(tau),
|R_bg(tau)| <= C_bg tau^(2h),       e0>0.              (3)
```

The localization satisfies `c=1` and `grad c=0` on the late core. With

```text
tau_* = min(q0,(e0/(2 C_bg))^100,tau0/2,r_L^2/(8 X_in)) > 0,
```

(3) gives

```text
|u_theta(x_tau,1-tau)| >= (e0/2) tau^(-101/200)
for 0<tau<=tau_*.                                      (4)
```

This binds the path to the same Q162 cutoff sequence, support hole,
localization and pressure gauge; it does not identify distinct canonical and
physical primary families.

## 3. Legal order of limits

The terminal contradiction does not require `C_T` or `L_T` in (2) to remain
bounded as `T` approaches one. For each fixed protected `tau`, set

```text
t_tau = 1-tau,      T_tau = 1-tau/2.
```

Then `0<t_tau<T_tau<1`. First recover pressure and send `R` to infinity in
(2) at fixed `T_tau`. This gives the exact identity
`v(x_tau,t_tau)=u(x_tau,t_tau)`. Only after that radius limit is complete is
`tau` sent to zero. A globally smooth competitor is bounded on one compact
space-time neighborhood of `(0,1)`, while (4) diverges because `A>0`.
Therefore, relative to the retained continuum providers, no globally smooth
same-force zero-datum finite-energy competitor exists.

This closes P7 as a conditional consequence. It does not construct a future
velocity or pressure and is not an unconditional forced-breakdown or
Navier--Stokes blow-up theorem.

## 4. Remaining foundation gate

The P1--P7 dependency frontier is now closed only **relative to** its named
analytical providers. The frozen comparison foundation still explicitly marks
`all_retained_analytic_dependencies_discharged=false`; Q1 binds the exact core,
roots and completed leading pressure on one object, but also states that it
does not independently reprove every provider theorem.

The next audit must therefore start below the implication chain rather than
repeat the exclusion. On the selected Q1--Q162 object, independently establish
for the all-index core map `T` and its coefficient norm `X`

```text
T(B_r) subset B_r,
||T(Phi)-T(Psi)||_X <= lambda ||Phi-Psi||_X,   lambda<1,       (5)
```

recover the original core equations from its fixed point, and prove in the
same gauge

```text
Pi_axis + C_p,infinity = 0.                               (6)
```

The existing Q1 join cites providers for (5)--(6) and checks their interfaces;
a fresh end-to-end derivation must determine whether those providers really
discharge this gate. Passing regression tests or matching dependency hashes
cannot substitute for that derivation. Afterward the finite defining,
existence, inverse and gap maps must receive the same treatment before any
unconditional theorem assessment.

## Reproduction

```bash
python3 research/navier_stokes_cascade/results/core_current_head_q163_comparison_core_exclusion_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q163_comparison_core_exclusion_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q163_comparison_core_exclusion_join/evidence.json \
  --output q163_replay.json
```
