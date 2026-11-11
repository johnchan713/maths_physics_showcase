# Q162: same-physical diagonal, endpoint and localized-force join

## Scope and correction

Q153--Q155 recorded the complete-potential diagonal, punctured endpoint and
localized-force conclusions conditionally because their physical target,
generated-source and fixed-stage tail premises had not then been closed on one
object. Q160/Q161 now close those mutually dependent premises by a simultaneous
finite-cycle induction, with the same labels, support holes, cutoff sequences,
fixed potential and pressure gauges. This checkpoint replays the downstream
chain on that exact object. It does not change the frozen Q153--Q155 evidence.

The conclusion is a continuum composition of the named analytical providers.
The tests reproduce exact rational identities and dependency hashes; they are
not a substitute for those estimates.

## 1. Complete-potential diagonal

At finite stage `j`, retain the fixed-gauge representative

```text
Delta u_j = curl A_j + b_j e_theta,       pi_j = fixed-gauge pressure.
```

For the complete diagonal, choose the potential cutoff scale only after all
finite constants required at that index have been computed, and set

```text
V_j = curl(chi(a_j q) A_j) + chi(a_j q) b_j e_theta,
P_j = chi(a_j q) pi_j.
```

This construction is divergence free term by term. The certified diagonal
provider proves that the scale `a_j` can absorb every finite stage constant for
all previously requested jets; a single diagonal choice handles the countable
list of fixed jet requests. Local finiteness, rather than convergence of a raw
uncut series, defines the resulting field.

The comparison stage is fixed before `q -> 0`. If `U` is that fixed stage,
`e` the locally finite tail and `pi` its pressure tail, the exact residual
difference is

```text
e_t - Delta e + grad pi
  + (U dot grad)e + (e dot grad)U + (e dot grad)e.
```

Thus the diagonal estimate retains time differentiation, unit viscosity,
complete pressure memory, both nonlinear crosses and the new square. It gives,
for every fixed Cartesian jet and every `N`, an `O(q^N)` corner bound. The
cutoff sequence is chosen once; it is not reselected after taking a limit.

## 2. Punctured closed-time endpoint

The implicit chart has lower Jacobian bound `1-2h=99/100` for the selected
`h=1/200`. On a punctured terminal compact only finitely many dyadic boxes,
bands and construction stages meet the point. The enlarged band convention is
kept exactly:

```text
1/4 <= q/Q <= 4.
```

For example, `c=2^-24`, `q_max=2^-20` and first band `ell_0=18` give the closed
window `18 <= ell <= 26`. With `q_*=1/8`, the observation choice satisfies both
`16 q_max < q_*` and `4 q_max < 2^-ell_0`.

For a requested mixed derivative of total order `m`, the endpoint ledger asks
for finite estimates through

```text
A: m+4,    direct swirl: m+3,    pressure: m+2,
corner residual: m+1.
```

The potential, divergence-safe direct swirl, velocity and fixed-gauge pressure
therefore have compatible one-sided jets on each punctured terminal compact.
The exterior pressure retains the full time-Leibniz convolution

```text
S_j = sum_{k=0}^j binom(j,k) A_k A_{j-k},
A_j = (h)_j (1+h)_j.
```

In particular `S_2 > 2 A_2`; replacing it by two copies of the velocity
coefficient would lose pressure memory. Also

```text
a_{j+1}/a_j = 4(h+j)(1+h+j)/((j+1)r^2),
```

which grows linearly in `j`. Hence the exterior time Taylor series has zero
radius at fixed nonzero `r`. Only smooth jet matching is licensed. This
argument gives no regular terminal velocity or pressure at the origin.

## 3. Exact localization and force extension

Let `c` be the axisymmetric space-time cutoff and put

```text
u = curl(c A) + c B,       p = c p0,
d = grad(c) cross A.
```

The exact unit-viscosity localized residual has twelve retained groups:

```text
c R_1(u0,p0), c_t u0, d_t,
-(Delta c)u0, -2 sum_i c_i partial_i u0, -Delta d,
p0 grad c, c(c-1)(u0 dot grad)u0,
c(u0 dot grad c)u0, c(u0 dot grad)d,
c(d dot grad)u0, (d dot grad)d.
```

The apparent thirteenth term vanishes because
`d dot grad(c)=0`. This ledger keeps the pressure-cutoff term, both nonlinear
crosses and the cutoff-curl square. The localized velocity is divergence free
and has zero initial datum.

The punctured terminal jets from Section 2 and the arbitrary-power Cartesian
corner flatness from Section 1 give smooth compact terminal jets for the
localized force. A single Borel schedule is chosen after the finite norm at
each index; every fixed mixed derivative then has a geometric tail. This
constructs one compact smooth global force matching the exact preterminal
residual. It does **not** construct a future velocity or pressure and does not
repair the terminal singular origin.

Writing

```text
F(t) = integral_0^t ||f(s)||_2 ds,
E(t) = ||u(t)||_2^2,       D(t) = ||grad u(t)||_2^2,
```

the retained energy estimate is

```text
E(t) + 2 integral_0^t D(s) ds <= F(t)^2.
```

For the diagnostic value `F=7/3`, the separate consequences are
`sup E <= 49/9` and `integral D <= 49/18`. Their sum is not asserted as a
single simultaneous budget.

## Result and remaining obligation

This closes frontier item P6 on the same physical field, potential/pressure
gauges, cutoff sequence and support data that close P3--P5. It establishes a
locally finite divergence-free preterminal field with all fixed complete
residual jets flat, compatible punctured endpoint jets, exact localization,
zero initial velocity, a smooth compact force, and the preterminal energy
estimate.

It does not establish a regular terminal field, a future solution, a globally
improved residual, competitor exclusion or Navier--Stokes blow-up. P7 remains:
replay fixed-gauge pressure recovery, expanding-ball comparison and the
protected-core lower bound on this exact localized object, in the order
`T<1`, then ball radius to infinity, then `tau` to zero.

## Reproduction

```bash
python3 research/navier_stokes_cascade/results/core_current_head_q162_diagonal_endpoint_force_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q162_diagonal_endpoint_force_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q162_diagonal_endpoint_force_join/evidence.json \
  --output q162_replay.json
```
