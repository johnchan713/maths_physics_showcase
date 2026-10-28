# Current log-reserved corner and endpoint join

Parent: `a08d07d296ae3da7b77a3d7578145f149ea4cf01`.
Status: `current-log-reserved-corner-and-punctured-endpoint-joined-conditionally`.

The parent instantiates the complete six-group residual estimate on the
current leading profile, full-path maps, same-choice projectors and weighted
compact/slow inverses. This checkpoint transfers that exact estimate to the
Cartesian terminal corner and joins it to the already frozen punctured
one-sided endpoint jets of the same current field.

The result remains conditional on the frozen core/root, regional/frame and
all-fixed-order physical bounds. Exact tests freeze exponent arithmetic,
derivative requests and provenance. They do not prove those continuum
inputs numerically.

## 1. Transfer the corrected current stage to the corner

On the closed terminal chart,

\[
 T+\eta^2=1,\qquad |\eta|\le1,\qquad
 L=1-2h\eta^2\ge1-2h>49/50.                         \tag{1}
\]

For each fixed residual jet `m`, target power `N`, and finite logarithmic
degrees, the parent chooses one finite comparison stage `J`. The current
source request at this same `J` includes the two additional derivatives
required by the extra `A0` and `A1` integrations. Its compact-source,
regular-residual and cutoff-tail reserves remain disjoint. Hence

\[
 \sup_{K,|\eta|\le1}|D^m\mathcal R_1(u_*,p_*)|
 \le C_{m,N,K}q^{N+1}.                                \tag{2}
\]

The supremum uses the normalized support margin and upper weighted edge
seminorms. It does not use a positive physical support-distance floor or a
positive pointwise edge-weight floor. Complete pressure memory, all six
covariance entries, actual axial velocity, time, viscosity, both cross
orientations and the correction square remain in (2).

With

\[
 \rho=\tau+|z|^{1/D},\qquad D=1/2-h,
\]

the unchanged implicit-coordinate estimate gives `q<=4 rho`. Therefore

\[
 |D^m\mathcal R_1(u_*,p_*)|
 \le 4^{N+1}C_{m,N,K}\rho^{N+1}.                      \tag{3}
\]

Since `N` is arbitrary, every fixed current residual jet is Cartesian flat
at the terminal corner. The stage is selected after `(m,N)` is fixed; there
is no substitution `J=J(q)` and no stage-uniform constant.

## 2. Join the current punctured endpoint ledger

For every punctured terminal compact `|x|>=d>0`, the current endpoint
checkpoint conditionally supplies uniform one-sided jets of the complete
potential `A`, the separate divergence-safe direct swirl `B`, and the full
fixed-gauge pressure `p`. The exterior retains the exact heat swirl and the
complete pressure convolution. No terminal regularity of velocity or
pressure is asserted at `x=0`.

The precise inherited analytical input is

\[
 \sup_{E\times[1-\delta,1)}
 |D_x^\alpha\partial_t^j(A,B,p)|<\infty               \tag{4}
\]

for every fixed order on each punctured compact. Equation (4) is conditional
on uniform all-fixed-order physical bounds up to `|eta|=1` on every positive
`q` compact. Interior smoothness and corner residual flatness alone do not
prove it; if the inherited physical estimates are not uniform there, (4) is
the exact blocking inequality.

For a future localized-force derivative with `a` spatial and `j` time
derivatives, put `m=a+j`. The exact finite Cauchy requests are

\[
 A:m+4,\qquad B:m+3,\qquad p:m+2,                     \tag{5}
\]

and the corner residual request is `m+1`. The current endpoint ledger gives
exactly (5). Raising the request may raise the finite comparison stage, but
does not select a new field, cutoff sequence, support, or pressure gauge.

Thus the compatible closed-time input package consists of current punctured
`A,B,p` jets, exact exterior heat/pressure jets, and complete residual
flatness at the origin. These roles remain distinct: the endpoint field jets
do not regularize the origin, and the corner residual estimate does not
supply away-from-origin endpoint limits.

## 3. Scope and next obligation

This checkpoint closes the corrected corner and punctured-endpoint transfer
on the exact current objects, conditionally on (4) and the frozen continuum
contracts. It does not replay the twelve-term localization identity,
construct or extend the compact force, prove the energy estimate, recover
pressure for comparison, exclude competitors, or prove blow-up.

The next obligation is to instantiate the twelve localization terms on this
joined package, retaining unit viscosity, every cutoff derivative, the full
pressure cutoff term, both nonlinear crosses and the curl square. Only after
smooth terminal force jets are revalidated should future force extension,
energy and expanding-ball comparison be advanced.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_log_reserved_corner_endpoint_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_log_reserved_corner_endpoint_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_log_reserved_corner_endpoint_join/evidence.json \
  --output core_current_log_reserved_corner_endpoint_join_fresh.json
```
