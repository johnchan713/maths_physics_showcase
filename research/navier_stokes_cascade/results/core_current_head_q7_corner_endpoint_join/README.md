# Current-head Q7 Cartesian-corner and punctured-endpoint join

Parent: `ae11c819e913e349d72855727799563cc6574933`.
Status: `current-head-Q7-cartesian-corner-and-punctured-endpoint-joined-conditionally`.

Q6 binds the exact Q1--Q5 field to the full six-group local residual replay.
This checkpoint transfers that exact residual to the closed Cartesian
terminal corner and joins it to the punctured one-sided endpoint ledger of
the same locally finite field.  It retains the selected parameters
`h=1/200`, `kappa=1/100000`, complete pressure/covariance memory, actual
axial velocity, the enlarged compact-source requests, supports, gauges and
cutoff sequence.

The result remains conditional on the endpoint inequality in Section 3.
Exact tests freeze exponent arithmetic, derivative requests and provenance;
they are not a numerical or continuum proof of that inequality.

## 1. Selected Q6 residual on the closed chart

On the terminal chart

\[
 T+\eta^2=1,\qquad |\eta|\le1,\qquad
 L=1-2h\eta^2\ge1-2h={99\over100}.                    \tag{1}
\]

For each fixed residual jet `m` and target `N`, Q6 chooses one finite stage
`J`.  Its complete minimum exponent `mu_(m,N)` obeys

\[
 \sup_{K,|\eta|\le1}|D^m\mathcal R_1(u_*,p_*)|
 \le C_{m,N,K}q^{\mu_{m,N}},\qquad \mu_{m,N}\ge N+1. \tag{2}
\]

The same stage carries the two additional `A0` and `A1` source derivatives.
All six residual groups and the disjoint compact-source, regular-residual and
cutoff-tail accounts remain present.  The supremum uses a fixed compact `X`
interval and the positive normalized support margin.  It does not use a
positive physical support-distance floor or a positive pointwise edge-weight
floor, either of which would be scale-incompatible at the closed edge.

With

\[
 \rho=\tau+|z|^{1/D},\qquad D={1\over2}-h,
\]

the inherited implicit-coordinate estimate is `q<=4 rho`.  Therefore

\[
 |D^m\mathcal R_1(u_*,p_*)|
 \le 4^{N+1}C_{m,N,K}\rho^{N+1}.                      \tag{3}
\]

Since `N` is arbitrary, every fixed Q6 residual jet is Cartesian flat at the
terminal corner.  The stage is selected after `(m,N)` is fixed; no
stage-uniform constant or substitution `J=J(q)` is used.

## 2. Exact punctured endpoint requests

On each punctured terminal compact `|x|>=d>0`, retain the complete potential
`A`, the divergence-safe direct swirl `B`, and the full fixed-gauge pressure
`p`.  For a future localization derivative with `a` spatial and `j` time
derivatives, set `m=a+j`.  The exact finite Cauchy requests are

\[
                 A:m+4,\qquad B:m+3,\qquad p:m+2,     \tag{4}
\]

and the corner residual request is `m+1`.  The current endpoint ledger
supplies precisely (4).  Raising the request changes the finite comparison
stage, never the selected field, cutoff sequence, support or pressure gauge.

Complete moving pressure, both background-pressure tails, all six covariance
entries, fresh axial velocity and the exact exterior heat/pressure jets are
retained.  Endpoint field jets and corner residual flatness have different
roles: neither implies terminal regularity of velocity or pressure at the
origin.

## 3. Precise remaining premise

For every punctured terminal compact `E` and every fixed order, the endpoint
join still requires

\[
 \boxed{
 \sup_{E\times[1-\delta,1)}
 |D_x^\alpha\partial_t^j(A,B,p)|<\infty .}             \tag{5}
\]

This is the exact obstruction.  Interior smoothness for `t<1` and the flat
residual estimate (3) do not prove (5).  In particular, Q3 explicitly leaves
historical first-band coverage unspecified; the later-tail frame bounds do
not silently supply a uniform all-fixed-order endpoint atlas through
`|eta|=1`.

## 4. Correction to the older parameter-family wording

The earlier corner replay also sampled `h=9/1000`.  That is a useful
arithmetic diagnostic, but it is not the Q6 object.  The current-head join
uses only `h=1/200` and `kappa=1/100000`.  No conclusion for a family of
nearby parameters is promoted from the diagnostic sample.

## 5. Scope and next obligation

This closes the Q6 Cartesian-corner transfer and the same-field punctured
endpoint ledger conditionally on (5).  It does not replay localization,
construct a compact future force, prove the energy estimate, recover pressure
for comparison, exclude a competitor or prove blow-up.

The next justified obligation is the exact twelve-term localization replay
on this Q7 package.  It must retain unit viscosity, every potential cutoff
derivative, the full fixed-gauge pressure cutoff, both nonlinear cross
orientations and the curl-commutator square.  Smooth terminal force jets may
be claimed only after (5) is supplied for the actual endpoint providers.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q7_corner_endpoint_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q7_corner_endpoint_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q7_corner_endpoint_join/evidence.json \
  --output fresh.json
```
