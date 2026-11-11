# Q167 current-head Cartesian-corner and punctured-endpoint audit

Parent: `8d114de77700e6cef29d40503b4f94fd85c8db19`.
Status: `Q167-current-head-cartesian-corner-and-punctured-endpoint-audited-conditionally`.

Q166 closed the complete six-group local residual frontier conditionally and
named Q7's Cartesian-corner and endpoint premise next.  Q7 correctly left a
uniform punctured closed-time inequality open at its own head.  The live
branch later supplied that inequality in Q154 and bound Q153--Q155 to the
same physical field in Q162.  Q167 audits the full identity chain and joins
those results directly to Q166.

This is an analytical composition relative to named continuum providers.
Exact arithmetic, frozen hashes and passing tests do not prove the underlying
continuum seminorms.

## 1. Same physical object

The evidence path is

```text
Q166 -> Q165 -> Q164 -> Q163 -> Q162 -> Q154 -> Q7,
```

and Q166 and Q7 use the byte-identical Q6 evidence.  Q162 explicitly freezes
the same physical field, complete potential, divergence-safe direct swirl,
fixed pressure gauge, supports and diagonal cutoff sequence.  Therefore
Q154's endpoint theorem applies to the Q164--Q166 object; it is not a theorem
for a parallel reconstruction.

## 2. Closed Cartesian corner

On the selected chart

\[
 h={1\over200},\quad \kappa={1\over100000},\quad
 L=1-2h\eta^2\ge {99\over100}.
\]

For each fixed residual jet `m` and power `N`, Q166 chooses one finite stage
with complete exponent `mu>=N+1`.  Q7's exact implicit-coordinate estimate

\[
 q\le4\rho,\qquad \rho=\tau+|z|^{1/D},
\]

gives

\[
 |D^m\mathcal R_1|\le4^{N+1}C_{m,N,K}\rho^{N+1}.
\]

As `N` is arbitrary, every fixed residual jet is Cartesian flat at the
terminal corner.  The proof retains time, unit viscosity, full pressure,
both nonlinear crosses, the new square, enlarged `A0/A1` requests and all
physical memory.  It uses neither a positive physical support-distance floor
nor a positive pointwise edge-weight floor.

## 3. Q7's endpoint inequality is now supplied

The exact Q7 premise is, for every punctured terminal compact `E` and every
fixed `alpha,j`,

\[
 \sup_{E\times[1-\delta,1)}
 |D_x^\alpha\partial_t^j(\mathcal A,\mathcal B,p)|<\infty. \tag{1}
\]

Q154 proves (1) on the same field.  The full accepted atlas

\[
 {1\over4}\le {q\over Q}\le4
\]

gives the closed window `ell=18,...,26` for the representative diagnostic
compact.  Its observation cylinder satisfies both `16 qmax<qstar` and
`4 qmax<2^-ell0`, so no unspecified historical first band is imported.
On a punctured compact only finitely many bands, slow boxes and diagonal
stages are active after the derivative order is fixed.  Finite products,
complete curls, constrained path variation, pressure/moment maps and compact
inverses then preserve the one-sided bounds supplied by the physical
providers.

The exact future localization request remains

\[
 \mathcal A:m+4,\qquad \mathcal B:m+3,\qquad p:m+2,
 \qquad \mathcal R:m+1.                                \tag{2}
\]

## 4. Complete exterior pressure and the sharp limitation

In the heat exterior, with `A_j=(h)_j(1+h)_j`, the fixed-gauge pressure jet
contains the full convolution

\[
 S_j=\sum_{k=0}^j {j\choose k}A_kA_{j-k}.
\]

In particular `S_2-2A_2=2A_1^2>0`; replacing the convolution by two velocity
copies would drop pressure memory.  Every fixed one-sided jet has a compatible
limit, but the successive swirl Taylor-coefficient ratio is

\[
 {4\over r^2}{(h+j)(1+h+j)\over j+1},
\]

which grows linearly.  For nonzero exterior amplitude, the velocity and
pressure time Taylor series have zero radius.  Smooth jet gluing is required;
analytic continuation or zero extension is invalid.

## 5. Remaining frontier

Q167 closes the conditional Q7 corner and punctured-endpoint frontier.  It
does not regularize the terminal origin, replay localization, construct a
future velocity or pressure, prove compact forcing, energy, comparison,
exclusion or blow-up.

The next obligation is to audit Q8's exact twelve-term axisymmetric
localization, zero initialization and compact-force extension on this same
object.  That replay must retain unit viscosity, every potential-cutoff
derivative, full pressure, both nonlinear crosses and the cutoff-curl square.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q167_corner_endpoint_frontier_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q167_corner_endpoint_frontier_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q167_corner_endpoint_frontier_audit/evidence.json \
  --output q167_fresh.json
```
