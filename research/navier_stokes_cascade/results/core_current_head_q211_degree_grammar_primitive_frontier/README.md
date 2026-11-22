# Q211: polynomial degree grammar and the weighted primitive frontier

Parent: Q210, `c4c86e6926d1828605aa9eca600c4e9d85f2d6c2`.
The live branch, Q210 and its active CI, Q192--Q202, the regular-source and
edge-envelope records and their frozen evidence were read first. Historical
evidence, the selected profile, roots, pressure gauge, supports, cutoff
sequences and unchanged `C^-20` activation width remain fixed.

Q210 correctly proves that per-band smoothness does not imply a normalized
uniform all-band estimate. Its proposed next step, however, starts again with
the initial primary-plus-mean source and one sourced/generic cycle. Those
conditional steps already exist at Q196--Q202. This checkpoint avoids
repeating them. It proves the finite polynomial-degree algebra behind those
steps and lowers the independent frontier to the actual weighted primitive
degree array which that algebra consumes.

## 1. Finite polynomial degrees form a closed algebra

For a fixed finite request write `a_n` for a polynomial `S` degree controlling
all jets of one factor through total order `n`. Leibniz' rule gives the exact
product majorant

\[
 (a\star b)_n=\max_{0\le i\le n}(a_i+b_{n-i}).       \tag{1}
\]

If `g` has a positive value floor, degree-zero value bound, and derivative
degrees `g_n`, differentiating `g g^{-1}=1` gives

\[
 r_0=0,\qquad
 r_n=\max_{1\le i\le n}(g_i+r_{n-i}).                \tag{2}
\]

The binomial coefficients change constants but not degrees. The same finite
maximum principle applies to Bell quotients, finite compositions and the
two-edge weighted primitives once their exact moment condition is retained.
Thus no finite sum, product, reciprocal with a positive floor, quotient or
weighted primitive creates superpolynomial `S` growth from polynomial
inputs.

For a zero-data sourced Duhamel solve, let `c_i` be coefficient degrees and
`a_i` source degrees. The differentiated forcing and solution degrees obey

\[
 f_n=\max_{0\le i\le n}(c_i+a_{n-i}),                \tag{3}
\]

\[
 t_0=f_0+1,\qquad
 t_n=\max\left\{f_n+1,
    \max_{1\le i\le n}(c_i+1+t_{n-i})\right\}.      \tag{4}
\]

Q192 supplies the conservative coefficient ledger

\[
 c_0=0,\qquad c_i=2i+4\quad(i\ge1),                 \tag{5}
\]

and Q193 supplies homogeneous amplitude degree `7i`. For any finite source
array through order `R`, choose

\[
 D_R=\max_{0\le i\le R}(a_i-7i,0),
 \qquad a_i\le7i+D_R.                               \tag{6}
\]

Because `2i+5<=7i` for every `i>=1`, (3)--(4) give the explicit sufficient
bound

\[
 f_i\le7i+D_R,qquad t_i\le7i+D_R+1.                \tag{7}
\]

This is an all-fixed-order degree estimate after the finite primitive source
array is supplied. It does not evaluate the optimized offset `D_R` or the
associated continuum constants.

## 2. The initial and first sourced cycles are already present

The live route contains the requested conditional work:

* Q196 derives the complete primary-plus-initial-mean quotient bound at
  log-renormalized order `1-kappa`, including smooth endpoint extension, and
  obtains a uniform pure-power bound after one `kappa` reserve.
* Q197 solves the first constrained sourced wave with the `n dot f` pressure
  term, recomputes the complete six-entry covariance, compact pressure,
  axial moment and current defects.
* Q198--Q200 correct the generic signed target and bind it to the physical
  same-family weighted state-difference route.
* Q201 recomputes post-signed fast and slow means with axial velocity,
  `Pi1/Pi2`, complete pressure and the whole moving product.
* Q202 joins the next constrained sourced cycle on the corrected physical
  domain.

Equations (1)--(7) explain why those finite grammars retain some finite
polynomial `S` degree. Repeating the initial source or first sourced cycle
would not independently validate their primitive inputs. Their conclusions
remain conditional; tests and source hashes do not prove their continuum
suprema.

## 3. Precise lower obstruction

The actual independent input is the weighted primitive degree array. For the
finite request `R`, define schematically

\[
 D_R^{prim}=\max_{|I|\le R}\inf\left\{d\ge0:
 \sup_{\ell\ge\ell_*}S^{-d}{\cal N}_I(\ell)<\infty
 \right\}.                                           \tag{8}
\]

Here `N_I` is the actual edge-weighted norm of the current target, original
root, frame and mean data on the physical band. The required current-target
identity is

\[
 X_e^{new}-X_e^{old}=-I_eQ_e\overline{\Delta E_e}.  \tag{9}
\]

A direct proof of (8)--(9) must retain all six covariances, radial and axial
velocity, complete pressure, both `Pi1/Pi2` tails, the derivative of the whole
moving pressure product, exact zero weighted totals, both flat-edge oriented
primitives and the reserved interior support.

Q199 proves that unweighted fixed-jet bounds, exact moments and compact
support do not imply (8). Q200 supplies a different conditional route: take
the same-family state difference before estimating, so each covariance
difference retains the full `zeta` weight and its zero weighted total. The
degree algebra above can propagate that input; it cannot manufacture it.

The first nonrepeated obligation is therefore to derive (8) for the first
physical current target directly from the explicit state-difference,
pressure and moment maps, and only then use (1)--(7) for later cycles. Q211
does not certify the actual primitive degree array, a uniform actual-pulse
hierarchy, the complete residual, an admissible force, competitor exclusion
or Navier--Stokes blow-up.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q211_degree_grammar_primitive_frontier/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q211_degree_grammar_primitive_frontier/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q211_degree_grammar_primitive_frontier/evidence.json \
  --output fresh.json
```
