# Q18 correction: harmonic sources do not inherit a full edge weight

Parent: `c5752be79a1370765eece93b1b5e91e8355b7767` (Q17).
The source-to-representative audit finds a mismatch in the claimed input
space, before Q14's derivative ceilings can justify Q16's component bounds.
This is a rigorous obstruction to an inherited implication, not a proof
that the actual construction fails. Historical records remain frozen.

## 1. The two distinct source contracts

`regular_source_foundation/README.md`, sections 1--3, defines the harmonic
class W_alpha with weight **sqrt(zeta) P**, finite edge poles and polynomial
S costs. Its complete regular linear remainder is in
W_(alpha+1/2-kappa). Its wave-wave products have weight zeta P^2 and may
therefore enter the full-weight mean class when their angular output is
zero. It explicitly retains pressure gradients, mixed viscosity, and the
complete curl remainder in the linear estimate.

In contrast, `core_log_renormalized_actual_source_induction/README.md`,
(2), and Q5 (10) assert a zeta^-1 normalized bound for **every** generated
source, including the complete primary linear remainder. The associated
`source_induction.py` checks epsilon rates but has no radial weight ledger.
A finite sum/product graph does not turn a half weight into a full weight.
The phrase "every primitive leaf has the full weight" cannot substitute
for the earlier, weaker W contract.

For a fixed band, write the radial part as

    zeta(x)=exp(-w^2/x^2-4/(L-x)^2), 0<x<L, w=C_setup^-20>0.

For theta>=0, let E_theta denote fixed jets bounded by
C zeta^theta d^-p, with finite p and d=min(1,x,L-x), allowing the existing
band powers, logs and Gaussian factors separately. At a fixed jet order:

| operation | justified radial exponent |
|---|---:|
| coefficient derivative | theta (with a larger finite pole) |
| finite sum | minimum of input exponents |
| product | sum of input exponents |
| multiplication by bounded background | theta |
| angular or auxiliary averaging | theta, absent a separate cancellation |
| signed target / original sqrt(zeta) root | theta - 1/2 |

The derivative row is a statement about the retained all-jet class, not
an inference from a value bound. No general inverse or square-root
operation is inferred from this table without its quantitative hypotheses.
For 0<zeta<=1, E_theta embeds in E_tau when theta>=tau.

Thus the source groups linear + wave-wave + mean-wave have exponents
1/2, 1, and 3/2 respectively if the mean has full weight. Their sum is
only guaranteed in E_(1/2). A stronger mean-wave weight cannot improve
the linear summand. Higher epsilon order and stronger radial weight are
independent properties.

## 2. Exact obstruction, even with arbitrary finite poles and logs

Take F(x)=sqrt(zeta(x)). It has smooth flat zero extension and all fixed
jets in E_(1/2), since each derivative of the exponential is itself times
a finite sum of inverse-distance powers. Nevertheless, for every finite
p>=0,

\[
 \zeta^{-1}d^p|F|=d^p\exp\{w^2/(2x^2)+2/(L-x)^2\}
 \longrightarrow\infty\quad(x\downarrow0).             \tag{1}
\]

For an explicit lower bound use x_n=w/n, for sufficiently large integers
n so d=x_n. Then

\[
 \log(\zeta(x_n)^{-1}d(x_n)^p|F(x_n)|)
 \ge n^2/2-p\log n+p\log w
 \ge n^2/2-pn+p\log w.                                \tag{2}
\]

For n>=4p this is at least n^2/4+p log w, hence diverges.
The width stays fixed at its original value. The w=1 rational checks are
fixtures for (2), not a replacement of that width in the construction.
For any fixed accepted band, epsilon^sigma, a finite logarithmic
renormalizer and a finite band constant are positive constants independent
of x. None can fix (1), however much epsilon reserve is spent. Increasing
the finite derivative/pole request likewise cannot supply the missing
exponential weight. Q14's additional 400009 derivatives do not repair it.

The same example times a nonzero angular harmonic and a supported smooth
slow/auxiliary factor has zero angular mean. Choose an interior point where
the fixed Gaussian and cutoff are nonzero; there the radial obstruction
is unchanged. Zero harmonic mean therefore does not strengthen the radial
weight of a nonzero harmonic source. Conversely a mean/covariance made
from two half-weight waves has the full weight by multiplication.

These are manufactured admissible class elements. We have NOT shown that
the actual linear remainder equals this example or saturates E_(1/2).
A genuine additional cancellation in the actual full linear remainder
could prove a stronger weight, but none follows from the cited mapping
estimate. This corrects the justification and scope of Q5's blanket
full-weight assertion, not its rational epsilon arithmetic.

## 3. Required typed repair and the precise remaining inequality

Keep nonzero harmonic source inputs in W_b with their sqrt(zeta) P
weight and whole-path f/P derivative contract. Keep compact mean,
covariance and signed-target inputs in their separately proved full-weight
classes. Wave-wave products can supply that weight to the mean system;
a single linear wave term cannot be silently substituted for such a product.
Retain all six covariance entries, complete pressure, gamma_d-A1 gamma_d,
the moving c_rho P product and both Pi1/Pi2 tails at every recomputation.

The signed root step makes this distinction essential. A full-weight
signed target divided by its original half-weight root has half weight
and flat extension. A merely half-weight target divided by that root
can be constant and need not extend by zero. Smoothness/compactness of the
numerator alone is insufficient. This agrees with the earlier sharp
threshold in `core_finite_support_closure`, section 2.

The corrected minimum source obligation is, for each fixed stage and
corrected finite request R,

\[
 \sup_{\ell\ge\ell_*}\epsilon_\ell^{-b_j}
 (1+|\log\epsilon_\ell|)^{-K_{j,R}}
 \max_{|I|\le R}\sup_{\Omega_{j,\ell}}
 \zeta^{-1/2}P^{-1}d^{p_I}|D^I F^{\rm harmonic}_{j,\ell}|<\infty.
                                                               \tag{3}
\]

Here P is the Gaussian envelope, not the scalar radial pressure defect.
Derivative bounds for the lifted f/P and full endpoint extension remain
additional inverse requirements, not consequences of this value notation.
For each *signed target* Xi the separately necessary retained full-weight
bound uses zeta^-1, without borrowing (3). If one instead retains Q5's
uniform full-weight assertion for all sources, one must explicitly prove
its stronger version of (3) for the actual **complete linear** remainder;
(1) proves that the W estimate alone cannot do so.

Q15's half-weight flat-edge conversion, Q16's conditional physical
exponents and Q17's projected counts remain usable. They do not prove
(3), its signed-target counterpart, or the typed propagation at every
node. The next checkpoint must propagate these distinct weights through
the actual sourced/signed/mean recursion on Q14's expanded domains.
No actual all-component bound, endpoint closure, unconditional residual
improvement, force theorem or blow-up follows here.

## 4. Validation

The exact checks exercise the typed algebra, a nonzero-harmonic failure
control, signed-root threshold and the divergent lower bound with arbitrary
finite poles. They are diagnostics for the analytical argument, not sampled
continuum suprema. The Q17 evidence and inherited hash graph are replayed.

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q18_edge_weight_obstruction/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q18_edge_weight_obstruction/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q18_edge_weight_obstruction/evidence.json --output FRESH.json
```
