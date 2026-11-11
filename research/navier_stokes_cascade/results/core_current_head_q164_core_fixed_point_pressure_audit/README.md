# Q164 current-head core contraction and leading-pressure audit

Parent: `6d0dacc5b2ec662857cb40eafdc750c22bcc3c84`.
Status: `Q164-all-index-core-contraction-and-leading-pressure-audited`.

Q163 named an independent all-index audit of the core fixed-point map,
original-equation recovery and same-gauge leading pressure as its next gate.
That frontier was too pessimistic relative to the live repository.  The
independent all-index source ledger and contraction had already been
rederived in `core_fixed_point_ledger_recheck`; `core_contraction_audit`
recovers the original equations and regular forward stresses; and the
pressure identification and same-field reconciliation prove the leading
same-gauge pressure identity.  Q164 checks their exact inequalities and binds
those records to the current Q1--Q163 object.

This is an analytical checkpoint for the selected exact pressure family.  It
does not alter frozen evidence or claim that finite tests prove the continuum
estimates.  It also does not identify the complete corrected pressure with
the leading centrifugal pressure: all positive-order background, wave, mean
and localization pressure terms remain.

## 1. All-index coefficient algebra

The independently rederived coefficient norm uses

\[
 a_{nk}=\frac{20^{-n}\rho^{-k}k!\binom{n+k}{k}}
 {(n+1)^2(k+1)^2}.
\]

For a Leibniz split, the binomial factors satisfy

\[
 \binom{i+j}{j}\binom{n-i+k-j}{k-j}\leq\binom{n+k}{k}.
\]

Each squared-denominator convolution is below 16, so applying it in the
radial and angular derivative indices gives the all-index product estimate

\[
 \|FG\|_\rho\leq A_0\|F\|_\rho\|G\|_\rho,
 \qquad A_0=256.
\]

The exact radial shifts give `D0=80/rho`.  For the worst integrated mixed
word the remaining allocation is

\[
 \frac{(i+1)j}{(i+j+1)(i+j+\nu)}\leq1,
 \qquad \nu\in\{1,2\}.
\]

Thus every integrated derivative word is bounded without asserting that a
bare parameter derivative is bounded.  For
`T=J2*(multiplication by chi)/2`, repeated application raises the radial
degree and gives

\[
 \|T^k\|\leq \frac{30720^k}{k!(k+1)!},\qquad
 \|(1+T)^{-1}\|\leq e^{30720}.
\]

The finite index samples in the regression code check arithmetic only; the
displayed combinatorial estimates carry the all-index quantifiers.

## 2. Complete source and pressure ledger

Direct substitution into the normalized core equations yields exactly 21
source terms.  Their exact coefficient sums are

| Group | Exact sum |
|---|---:|
| `R1` fixed linear | `19237/50` |
| `R1` ordinary quadratic | `952/25` |
| `R1` zeta quadratic | `10 M` |
| `R2` fixed linear | `442` |
| `R2` ordinary quadratic | `32` |
| `R2` pressure | `16+10+8=34` |

The pressure primitive contains four factors `g*g*Phi*Phi`, hence three
algebra products and one radial integration:

\[
 \|p\|\leq80A_0^3S^2.
\]

The three pressure source terms are `-4 A eta p`, `d p_eta`, and
`-2 eta Dp`.  Multiplying by their fixed coefficients costs one more
algebra product, so the source coefficient is

\[
 34\cdot80A_0^4=2720A_0^4.
\]

For two states, replacing either `Phi` factor contributes `2S`; the pressure
Lipschitz coefficient is therefore `5440 A0^4 S`.  Reducing either algebra
power, dropping the factor two, or removing any pressure row gives a false
ledger.

## 3. Invariant ball and strict contraction

Let

\[
 S=2+e^{30720}+40Z_b,
\]

and let `Fb` and `Lb` be the full displacement and Lipschitz numerators made
from all source terms above.  The exact map is

\[
 (\Phi,u)\mapsto
 \left(\Phi_0+\frac{(1+T)^{-1}J_2R_1}{2\Lambda},
 u_0+\frac{J_1R_2}{2\Lambda}\right).
\]

Discarding the factors `1/2` is conservative.  With
`epsilon=10^-16`, `H=K^(1/16)>1`, and

\[
 \Lambda=10^{12}H^{128}\epsilon^{-2},
\]

the selected parameter schedule gives both displacement and Lipschitz bounds

\[
 \frac{E\max(S,F_b,L_b)}{\Lambda}
 \leq\frac{\epsilon^2}{10^{18}H^{127}}
 =\frac{10^{-50}}{H^{127}}<10^{-50}.
\]

The closed unit ball is invariant and the map is a strict contraction in the
complete infinite coefficient space.  Banach's theorem supplies the unique
fixed point.  A finite triangular Taylor iteration is only a diagnostic and
is not used for convergence.

## 4. Original equations and forward stress

Absolute convergence on a common smaller analytic neighborhood permits the
radial differential operators to be applied to the fixed point.  Undoing the
normalization gives exactly

\[
 -2L\frac{X\phi_{XX}+2\phi_X}{\phi}=S_q,
 \qquad -2L(XU_{XX}+U_X)=S_n,
 \qquad \Pi_X=\frac{\phi^2}{C^2}.
\]

The regular forward candidates

\[
 Q=-2L\phi_X/\phi,\qquad N=-2LU_X
\]

satisfy the defining source equations.  Their homogeneous differences from
the forward primitives are a parameter constant divided by `XH` and by `X`.
Axis regularity forces both constants to vanish.  Hence the leading core
stress is exactly `T0=0`; no small activation factor is used to hide an
approximate core defect.

## 5. Same-gauge leading pressure

For the difference `delta` between the core pressure datum and the completed
leading backward pressure, the five exact leading rows reduce the remaining
angular equation to

\[
 4A\eta\delta-(1-\eta^2)\delta'=0,
 \qquad A=101/200>0.
\]

Every nonzero homogeneous solution

\[
 \delta=C(1-\eta^2)^{-2A}
\]

is unbounded at the endpoints.  Smoothness forces `C=0`, proving
`Pi_axis+C_p,infinity=0`.  Together with
`Pi_X=E_0^2/(2X)` on the protected core, this identifies the local core
pressure with the completed global leading backward pressure on the same
field.

This conclusion is deliberately limited to the leading pressure.  The full
pressure still retains all later background, wave, mean and localization
terms, including those entering the localized forced pressure equation.

## 6. Result and remaining frontier

Q164 closes Q163's U1 obligation by binding the existing independent
all-index rederivation, equation recovery and leading-pressure proof to one
current-head object.  No omitted core source or false pressure constant was
found.

The next obligation is the independent finite-map audit.  It begins with
Q2's positive core, determinant, orientation, same-tilt cone and factored
flat-edge continuum infima over the actual finite regional cover.  Q3 then
requires the original full-frame primitive estimates, complete path and phase
jets, and accepted physical domain.  Until those and the later chain are
independently closed, no unconditional pulse bound, complete residual
improvement, admissible-force theorem, competitor exclusion or blow-up
result is claimed.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q164_core_fixed_point_pressure_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q164_core_fixed_point_pressure_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q164_core_fixed_point_pressure_audit/evidence.json \
  --output /tmp/q164_fresh.json
```
