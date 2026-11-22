# Q204: global ledger, uniform primary and covariance reconciliation

Parent: Q203, `fd7d3b924b1a731bcc02a68507fee685642b214a`.
The live branch, Q203, its active CI, the global pressure/moment recheck,
coefficient-support recheck, Q148--Q153 and the Q161/Q180/Q202 physical
cycle were read first. Historical evidence, the finalized profile, physical
family, pressure gauge, roots, supports, accepted bands and unchanged
`C^-20` width are preserved.

Q203 names the global pressure/moment, coefficient/reserved-support and
uniform actual-pulse ledgers next. Those obligations were developed earlier
in the branch, but their exact relation to Q203's corrected physical domain
was not recorded at the current head. This checkpoint verifies the frozen
same-object chain, joins the completed ledgers, and corrects the remaining
frontier: the implemented rule-defined Borel schedule already discharges the
apparent coefficient-scale inequality and gives uniform primary-pulse bounds
on every used band. It does not give uniform constants for every generated
cycle or absolute pulse decay.

All analytical conclusions remain conditional on the named continuum
providers. Exact identities, rational tests and evidence hashes certify the
composition and bookkeeping; they are not a continuum proof of the provider
estimates.

## 1. Frozen same-object chain

The dependency hashes give the exact chain

```text
Q203 -> Q202 -> Q180 -> Q161 -> Q153 -> Q152 -> Q149 -> Q148
     -> coefficient-support recheck -> global pressure/moment recheck.
```

Every arrow uses the current frozen child digest. Thus the global rows,
coefficient supports, physical primary pulses, complete covariance/mean
updates and Q203 diagonal refer to the same finalized background and physical
family. No canonical-family transfer or new pressure representative is used.

## 2. Global pressure, moments and exterior stresses

The finite matching rows are ordered as `(M,I,J,S,Cp)`, while the global heat
row is the renormalized moment

\[
 R=\int_0^\infty(H-H_{pow})\,dX
  =\lim_{X\to\infty}\left[I(X)-{XH_{pow}(X)\over1-h}\right]. \tag{1}
\]

The raw `I(infinity)` diverges and is never set to zero. For two profiles
with the identical future heat tail and identical power subtraction,

\[
 \Delta R=\Delta I(X_{rep})+
 \int_{X_{rep}}^\infty\Delta(H-H_{pow})\,dX.          \tag{2}
\]

Therefore the exact finite endpoint match gives `Delta R=0` only because
both common-tail facts are retained. The same repaired physical profile has

\[
 M_\infty=R=J_\infty=S_\infty=0,
 \qquad \Pi_{axis}+C_{p,\infty}=0                    \tag{3}
\]

as functional identities in `eta`. The original gauge then gives

\[
 \Pi(X)=-\int_X^\infty {E(x)^2\over2x}\,dx.          \tag{4}
\]

Intermediate heat, shape and compact edits may have nonzero pressure memory;
none is reset at an interface. All five defects and their angular derivatives
are needed to remove radial velocity and both exterior stress memories. In
particular, an axial pressure defect with zero value but nonzero `eta`
derivative produces a term linear in `X`.

## 3. Finite coefficients and reserved means

At each fixed positive coefficient order, the ordinary five-row correction
is the linear solve

\[
 c_n=-B^{-1}z_n,
 \qquad z_n=(M,(M-J)/\lambda,I,-S,C_p),              \tag{5}
\]

with the continuum determinant bound `||B^-1||<800`. It needs no small-target
condition. The exact right supports are

| Order | Velocity | Pressure | Stress |
|---|---|---|---|
| `1` | `Xvel` | `Xv` | `Xb` |
| `2` | `Xvel` | `Xv` | `Xv` |
| `n>=3` | `Xvel` | `Xvel` | `Xvel` |

Positive-order velocity is zero on the reserved I4 interval, including
inside potential-cutoff transitions, because `M_n` and all its fixed angular
jets vanish there. Pressure is reconstructed rather than treated as a
velocity component. On I4 it remains

\[
 \Pi_1=K_1(\eta)-{\mu^2\over4X},\qquad
 \Pi_2=K_2(\eta),\qquad \Pi_n=0\quad(n\ge3).          \tag{6}
\]

Both `K1` and `K2` carry later-source memory. For any `Xa<Xb` in I4, their
unknown constants cancel from the exact two-point difference, yielding

\[
 \sup_{X\in I4}|q^{2A}p_B-\Pi_0|
 \ge { |\chi(c_1q)|q^{2h}\mu^2\over8\rho_3}
       (e^{-6}-e^{-11}).                              \tag{7}
\]

The factor `1/8` is sharp even after arbitrary `X`-independent cancellation.
Equation (7) is a rigorous obstruction to deleting the reserved pressure or
replacing the full mean pressure by `Pi0`; it is not a lower bound for the
full Navier--Stokes residual.

## 4. Correction: the Borel rule closes the primary schedule

Q148 displayed the sufficient row

\[
 C_{n,m}a_n^{-nh}\le2^{-n},\qquad m\le n.            \tag{8}
\]

Taken against a pre-existing numerical scale list, fixed-order finiteness
would not imply (8). The live construction does not use such a list. After
all finite seminorms at order `n` are known, it selects `b_n=log a_n`
strictly above

\[
 1+\max\left\{b_{n-1}+\log2,\ B+\log2,
 {\log C_{n,m}+n\log2\over nh}:0\le m\le n,
 \text{all finite physical rows}\right\}.            \tag{9}
\]

Hence on `q<=a_n^-1`,

\[
 C_{n,m}q^{nh}\le2^{-n}e^{-nh}<2^{-n}.               \tag{10}
\]

The locally finite background keeps every pressure coefficient and
cutoff-curl commutator. No uniform coefficient bound in `n` and no
convergence of the uncut formal series is used.

On the same finalized background, Q149's physical primary theorem supplies,
for every used band `Q=2^-ell` and fixed physical orders `a,b`,

\[
 |\nabla^a\partial_t^b u_p|
 \le C_{a,b}Q^{-(a+1)s_x-bs_t}S^B\sqrt\zeta\,\delta^{-B}P, \tag{11}
\]
\[
 |\nabla^a\partial_t^b p_p|
 \le C_{a,b}Q^{-(a+2)s_x-bs_t}S^B\sqrt\zeta\,\delta^{-B}P, \tag{12}
\]
\[
 |\nabla^a\partial_t^b\mathcal A_p|
 \le C_{a,b}Q^{-as_x-bs_t}S^B\sqrt\zeta\,\delta^{-B}P,   \tag{13}
\]

where `s_x=1/2+h/2`, `s_t=1+3h/2`. The constants may depend on the fixed
derivative and profile but not on the running band or pulse label. The
potential and pressure have fixed closed transverse, pulse, radial and slow
supports; curl does not enlarge them. This is a uniform primary-pulse
size/derivative/support theorem, not absolute smallness as `Q` tends to zero.

## 5. Complete covariance and current sourced cycle

The inherited physical family has its own exact same-family stress assembly

\[
 H_\beta y_\beta=T_\beta,\qquad
 a_\beta T_\beta=\Sigma_0,\qquad
 \sum_\beta\eta_\beta^2=1.                            \tag{14}
\]

Thus its selected principal covariance equals `Sigma0` without a transfer
from the distinct canonical family. The complete-curl remainder has an
additional margin `1/2` after the required epsilon division.

The initial and dynamic ledgers retain all six covariance entries, radial
and axial velocity, every curl remainder, `Pi1/Pi2`, compact pressure and the
whole derivative of the moving product `c_rho P`. The current post-signed
source has exact order

\[
 D={3\over2}-2\kappa,                                 \tag{15}
\]

and the constrained inverse retains both moving-normal and source-normal
pressure. The simultaneous next cycle reaches

\[
 B_1=2-5\kappa,\qquad C_1={5\over2}-5\kappa           \tag{16}
\]

with target slack `2*kappa`, relative to the named continuum inputs.

Q203 then retains time, unit viscosity, complete pressure, both nonlinear
cross orientations, the square, every scale-cutoff curl and pressure
derivative in the six-group diagonal residual. Its localized residual has
the exact twelve groups and the complete pressure-cutoff term.

## 6. Precise remaining limits

The correction in Section 4 means that unevaluated giant coefficient
seminorms prevent a numerical scale certificate but do not logically block
the rule-defined background or the uniform primary theorem.

Three limitations remain and are not interchangeable:

1. Equation (7) forbids deleting reserved pressure memory.
2. The available zero-jet physical upper-bound powers are negative:

   \[
   -s_x=-{201\over400},\qquad -2s_x=-{201\over200},    \tag{17}
   \]

   for primary velocity and pressure, while Q202 gives

   \[
   -A+hD=-{4975001\over10000000},\qquad
   -2A+h(D+1/2)=-{10000001\over10000000}.             \tag{18}
   \]

   These estimates do not certify absolute decay. They are upper-bound
   failures, not lower bounds on the physical pulses.
3. The complete sourced/signed hierarchy has fixed-jet estimates at every
   finite cycle, but its constants may depend arbitrarily on the cycle,
   derivative and requested power. No single uniform generated-hierarchy
   constant has been proved.

Consequently this checkpoint establishes the global ledger composition and
uniform primary family conditionally, but no unconditional complete pulse
hierarchy, global residual improvement, regular terminal origin, future
velocity/pressure, competitor exclusion or blow-up.

The next obligation is to audit the actual sourced and signed pulse
value/derivative/closed-support constants through every finite physical cycle
and determine whether the stage-dependent bounds suffice directly or need a
new hierarchy majorant. Complete pressure, axial velocity, covariance, unit
viscosity, cutoffs and every interaction must remain in that audit.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q204_global_ledger_uniform_primary/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q204_global_ledger_uniform_primary/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q204_global_ledger_uniform_primary/evidence.json \
  --output q204-fresh.json
```
