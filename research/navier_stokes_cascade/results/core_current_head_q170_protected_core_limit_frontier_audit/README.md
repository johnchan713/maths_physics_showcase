# Q170: protected physical core and terminal-limit audit

Parent: `388c35b51044b14738b381d86fd1532a0783a7fd` (Q169).
Status: `Q170-current-head-protected-core-and-limit-order-audited-conditionally`.

This checkpoint audits Q157/Q163's protected core path and legal terminal
limit on the exact Q164--Q169 object.  The evidence chain runs from Q169
through Q168--Q164 to Q163, and Q163 contains the exact Q157 evidence.  Q157
and Q168 share Q155 byte for byte; Q157 and Q169 share Q156 byte for byte.
The field, Borel background, physical diagonal, force, localization, support
hole and pressure gauge are therefore the same.  Historical evidence is not
modified.

The conclusion remains conditional on the named continuum providers.  Tests
check exact arithmetic, scope and provenance; they do not replace those
analytical premises.

## 1. Fixed physical path and positive leading coefficient

Keep

\[
 h={1\over200},\qquad A={1\over2}+h={101\over200},
 \qquad X_a={4\over101},\qquad X_{\rm in}={1\over101}.
                                                               \tag{1}
\]

The support-hole margin is exactly

\[
 X_a-X_{\rm in}={3\over101}>0.                         \tag{2}
\]

For \(\tau=1-t>0\), choose

\[
 x_\tau=(\sqrt{2X_{\rm in}\tau},0,0).                 \tag{3}
\]

At \(z=0\), the implicit chart gives

\[
 q=\tau,\qquad\eta=0,\qquad X=X_{\rm in},
 \qquad r^2=2X_{\rm in}\tau.                         \tag{4}
\]

On this path the radial term in the Cartesian core ansatz has zero
azimuthal component.  The leading normalized swirl is

\[
 q^Au_{\theta,0}(x_\tau,1-\tau)
 =E_0(X_{\rm in},0)=:e_0,
\]

where

\[
 e_0=\sqrt{2X_{\rm in}}{\varphi_*\Phi\over C}>0.      \tag{5}
\]

Positivity follows from the retained exact-core lower bound
\(\Phi>33/125\), the positive exponential \(\varphi_*\), and \(C>0\).  This
is a velocity coefficient identity; stress-growth diagnostics are not used
as a substitute.

## 2. Exact annular vanishing versus background smallness

The two correction mechanisms must remain distinct.

Every fixed-stage physical potential, direct swirl, pressure difference and
synchronized stage-zero annular term has closed support in

\[
 X_a\le X\le X_b.                                    \tag{6}
\]

Multiplication by the scale cutoff, differentiation, curl, and the cutoff
commutator do not enlarge this support.  Since \(X_{\rm in}<X_a\), every
annular diagonal correction vanishes identically in a neighborhood of (3).

The positive-order Borel background does not vanish.  Q149 binds the
physical primary to the same finalized leading profile and rule-defined
background used by the core provider.  Its finite-head/geometric-tail
estimate gives

\[
 \tau^Au_\theta(x_\tau,1-\tau)=e_0+R_{\rm bg}(\tau),
 \qquad |R_{\rm bg}(\tau)|\le C_{\rm bg}\tau^{2h}.    \tag{7}
\]

No uniform-in-coefficient-order constant or canonical/physical pulse
identification is inserted.

## 3. Exact positive threshold

Let \(q_0,C_{\rm bg}>0\) be the range and constant in (7).  Since
\(1/(2h)=100\), define

\[
 \rho=\min\left\{1,{e_0\over2C_{\rm bg}}\right\},
 \qquad q_{\rm core}=\min\{q_0,\rho^{100}\}.          \tag{8}
\]

For \(0<\tau\le q_{\rm core}\), monotonicity gives

\[
 C_{\rm bg}\tau^{1/100}
 \le C_{\rm bg}\rho\le {e_0\over2}.                 \tag{9}
\]

Equation (9) is the exact estimate.  The executable ledger records it using
the rational base \(\rho\); it does not evaluate a rational number to a
fractional power numerically.

Let \(\tau_0>0\) be the time-cutoff width and \(r_L>0\) the spatial plateau
radius.  Put

\[
 \tau_*=min\left\{q_{\rm core},{\tau_0\over2},
                  {r_L^2\over8X_{\rm in}}\right\}>0. \tag{10}
\]

The last restriction yields

\[
 r_\tau^2=2X_{\rm in}\tau\le {r_L^2\over4},          \tag{11}
\]

while the time restriction places the path on the late unit plateau.  Thus
\(c=1\), \(\nabla c=0\), and localization preserves (7).  Combining
(7)--(11) gives

\[
 \boxed{
 |u_\theta(x_\tau,1-\tau)|
 \ge {e_0\over2}\tau^{-101/200}}
 \qquad(0<\tau\le\tau_*).                            \tag{12}
\]

## 4. Fixed-time, radius, then terminal limit

Let \(v,P\) satisfy the ordinary Q169 comparison hypotheses: local
smoothness, the same compact smooth force and zero datum, and a uniform
\(L^2\) bound on every fixed preterminal interval.  For each fixed
\(0<\tau\le\tau_*\), set

\[
 t_\tau=1-\tau,
 \qquad T_\tau=1-{\tau\over2}.                       \tag{13}
\]

Then

\[
 0<t_\tau<T_\tau<1.                                  \tag{14}
\]

First fix \(T_\tau\) and recover the pressure gradient.  Next send the
comparison radius \(R\) to infinity at this fixed time horizon.  Q169 then
gives the exact identity

\[
 v(x_\tau,t_\tau)=u(x_\tau,t_\tau).                  \tag{15}
\]

Only after the radius limit is finished is \(\tau\downarrow0\) allowed.
No bound on comparison constants uniform as \(T\uparrow1\) is needed.

The path tends to \((0,1)\), while (12) diverges because
\(101/200>0\).  A globally smooth competitor would have bounded velocity on
one compact neighborhood of \((0,1)\), contradicting (12) and (15).
Therefore, relative to every named continuum provider, no globally smooth
same-force zero-datum competitor in the fixed-interval finite-energy class
exists.

## 5. Scope and next obligation

This is a conditional same-force competitor-exclusion implication.  It does
not construct velocity or pressure after time one, prove that every retained
provider theorem is independently discharged, establish an unconditional
admissible force, or prove Navier--Stokes blow-up.  Dependency hashes and
passing tests certify interfaces and regressions only.

Q163's old next step was an independent core audit; Q164 has already
completed that checkpoint, so repeating it would be stale.  The next
justified obligation is a consolidated provider-discharge audit of Q164--Q170.
It must reconcile the older comparison-foundation unresolved flag with the
later independent audits and prove, or state the precise obstruction to, the
first continuum premise still not independently discharged.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q170_protected_core_limit_frontier_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q170_protected_core_limit_frontier_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q170_protected_core_limit_frontier_audit/evidence.json \
  --output q170-fresh.json
```
