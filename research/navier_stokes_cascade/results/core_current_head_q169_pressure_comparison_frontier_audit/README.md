# Q169: fixed-gauge pressure and expanding-ball comparison audit

Parent: `06185eb484f89e9b767f2221e30b4dc093b083f0` (Q168).
Status: `Q169-current-head-fixed-gauge-pressure-comparison-audited-conditionally`.

This checkpoint audits Q156's pressure recovery and expanding-ball comparison
on the exact Q164--Q168 localized field and force.  Q156 already occurs in
Q164's ancestry through Q163, but Q168's forward frontier did not explicitly
bind or rederive it.  The frozen evidence chain is

\[
 \mathrm{Q168}\to\mathrm{Q167}\to\mathrm{Q166}\to\mathrm{Q165}
 \to\mathrm{Q164}\to\mathrm{Q163}\to\mathrm{Q156},
\]

and Q163 and Q168 contain identical Q162 evidence, while Q156 and Q168
contain identical Q155 evidence.  Thus the field, force, supports,
localization and pressure gauge are the same.  Historical records are not
modified.

The result is conditional on the named continuum providers.  The executable
ledger verifies exact arithmetic, scope and provenance; it is not a
continuum proof of those providers.

## 1. Difference equation with full pressure memory

Fix \(T<1\).  Let \(u,p,f\) be the Q168 localized field, complete
fixed-gauge pressure and compact smooth force.  Let \(v,P\) be locally smooth
on \(\mathbb R^3\times[0,T]\), solve Navier--Stokes with the same force and
zero datum, and satisfy

\[
 \sup_{0\le t\le T}\|v(t)\|_2<\infty.                \tag{1}
\]

Put

\[
 w=v-u,\qquad \pi=P-p,\qquad
 g=v\otimes v-u\otimes u
   =w\otimes w+w\otimes u+u\otimes w.                \tag{2}
\]

The individual pressure equations retain the complete force divergence:

\[
\begin{aligned}
 -\Delta p&=\partial_i\partial_j(u_i u_j)-\nabla\cdot f,\\
 -\Delta P&=\partial_i\partial_j(v_i v_j)-\nabla\cdot f.
                                                               \tag{3}
\end{aligned}
\]

Subtracting cancels the common force exactly once:

\[
 w_t+\operatorname{div}g=\Delta w-\nabla\pi,
 \qquad \operatorname{div}w=0,
 \qquad-\Delta\pi=\partial_i\partial_jg_{ij}.        \tag{4}
\]

Neither individual pressure is replaced by an unforced Riesz formula, and
\(\nabla\cdot f\) is not subtracted twice.

## 2. Removing the harmonic pressure gradient

On \([0,T]\), \(w\) is uniformly \(L^2\) and \(g\) uniformly \(L^1\).  Set

\[
 \pi_*=\sum_{i,j}R_iR_jg_{ij}.                       \tag{5}
\]

The Fourier transform of an \(L^1\) stress is bounded.  Since
\((1+|\xi|^2)^{-2}\) is integrable in three dimensions,
\(\pi_*\in H^{-2}\).  This uses an exponent strictly above \(3/2\); the
critical endpoint is not invoked.

Poisson's equation alone only says that \(\pi-\pi_*\) is harmonic.  For a
compact smooth time test \(a\), the momentum equation gives

\[
 \int a\nabla\pi\,dt
 =\Delta\int aw\,dt+\int a'w\,dt
  -\operatorname{div}\int ag\,dt\in H^{-3}.          \tag{6}
\]

The same is true for \(\int a\nabla\pi_*\,dt\).  Their difference is a
harmonic \(H^{-3}\) distribution.  Its Fourier transform is a weighted
\(L^2\) function supported at the single point zero, and therefore vanishes.
Consequently

\[
 \nabla\pi=\nabla\pi_*                               \tag{7}
\]

as a space-time distribution.  A time-dependent spatial constant may remain
in \(\pi\), but its flux against \(w\nabla\chi_R\) is zero by
incompressibility.  No decay or temperedness assumption on \(P\) is needed.

The finite-energy hypothesis is essential.  Without it,
\(v=t e_1,P=-x_1\) and \(u=p=0\) have the same zero datum and force but a
nonzero harmonic pressure gradient.

## 3. Compact mixed pressure

Split

\[
 g=w\otimes w+h,
 \qquad h=w\otimes u+u\otimes w.                    \tag{8}
\]

The Q168 velocity has fixed compact support \(K\subset B_L\), so \(h\) is
supported in \(K\) and uniformly \(L^1\).  For
\(R\ge\max(1,2L)\) and \(R\le|x|\le2R\), kernel separation yields

\[
 |\pi_h(x,t)|\le C H_T R^{-3}.                       \tag{9}
\]

Choose \(\chi_R=\phi(x/R)^8\), where \(\phi=1\) on the unit ball and is
zero outside the radius-two ball.  Pairing (9) with
\(|\nabla\chi_R|\le CR^{-1}\), using the annular volume cost \(R^{3/2}\)
and the global \(L^2\) bound for \(w\), gives

\[
 \left|\int\pi_h w\cdot\nabla\chi_R\right|
 \le C_T R^{-3-1+3/2}=C_T R^{-5/2}.                 \tag{10}
\]

No isotropy or vanishing mixed moment is assumed.

## 4. Quadratic pressure commutator

Define

\[
 E_R=\int\chi_R|w|^2,
 \quad D_R^2=\int\chi_R|\nabla w|^2,
 \quad B_R=\|\phi_R^4w\|_6.
\]

Local Sobolev gives

\[
 B_R\le C(D_R+R^{-1}W_T).                            \tag{11}
\]

For \(\pi_q=\sum R_iR_j(w_iw_j)\), split

\[
 \phi_R^4\pi_q
 =\sum R_iR_j(\phi_R^4w_iw_j)
  +\sum[\phi_R^4,R_iR_j](w_iw_j).                   \tag{12}
\]

The commutator kernel is bounded by

\[
 C|z|^{-3}\min\{|z|/R,1\}.
\]

At exponent \(4/3\), its radial norm to the \(4/3\) power is

\[
 R^{-4/3}\int_0^Rr^{-2/3}\,dr
 +\int_R^\infty r^{-2}\,dr
 ={3\over R}+{1\over R}={4\over R}.                 \tag{13}
\]

Hence the commutator norm is \(O(R^{-3/4})\).  Riesz boundedness for the
first term and interpolation for the remaining factor give

\[
 \left|\int\pi_qw\cdot\nabla\chi_R\right|
 \le C_T\left(R^{-1}B_R^{3/2}
              +R^{-7/4}B_R^{3/4}\right).             \tag{14}
\]

The eighth cutoff power is minimal for this factorization:
\(k/2-1\ge3k/8\) is equivalent to \(k\ge8\).  The seventh power fails the
required weight comparison.

## 5. Expanding-ball estimate and legal limit

The complete weighted difference identity is

\[
\begin{aligned}
 {1\over2}E_R'+D_R^2={}&
 -\int\chi_R(w\cdot\nabla)u\cdot w
 +{1\over2}\int|w|^2\Delta\chi_R\\
 &+{1\over2}\int|w|^2v\cdot\nabla\chi_R
 +\int\pi w\cdot\nabla\chi_R.                       \tag{15}
\end{aligned}
\]

For sufficiently large \(R\), \(u=0\) on the cutoff annulus, so \(v=w\)
there.  This support separation is made before estimating the transport
flux.  After substituting (11), the two nonconstant dissipation/radius rows
are

| term | \(D_R\) power \(a\) | radius power \(b\) | Young decay | Sobolev offset |
|---|---:|---:|---:|---:|
| quadratic pressure and transport | \(3/2\) | \(1\) | \(4\) | \(5/2\) |
| pressure commutator | \(3/4\) | \(7/4\) | \(14/5\) | \(5/2\) |

Here Young's inequality sends \(R^{-b}D_R^a\) to a dissipation term plus
\(O(R^{-2b/(2-a)})\).  Together with (10), every such error decays faster
than \(R^{-2}\).  The coarse cutoff-Laplacian bound is the limiting term:

\[
 \boxed{
 {1\over2}E_R'+{1\over2}D_R^2
 \le\|\nabla u(t)\|_\infty E_R+C_TR^{-2},
 \qquad E_R(0)=0.}                                   \tag{16}
\]

No derivative, dissipation or boundary-flux bound for the competitor enters
(16).  The \(R^{-2}\) rate is the result of this coarse uniform-\(L^2\)
estimate, not a universal optimal rate for each fixed solution.

Let \(L_T=\sup_{[0,T]}\|\nabla u\|_\infty<\infty\).  Gronwall gives

\[
 E_R(t)\le2C_TT\exp(2TL_T)R^{-2}.                   \tag{17}
\]

Fix \(T<1\) before sending \(R\to\infty\).  Every fixed ball eventually
lies where \(\chi_R=1\), so (17) forces \(w=0\) there.  Since \(T<1\) is
arbitrary,

\[
 v=u\quad\hbox{on }\mathbb R^3\times[0,1).          \tag{18}
\]

The constants in (17) need not remain bounded as \(T\uparrow1\), and the
cutoff family need not be monotone.

## 6. Scope and next obligation

Q169 closes fixed-gauge pressure recovery and preterminal comparison
conditionally on the current object.  It does not use the protected-core
lower bound, take a terminal contradiction, construct future velocity or
pressure, or prove unconditional exclusion or blow-up.

The next justified checkpoint is to audit Q157/Q163's quantitative physical
core path, localization plateau and order of limits on the Q164--Q169 object.
For each fixed protected \(\tau\), that audit must choose \(T_\tau<1\), finish
the radius limit at fixed \(T_\tau\), and only then let \(\tau\downarrow0\).

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q169_pressure_comparison_frontier_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q169_pressure_comparison_frontier_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q169_pressure_comparison_frontier_audit/evidence.json \
  --output q169-fresh.json
```
