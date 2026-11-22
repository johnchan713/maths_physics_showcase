# Q210: finite primitive unification and the uniform all-band obstruction

Parent: Q209, `218ba251c6a9c79d56623da304db1e28784d6ad8`.
The live branch, Q209, its active CI, Q207's stage-602 trace, Q165, Q184,
Q196, the weighted-stress and frame/gap records, Q189, Q192 and the direct
compact/slow reconstruction were read first. Historical evidence, the
selected profile, roots, pressure gauge, supports, cutoff sequences and the
unchanged `C^-20` activation width remain fixed.

Q209 correctly separates the background and generated-stage cutoffs. Its
remaining frontier nevertheless phrases a separately displayed primitive
majorant `A_R`, gap minimum and inverse archive as the first obligation. The
live all-fixed-order contracts already imply that finite bundle for the
concrete stage-602 request, relative to their named continuum premises. No
numerical materialization is needed. What they do not independently prove is
the normalized uniform supremum over the infinite accepted-band tail. This
checkpoint proves the finite unification, gives a smooth exact countermodel
to the invalid per-band-to-uniform inference, and isolates the surviving
inequality. This narrows Q209's frontier; it does not contradict Q209's own
warning that no numerical high-jet norm had been evaluated.

## 1. A finite request has one finite conditional majorant

Q207's retained first request is

\[
 j=602,\qquad R=36\,479\,350\,300,
 \qquad q_{pole}=218\,876\,101\,818.                 \tag{1}
\]

Q184's corrected source contract says that for each fixed cycle, harmonic
and multi-index `I` there are finite `C_I,d_I,r_I,s_I` such that

\[
 \left|D^I(f_{j,m}/P)\right|
 \le C_I\epsilon^{B_j}S^{d_I}\sqrt\zeta\,
             y_a^{-r_I}y_b^{-s_I},\qquad S=\ell^2.   \tag{2}
\]

Q196 proves the corresponding initial-source row with its corrected order,
and Q192 propagates any supplied primitive array through the finite frame
graph without changing the field. The set of multi-indices with
`|I|<=R+1` in the fixed finite coordinate dimension is finite. Therefore,
conditional on (2), one may define

\[
 A_R=\max_{|I|\le R+1} C_I<\infty,                   \tag{3}
\]

while retaining the individual polynomial and edge weights
`S^-d_I y_a^r_I y_b^s_I`. The maximum may be enormous. Finiteness of a
maximum over a finite set does not require listing its entries or evaluating
their binary expansions.

The same finite-set argument applies to the other two parts of Q207's
conditional recipe. The same-band and leading-frame records give positive
value, stress-direction, normal and determinant floors on one accepted tail;
Q189 supplies the applicable compact-repair gaps on the same full-pressure
profile; the finite stage prefix uses finitely many accepted root and support
gaps. Their minimum is positive, conditional on those continuum records.
The direct compact/slow reconstruction proves a finite operator constant for
every fixed derivative request and pole, including both flat edges, the
physical characteristic shift, Fourier divisors, axial velocity and the
whole moving pressure product. The finite request archive therefore has a
finite maximum `K_R`.

Thus the existing provider statements imply

\[
 C_{602,0,0}
 =F_{602,0,0}(1+A_R,1+g_R^{-1},1+K_R)<\infty         \tag{4}
\]

without numerical evaluation. This is still a conditional theorem: (3)--(4)
do not independently prove (2), the gap premises or the operator inputs from
the raw selected maps.

## 2. Per-band smoothness does not prove the uniform tail

The unresolved quantifier in (2) is not finiteness on any individual band.
For a fixed nonzero `phi in C_c^infinity(Omega)` consider the scalar family

\[
 F_\ell=\epsilon_\ell^B e^{\sqrt\ell}\phi,
 \qquad \epsilon_\ell=2^{-h\ell},\quad B>0.          \tag{5}
\]

Every band is smooth, has the same compact interior support, and has all
endpoint jets equal to zero. Its unnormalized size tends to zero because

\[
 -hB(\log2)\ell+\sqrt\ell\longrightarrow-\infty.    \tag{6}
\]

Nevertheless, for every fixed polynomial degree `d`, the normalized size is

\[
 \epsilon_\ell^{-B}S^{-d}\|F_\ell\|
 =e^{\sqrt\ell}\ell^{-2d}\|\phi\|.                 \tag{7}
\]

On the integer subsequence `ell=n^4`, its logarithm is

\[
 n^2-8d\log n\ge n^2-4dn=n(n-4d)\longrightarrow\infty, \tag{8}
\]

using the elementary bound `log n<=n/2`. Hence no fixed polynomial degree
makes (7) uniformly bounded. This is an analytical counterexample to an
inference, not a counterexample to the actual selected Navier--Stokes source.
It also shows why absolute decay of every band does not repair the missing
normalized estimate.

## 3. Exact surviving inequality

The first independent source obligation is therefore, for every fixed
stage/harmonic/jet requested by the finite recurrence, a proof from the
explicit selected maps that some finite polynomial degree exists with

\[
 \boxed{
 \sup_{\ell\ge\ell_*}\epsilon_\ell^{-B_j}S^{-d_{j,I}}
 \sup_{\Omega_{j,\ell}}
 \zeta^{-1/2}y_a^{r_{j,I}}y_b^{s_{j,I}}
 \left|D^I(f_{j,m}/P)\right|<\infty .}              \tag{9}
\]

Equation (9) must use the same selected profile, roots, frame and gauges. It
must retain radial and axial velocity, all six covariances, `Pi1/Pi2`, the
whole moving pressure product, viscosity, cylindrical connections, cutoff
curls, both nonlinear orientations, the square and every supported row.
Q192 can propagate a supplied polynomial-in-`S` primitive majorant, and the
compact/slow theorem maps the resulting weighted source space boundedly; an
operator estimate cannot prove its own input (9).

Q184 and Q196 derive (9) conditionally from the named profile/frame
contracts. This checkpoint does not promote those contracts by source hashes
or tests. It corrects the requested action: evaluating `36,479,350,301`
derivatives or a giant cutoff number is not the mathematical next step.
The next justified step is a raw-map derivation of (9), beginning with the
explicit primary-plus-initial-mean source and then one sourced/generic cycle.

Q210 proves no numerical stage-602 constant, unconditional pulse hierarchy,
complete residual, admissible force, competitor exclusion or Navier--Stokes
blow-up.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q210_finite_bundle_uniform_tail/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q210_finite_bundle_uniform_tail/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q210_finite_bundle_uniform_tail/evidence.json \
  --output fresh.json
```
