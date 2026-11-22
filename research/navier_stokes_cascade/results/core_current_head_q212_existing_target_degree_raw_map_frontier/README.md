# Q212: existing target degree contract and the restored raw-map frontier

Parent: Q211, `63c6394c8c27520357514e5b9ab7d2fc1da8d0b5`.
The live branch, Q211 and its active CI, Q149/Q150/Q152, Q160, Q180, Q184,
the finite-stage class definition, Q192, Q209/Q210 and their frozen evidence
were read first. Historical evidence, the selected profile, roots, pressure
gauge, supports, cutoffs and unchanged `C^-20` activation width remain fixed.

Q211 correctly closes the finite polynomial degree algebra, but its proposed
next step asks again for the first same-family weighted current-target degree
array. Q160 already derives that exact state-difference theorem, and the
finite-stage definition states that membership in its `W/M/S` classes
includes a finite polynomial degree in `S=ell^2` at every fixed stage and
jet. Q180 corrects the sequential reading by coupling target and source, and
Q184 retains the polynomial loss explicitly. This checkpoint joins those
records and restores the independent frontier to the raw selected-map
premises on which they depend.

## 1. Restore the suppressed polynomial factor

Q160 displays, for every finite cycle `j` and fixed mixed jet `I`,

\[
 |D^I\Xi_j|\le C_{j,I}\epsilon^{C_j}\zeta
 y_a^{-p_{j,I}}y_b^{-q_{j,I}}.                       \tag{1}
\]

Read alone, (1) appears to assign degree zero in `S`. That is not the class
definition used by the proof. The finite-stage record defines `W/M/S` to
allow constants and polynomial `S` degrees depending on the stage and
requested derivative while remaining uniform over physical bands and labels.
Q184 later makes the same convention explicit when correcting the generated
source normalization. Thus the complete statement represented by (1) is

\[
 |D^I\Xi_j|\le C_{j,I}\epsilon^{C_j}S^{d_{j,I}}\zeta
 y_a^{-p_{j,I}}y_b^{-q_{j,I}},\qquad d_{j,I}<\infty. \tag{2}
\]

No value for `d_(j,I)` is supplied. In particular, this checkpoint does not
promote the displayed shorthand (1) to a degree-zero theorem.

For a finite derivative request `R`, the index set is finite, so

\[
 D_R^\Xi:=\max_{|I|\le R}d_{j,I}<\infty.             \tag{3}
\]

This is an existential maximum, not a numerical evaluation. The same finite
maximum argument applies to the finitely many root, frame and mean rows. Q192
then propagates the supplied degree array through products, positive-floor
reciprocals and the differentiated Duhamel graph.

## 2. The current-target primitive is already in that class

Q160 takes the cancellation before estimating:

\[
 X_e^{new}-X_e^{old}=-I_eQ_e\overline{\Delta E_e}.   \tag{4}
\]

Its pressure reconstruction retains all six covariance entries. It also
retains radial and axial velocity, compact pressure, `Pi1/Pi2` in the same
base equation, the derivative of the whole moving pressure product and both
complete tangential moment identities. Every physical wave product carries
`zeta`; current mean products preserve it because `zeta^2<=zeta`.

The inputs of the two oriented current inverses have exact zero weighted
totals. Consequently both edge primitives retain `zeta`, with finite
polynomial-`S` and inverse-edge costs. Reserved corrections remain on the
fixed interior support where `zeta` has a positive floor. Q180 performs this
step simultaneously with the generated source, avoiding the circular
sequential reading of Q160. Therefore the target/primitive part of Q211's
`D_R^prim<infinity` obligation is already closed **relative to the named
continuum providers**.

This conclusion is conditional. Source hashes and executable tests establish
that the same records are being joined; they do not prove their continuum
suprema.

## 3. Exact surviving obstruction

The conditional target theorem rests on the selected leading/background
geometry and finite existence, gap and inverse premises. The first
independent same-object obligation is again the Q209 raw-map inequality

\[
 \boxed{
 \max_{|I|\le R+1}\sup_{\mathcal D}
 w_I S^{-d_I}
 \left|D^I(F_0,U_0,\Pi_0,M_0,V_0,
                  E_1,U_1,W_1,T_1/\zeta)\right|
 \le A_R<\infty .}                                  \tag{5}
\]

It must hold on one actual accepted tail for the selected profile and first
coefficient, with complete pressure memory. On that same object one must also
prove

\[
 g_R:=\min\{\hbox{Jacobian, normal, root, determinant, variance,
 support gaps}\}>0,                                 \tag{6}
\]

and a finite maximum `K_R` of every requested compact/slow pressure and
moment inverse constant. Q210 conditionally unifies these finite rows; it
does not prove (5)--(6) from the raw selected maps. Per-band smoothness does
not imply their normalized all-band form.

Q212 therefore corrects Q211's next step without changing its degree algebra.
It establishes no evaluated high-jet degree array, unconditional actual-pulse
hierarchy, complete residual, admissible force, competitor exclusion or
Navier--Stokes blow-up. The next obligation is to prove (5), (6) and the
same-object inverse archive directly for the selected profile/first
coefficient family.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q212_existing_target_degree_raw_map_frontier/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q212_existing_target_degree_raw_map_frontier/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q212_existing_target_degree_raw_map_frontier/evidence.json \
  --output fresh.json
```
