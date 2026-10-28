# Current-head Q3 full frame, whole path and phase domain

Parent: `c977837a780f6830590d4dca4814cbc171ffb7b5`.
Status: `current-head-Q3-full-frame-whole-path-and-phase-domain-bound`.

Q1 fixes the exact pressure-family core, roots, repair and leading pressure.
Q2 puts its value, determinant, orientation, cone and factored-edge margins
on one explicit later tail. Q3 asks whether the original complete frame,
phase jets and whole-path primitive estimates hold on that same object and
physical domain with finite continuum constants.

The answer is affirmative on the Q2 later tail, with one correction. The
unverified numerical comparison `36 Cphi^2<=beta S^2` belongs to an older
sufficient normal-floor route. It is not a necessary condition for the
original full-frame proof. Q2 already supplies a positive reference-relative
normal floor on the complete enlarged charts. This checkpoint uses that
floor and derives the remaining finiteness from fixed-order Borel seminorms
and compact whole-path domains; it does not certify an unspecified historical
first band.

## 1. Compactness supplies the original finite frame constants

Fix a retained label and its closed lifted interval `0<=v<=L<=S`. The Q2
tail has finitely many regional/chart types and a derivative-independent
accepted value domain. For every fixed order, the real Borel construction
gives finite background seminorms before its cutoff scale is selected. Hence
the full phase coefficients

\[
 n=\nabla_*\Phi,\qquad n',\qquad K=K_{bg},\qquad
 d=\epsilon k^2|n|^2,
\]

and all their fixed mixed jets are finite on the closed path. Both axial
transport terms in `F,G` remain in `n,n'`; no leading-only phase is used.

Q2 supplies positive normal, shear, carrier/frame and canonical-root floors
on the same tail. Every frame, left inverse and phase coefficient is therefore
a smooth rational expression with denominators bounded away from zero. In
particular the original constants may be defined by the finite suprema

\[
 C_E=\sup L\|E\|,\qquad C_d=\sup L|d-d_{ref}|.          \tag{1}
\]

The endpoint frame cost is finite because `cmin>0` and the normal floors are
positive. Thus

\[
 H=C_E+C_d<\infty,qquad
 K_{frame}={(1+S_0)\max(1,c_{max})\over\min(1,c_{min})}<\infty. \tag{2}
\]

No small numerical upper bound is needed. The tangent propagator estimate

\[
 \|V_m(v,w)\|\le K_{frame}e^H{P(v)\over P(w)}
 e^{-(m^2-1)d_*(v-w)},qquad d_*={b_{min}^2\over4}>0,  \tag{3}
\]

is finite on every subinterval. The positive profile `P` is an exponential
of the finite varying rate; its endpoint and Duhamel factors cancel through
the displayed ratio, so no artificial lower bound on `P` is introduced.

This closes the Q3 finiteness premise analytically. Rational fixtures test
the interface formulas, but they are not numerical enclosures of the giant
selected profile.

## 2. Complete constrained ODE and pressure

For each nonzero harmonic, retain the original system

\[
 t'+Kt+m^2dt+ikm\,n\pi=-f,\qquad n\cdot t=0.           \tag{4}
\]

Differentiating the constraint gives

\[
 \pi={i\{n\cdot Kt-n'\cdot t+n\cdot f\}
          \over km|n|^2}.                              \tag{5}
\]

The homogeneous canonical primary has `f=0` and its reference-frame
entrance. The adjacent sourced inverse uses the same operator and retains
the unprojected `n dot f`. A purely normal source can have zero projected
velocity forcing and nonzero pressure, so projecting before (5) is invalid.
The moving-normal term is independently nonzero and cannot be dropped.

Continuous finite coefficients on `[0,L]` give Picard differences bounded
by `(ML)^j/j!`; their ratio is `ML/(j+1)`, eventually at most `1/2` for every
finite `M,L`. Existence, uniqueness and each fixed differentiated jet need
no contraction smallness and no derivative-dependent tail. Complete
pressure and potential jets inherit the same positive denominators.

## 3. Whole-path support and cutoff placement

The raw homogeneous pulse is solved on the entire path and is not declared
compact at the outgoing endpoint. A zero-entry forced response may also be
nonzero there. Temporal and transverse cutoffs are applied to the potential
and full pressure only after the raw solve:

\[
 \operatorname{curl}(\psi A)=\psi\operatorname{curl}A+
 \nabla\psi\times A,
\]

\[
 \nabla(\psi\pi)=\psi\nabla\pi+\pi\nabla\psi.          \tag{6}
\]

Both commutators and every fixed cutoff derivative remain. Smooth transverse
zero extension follows only when the source vanishes along the whole path;
entrance flatness alone is insufficient.

## 4. Result and remaining frontier

The Q1-Q2 field, fixed-order Borel seminorms, positive denominator/root
floors, original frame formulas and canonical constrained ODE now form one
acyclic same-object Q3 chain on the explicit later tail. This binds the
full-frame primitive/path-domain premise without the obsolete numerical
`Cphi` gate or an unknown historical first band. It is an analytical
composition of the cited continuum providers, not an independent reproof of
every provider and not a proof by tests.

The next obligations are the fixed compensator/projector supports and jets,
then the compact and slow inverse weighted domains on this same field.
Generated harmonic sources, complete-wave covariance and fresh axial/full-
pressure means, and the entire viscosity/cutoff/interaction residual remain
downstream. No unconditional pulse, force, exclusion or blow-up result is
claimed.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q3_full_frame_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q3_full_frame_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q3_full_frame_join/evidence.json \
  --output fresh.json
```
