# Same-tail canonical full-path ODE and pressure join

Parent: `382c7f90c4055b1a625d7a5bd6e90b4f112fccd4`.
Status: `same-tail-canonical-full-path-ode-pressure-and-cutoff-contract-joined`.

The parent put the complete phase normal, background values, spectral gaps
and stress cone on one accepted tail of the rule-defined Borel background.
This checkpoint audits the third primitive contract from the corrected
foundation frontier. On that same tail, it joins the original canonical
nonzero-entrance homogeneous primary, the zero-entry sourced inverse needed
by its pressure interface, and the original cutoff placement on the whole
lifted path.

The audit corrects an ambiguity in the inherited frontier wording. A
homogeneous equation has `f=0`, so it cannot itself retain a nonzero
unprojected source. The correct contract has two adjacent statements:
the canonical primary solves the homogeneous equation with its original
reference-frame entrance, while every sourced mode uses the same constrained
operator and keeps `n.f` in its pressure. Conflating them either invents a
source in the primary or deletes a real pressure term from later modes.

This is the separately certified later-tail construction. It does not
certify an unspecified historical first band or identify the canonical
reference-mixing primary with the inherited actual-mixing primary.

## 1. Full constrained equation and exact pressure

Fix a retained label and its complete lifted interval `0<=v<=L`. All
integer carriers, representatives, signs, length and entrance choices are
frozen before continuous differentiation. With the full summed background,

\[
 \Phi=p\theta+p_zZ/\epsilon+x_0R-v(pF+p_zG),\qquad n=\nabla_*\Phi,
\]

and `K=K_bg`, put `d=epsilon k^2 |n|^2`. For every nonzero harmonic `m`,
the original constrained system is

\[
 t'+Kt+m^2dt+ikm\,n\pi=-f,\qquad n\cdot t=0.             \tag{1}
\]

Differentiating the constraint and taking the scalar product of (1) with
`n` gives, without projecting the source first,

\[
 \pi={i\{n\cdot Kt-n'\cdot t+n\cdot f\}\over
             km|n|^2}.                                  \tag{2}
\]

Substitution yields

\[
 t'=A_nt-m^2dt-\Pi_{n^\perp}f,
 \quad A_n=-K+{n(n^TK-(n')^T)\over|n|^2}.               \tag{3}
\]

Conversely, (2)--(3) imply (1), and

\[
             (n\cdot t)'=-m^2d(n\cdot t).               \tag{4}
\]

Thus tangent entrance data stay tangent. A pure normal source has zero
projected velocity forcing but generally nonzero pressure; this is an exact
obstruction to replacing (2) by a projected-source pressure surrogate.
The moving-normal term is equally required. Exact rational controls in this
directory verify both directions of the algebra and detect either omission.

For the homogeneous primary set `f=0`, `|m|=1`, and use the original
canonical entrance

\[
 t(0)=P(0)B_r(0)e_+,
\]

where `B_r=U M_r(s_ref)`. This is not the inherited entrance in
`B_a=U M_a(s_actual)`. The canonical covariance and roots are therefore
rebuilt from this primary, as the existing canonical provider requires.

## 2. The parent tail supplies every denominator and frame input

On the parent's common tail,

\[
 |n|\ge {31\over32}B_s,\qquad B_s^2\ge\beta_{phase}>0, \tag{5}
\]

so `|n|^-2<=1024/(961 beta_phase)<4/beta_phase` throughout each complete
path and every accepted chart enlargement. The first contract supplies all
fixed background jets of the exact-real Borel field; differentiating the
full phase supplies `n,n',K,d` at every fixed order. No leading-only field
is substituted into these coefficients.

The original frame calculation on this same field gives

\[
 B^\ell(A_nB-B')=\operatorname{diag}(\lambda,-\lambda)+E,
 \quad \|E\|+|d-d_{ref}|\le H/L,                       \tag{6}
\]

with the accepted-tail length rows already included in the parent's one
integer selector. The tangent energy estimate consequently holds on every
subinterval:

\[
 \|V_m(v,w)\|\le T_0{P(v)\over P(w)}
  e^{-(m^2-1)d_*(v-w)},\qquad 0\le w\le v\le L.         \tag{7}
\]

The proof uses the complete moving frame, not a constant leading matrix.
The same derivative-independent tail works for all fixed jet orders because
the background Borel construction provides finite fixed-order seminorms;
orders change constants, not the value-domain selector.

For the Gaussian-adapted ambient propagator let `h=n/|n|^2`, `Hn=h n^T`,
`A=A_n-dI`, and

\[
             \widehat A=A+(d+P'/P)H_n.                 \tag{8}
\]

Then `Ahat=A` on tangent vectors and

\[
             (n')^T+n^T\widehat A=(P'/P)n^T.           \tag{9}
\]

This is why normal extensions grow with the same `P` ratio and do not
introduce a reciprocal Gaussian endpoint loss.

## 3. Varying growth and all fixed mixed jets

The actual reference rate is

\[
 {P'\over P}=\lambda_0(1+s_\sigma^2)^{-1/2}
 -{\lambda_0(1+s_\sigma^2)\over(1+u^2)^{3/2}},\qquad
 s_\sigma=\sigma u(1/2+v/L).                            \tag{10}
\]

It is not a frozen constant. Its third and higher derivatives are generally
nonzero. The finite inverse-square-root composition recurrence from the
phase-coefficient provider bounds every fixed jet in (10). The exact series
control here expands at a rational square point and verifies a nonzero
higher jet; setting all higher rate jets to zero fails that identity.

For the canonical nonzero entrance, the original derivative recurrence
retains jets of `P(0)B_r(0)e+`. Duhamel and (7) give

\[
 |D^It(v,x)|\le P(v)U_{|I|}(S),                         \tag{11}
\]

on the whole closed interval. The factor `P(0)` in the entrance and every
`P(w)` in an inhomogeneous term cancel only against the matching ratio in
(7); no lower bound on `P` is assumed. Pulse derivatives follow from the
ODE, including all coefficient derivatives. Formulae (2) and

\[
             {\cal A}={i(n\times t)\over km|n|^2}       \tag{12}
\]

then give fixed mixed pressure and potential jets. For sourced zero-entry
modes, the differentiated Duhamel recurrence retains derivatives of
`(m^2-1)d`; favorable damping at value order does not justify dropping
those terms. Smooth entrance-compatible sources give the same fixed mixed
jet induction.

Existence itself is finite-dimensional: continuous coefficients on the
closed interval make Picard differences bounded by `(ML)^j/j!`, so the
series converges uniformly for every fixed mode. The same estimate gives
uniqueness. No contraction smallness or source-size condition is needed.

## 4. Whole-path cutoffs and support

The raw canonical homogeneous solution is not asserted to be compactly
supported in `v`. Nor is a zero-entry forced response asserted to vanish at
the outgoing endpoint. The construction first solves and estimates the raw
equations on the entire interval, then multiplies the vector potential and
the full pressure by the original temporal and transverse cutoffs before
taking the physical operators:

\[
 \operatorname{curl}(\psi{\cal A})
 =\psi\operatorname{curl}{\cal A}+\nabla\psi\times{\cal A},
\quad
 \nabla(\psi\pi)=\psi\nabla\pi+\pi\nabla\psi.           \tag{13}
\]

Every derivative of the fixed cutoffs is retained. The same applies to the
transverse product cutoff in the canonical Haar columns. Flat cutoff zero
extensions and whole-path transverse source support give smooth support
extension by differentiated ODE uniqueness. A source that vanishes only at
the entrance does not give transverse support later, and an interior source
may leave a nonzero raw exit. These distinctions prevent a false support
claim.

The full pressure in (13) is (2), including `n.f`. Cutting only velocity,
or cutting a pressure reconstructed from the projected source, would miss
the second commutator in (13). Exact product-jet controls detect both lost
terms. These cutoff contributions, viscosity and nonlinear interactions
remain entries in the later residual ledger; this checkpoint establishes
their definitions and fixed-jet/support bounds, not their final smallness.

## 5. Closed contract and remaining frontier

Equations (1)--(13), combined with the parent's same-tail bounds, close the
third primitive contract: the canonical primary exists with its correct
entrance and complete moving-normal pressure on the whole path; the adjacent
sourced inverse keeps unprojected pressure; varying growth, cutoff
derivatives and support placement are retained on the same accepted tail.
This is an analytical composition of the cited provider estimates. The
finite controls verify algebra and provenance but are not continuum
sampling of the selected profile.

The next primitive contract is the fixed compensator/projector family:
prove its actual supports and fixed jets on this same domain, then join the
compact and slow inverse weighted-domain estimates. Generated harmonic
sources, complete covariance and fresh axial/full-pressure means, and the
complete viscosity/cutoff/interaction residual remain open. No
unconditional pulse family, admissible force, residual improvement or
blow-up theorem follows.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_canonical_full_path_ode_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_canonical_full_path_ode_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_canonical_full_path_ode_join/evidence.json \
  --output /tmp/core_canonical_full_path_ode_join_fresh.json
```
