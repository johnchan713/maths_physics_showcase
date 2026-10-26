# Fixed-profile compensator and same-choice projector join

Parent: `5dfed175bd71a5c9537c71c8fde591115107cef5`.
Status: `fixed-profile-same-choice-projector-contract-joined`.

The parent closes the canonical full-path ODE, moving-normal pressure and
cutoff contract on one separately certified later tail. This checkpoint
audits the fourth primitive contract: the compact pressure compensator and
the two moment projectors on that same domain.

The audit corrects two inherited ambiguities. First, mass, support and a
second moment do not bound the derivatives of an arbitrary family of
profiles. The construction, however, chooses one smooth compact profile
once. Its fixed finite-order seminorms are finite and yield constants uniform
over the physical bands. No numerical formula for that profile is needed for
this asymptotic existence statement. Second, the repository did not identify
an explicit historical projector shape. We therefore make a new admissible
choice from the same fixed pressure profile and recompute every
projector-dependent target. This is not an equality with an earlier
unspecified projector or target.

## 1. Correct fixed-profile quantifier

Fix once and for all

\[
 \widehat\rho\in C_c^\infty((\ell,u)),\qquad
 \int \widehat\rho(x)\,dx=1,
 \qquad 0<\ell<u,
\]

inside the retained reserved annulus. On the existing scale chart put

\[
 \rho(R,Z,T)=s(Z,T)^{-1/2}
       \widehat\rho\!\left(R/\sqrt{s(Z,T)}\right),
 \qquad \tfrac14\le s\le4.                              \tag{1}
\]

For every fixed integer `m`, the finitely many seminorms
`H_j=||partial_x^j rhohat||_infinity`, `0<=j<=m`, are finite. The retained
scale-coordinate estimates and the finite chain rule then give a finite
constant

\[
 \max_{a+b+r\le m}\|\partial_Z^a\partial_T^b
             \partial_R^r\rho\|_\infty
 \le C_m(\ell,u;H_0,\ldots,H_m),                         \tag{2}
\]

uniformly in the band label. The same argument applies to the reference
density, the mismatch primitive and

\[
 c_\rho={1\over2}\int R^2\rho\,dR
        =s\,{1\over2}\int x^2\widehat\rho(x)\,dx.        \tag{3}
\]

The moving-density product calculation retains all six commutator terms and
the whole `epsilon partial_Z(Jz+c_rho P)`, including `c_rho,Z P`. Thus the
fixed profile supplies every finite compensator jet requested at this node.
Constants may depend on the profile and on `m`; there is no common bound over
all profiles, all orders or infinitely many construction cycles.

This quantifier is essential. The auxiliary families in the earlier gate
keep mass, support and second moment while their derivatives grow without
bound. Likewise a normalized moving projector may contain a factor
`1+(R-m_e)sin(NZ)/2`; its moment remains one but its `Z` derivative grows
like `N`. Those examples disprove a normalization-only inference. They do
not obstruct choosing and retaining one fixed shape as in (1).

## 2. Explicit same-profile projectors

Set

\[
 A_2=\int \widehat\rho(x)^2\,dx>0,
 \qquad
 \sigma_e(R,Z,T)=
 {R^{-e}\sqrt{s}\,\rho(R,Z,T)^2\over A_2},
 \quad e=1,2.                                            \tag{4}
\]

The positivity of `A_2` follows from unit mass. More quantitatively,
Cauchy--Schwarz on an interval of length `u-l` gives

\[
 1\le (u-\ell)A_2,
 \qquad A_2^{-1}\le u-\ell.                             \tag{5}
\]

Each projector is nonnegative and has the required exact weighted moment.
Indeed, with `R=sqrt(s)x`,

\[
 \int R^e\sigma_e(R,Z,T)\,dR
 =A_2^{-1}\int\widehat\rho(x)^2\,dx=1.                  \tag{6}
\]

The support at a fixed slow point is the actual moving interval

\[
       \sqrt{s}\,\ell\le R\le\sqrt{s}\,u.              \tag{7}
\]

Its common integration hull over the chart is `[ell/2,2u]`. The reserved
weight floor is used only on (7), where the field is supported, and not on
the rectangular product of the common hull with all slow variables. The
retained strict reserved margins keep (7) away from both radial edges.

Because `R` is bounded away from zero on this support, differentiating (4)
uses only finite Leibniz and chain-rule sums. For every fixed order `m`, the
existing provider gives a bound of the form

\[
 \|\sigma_e\|_m
 \le (u-\ell)\,B_{e,m}(\ell,u,s)\,
                  C_m(\widehat\rho)^2,                  \tag{8}
\]

where `B_e,m` contains only the already certified coordinate, dilation and
reciprocal-radial factors. Applying the existing logarithmic pullback and
the actual-support weight comparison once gives the required coefficient
bound in the Haar chart. Equations (2), (5) and (8) are finite and
band-uniform at each fixed order. They are profile-dependent tame bounds,
not numerical first-band evaluations.

## 3. Same-choice recomputation and support closure

Formula (4) is a new admissible choice. It cannot be inserted into cached
objects computed with an unspecified historical `sigma_e`. Starting from
the exact same pressure profile, the cited same-choice provider recomputes

\[
 h_e=E_e-\sigma_e M_e,
\]

so `int R^e h_e dR=0`, then recomputes the two-edge primitives, signed
amplitude and root input, curl, all six covariance entries and their full
pressure norms. The axial velocity channel and both background pressure
tails remain present. Old scalar defects are added to the new wave changes;
they are not replaced by them. All identity guards require this exact source,
projector family, root and pressure container.

The fixed smooth profiles have flat zero extensions, and the existing finite
support theorem propagates those extensions through the finite products,
cut potentials, full pressures, Haar averages and compact mean operations.
This is a qualitative finite-order support conclusion. It does not evaluate
the very high actual target seminorms or prove a smallness gain.

Together, (1)--(8) and the same-choice recomputation close the fourth
primitive contract on the parent's same accepted domain: the compensator and
projector profiles have their actual supports and finite fixed jets, both
moment identities are exact, and every downstream object is attached to the
same explicit choice.

## 4. Scope and next obligation

The exact controls in this directory use the finite-regularity polynomial
`30(x-1)^2(2-x)^2` only to verify rational normalization, dilation and
failure identities. It is not substituted for the selected smooth compact
profile. The analytical result is the fixed-profile composition argument
above. Passing the controls is not a continuum estimate for an actual giant
field.

The next obligation is to join the compact and slow inverse maps to these
same-choice coefficients on their complete weighted domains and whole-path
supports. Actual linked current-field high-jet bounds, strict epsilon
margins, fresh axial/full-pressure means, generated interactions, viscosity
and every cutoff term remain open. This checkpoint proves no uniformity over
profiles or orders, unconditional pulse family, admissible force, complete
residual improvement or blow-up theorem.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_fixed_profile_projector_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_fixed_profile_projector_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_fixed_profile_projector_join/evidence.json \
  --output /tmp/core_fixed_profile_projector_join_fresh.json
```
