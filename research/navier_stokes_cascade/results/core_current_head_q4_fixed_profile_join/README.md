# Current-head Q4 fixed profile, compensator and same-choice projectors

Parent: `8d9b31a46eb8ab13028e338fdabc0fc6ece8ec2a`.
Status: `current-head-Q4-fixed-profile-compensator-and-same-choice-projectors-bound`.

Q1 fixes the exact core, roots, repair and leading pressure. Q2 supplies the
regional and factored-edge margins on an explicit later tail. Q3 binds that
same object to the complete frame, whole lifted path, phase jets and
constrained pressure ODE. Q4 asks whether the compact pressure compensator
and moment projectors can be selected on this same field with fixed-order
band-uniform jets, exact moments and compatible supports.

They can. The conclusion uses one fixed smooth compact profile and a fresh
same-choice recomputation. It corrects two inherited shortcuts: moments and
support do not control derivatives of an arbitrary profile family, and a
moving interval enclosing a scaled compact support is not the support
itself.

## 1. Fixed-profile quantifier and complete moving pressure

Choose once and for all

\[
 \widehat\rho\in C_c^\infty((\ell,u)),\qquad
 \int\widehat\rho(x)\,dx=1,
\]

inside the reserved annulus, and set on the unchanged scale chart

\[
 \rho(R,Z,T)=s(Z,T)^{-1/2}
  \widehat\rho\!\left(R/\sqrt{s(Z,T)}\right),
 \qquad \tfrac14\le s\le4.                              \tag{1}
\]

For every fixed `m`, the finitely many seminorms
`H_j=||partial_x^j rhohat||_infinity`, `j<=m`, are finite. Q3 supplies the
fixed scale/frame coordinate jets and positive denominators on the closed
accepted paths. Finite chain and product rules therefore give

\[
 \max_{a+b+r\le m}\|\partial_Z^a\partial_T^b
    \partial_R^r\rho\|_\infty
 \le C_m(\ell,u;H_0,\ldots,H_m)<\infty,                 \tag{2}
\]

uniformly over the accepted later bands at that fixed order. The constants
may depend on the profile and `m`; no uniformity over arbitrary profiles,
all orders or infinitely many cycles is asserted.

The same argument applies to the reference density, mismatch primitive and

\[
 c_\rho={1\over2}\int R^2\rho\,dR
       ={s\over2}\int x^2\widehat\rho(x)\,dx.           \tag{3}
\]

The complete moving-density calculation retains all six product/commutator
terms. In particular the axial pressure moment remains

\[
 \epsilon\partial_Z(J_z+c_\rho P)
 =\epsilon\{J_{z,Z}+c_{\rho,Z}P+c_\rho P_Z\}.           \tag{4}
\]

The middle term is generally nonzero. Freezing `c_rho`, replacing total
pressure by only its new change, or dropping either old `Pi1/Pi2` pressure
tail changes the operator.

## 2. Explicit projectors and their actual supports

Let

\[
 A_2=\int\widehat\rho(x)^2\,dx,
 \qquad
 \sigma_e(R,Z,T)={R^{-e}\sqrt{s}\,\rho(R,Z,T)^2\over A_2},
 \quad e=1,2.                                           \tag{5}
\]

Unit mass and Cauchy--Schwarz give

\[
 A_2>0,\qquad 1\le(u-\ell)A_2,
 \qquad A_2^{-1}\le u-\ell.                            \tag{6}
\]

With `R=sqrt(s)x`, both required moments are exact:

\[
 \int R^e\sigma_e(R,Z,T)\,dR=1.                        \tag{7}
\]

If `K=supp(rhohat) \Subset (ell,u)`, then

\[
 \operatorname{supp}\sigma_e=\sqrt{s}\,K
 \Subset(\sqrt{s}\ell,\sqrt{s}u).                     \tag{8}
\]

The interval in (8) and the common hull `[ell/2,2u]` are enclosures. The
reserved weight floor is applied on the actual moving support, where the
strict Q2 margins keep it away from both radial edges. Since `R` is bounded
away from zero there, every fixed derivative of (5) is a finite Leibniz and
chain-rule sum. The logarithmic pullback and actual-support weight comparison
are each applied once, yielding finite band-uniform fixed-order projector
constants.

The fixed quantifier is essential. One can preserve mass, support and the
weighted moment while inserting an oscillatory factor whose `Z` derivative
grows like an arbitrary frequency. Thus normalization alone cannot prove a
projector jet bound; equations (1)--(2) do so because the shape is fixed.

## 3. Same-choice recomputation

Formula (5) is an explicit new admissible choice. It is not identified with
an unspecified historical projector. Starting with the unchanged Q1--Q3
field and pressure profile, the cited same-choice branch recomputes

\[
 h_e=E_e-\sigma_e M_e,
 \qquad \int R^e h_e\,dR=0,                             \tag{9}
\]

then recomputes both edge primitives, the signed amplitude and root input,
the complete cylindrical curl, all six covariance entries and their full
pressure norms. The axial velocity, current whole-path wave pressure,
moving compensator and both background pressure tails remain in the same
container. Old scalar defects are added to the new wave changes rather than
replaced by them.

The fixed profiles vanish flatly outside their actual supports. Under the
retained weighted-input contract, the finite support theorem propagates
these zero extensions through finite products, cut potential and pressure,
Haar averages and compact mean operations. This is a qualitative
fixed-order support result, not a numerical evaluation of the very high
target seminorms.

## 4. Result and remaining frontier

Equations (1)--(9), Q3's whole-path coefficient bounds, the reserved support
margins and fresh same-choice recomputation bind the fourth primitive
contract to the explicit Q1--Q3 later-tail object. The exact rational
controls verify the moment, support and dependency interfaces; they do not
certify continuum bounds by computation and their finite polynomial profile
is not substituted for `rhohat`.

The next obligation is to join the shifted compact maps, zero-Haar temporal
inverse and fixed five-profile slow response to these exact coefficients on
their complete weighted domains and whole-path supports. Operator
boundedness will not by itself prove the actual generated-source inequality

\[
 \sup_{\ell\ge\ell_*}\epsilon_\ell^{-b_j}
 \max_{|\alpha|\le R_j(m,N)}\sup_{\Omega_{j,\ell}}
 \zeta^{-1}\delta^{p_{j,\alpha}}|D^\alpha F_{j,\ell}|<\infty. \tag{10}
\]

That input, complete covariance and fresh axial/full-pressure means, and the
full viscosity/cutoff/interaction residual remain downstream. No
unconditional pulse, force, exclusion or blow-up result is claimed.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q4_fixed_profile_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q4_fixed_profile_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q4_fixed_profile_join/evidence.json \
  --output fresh.json
```
