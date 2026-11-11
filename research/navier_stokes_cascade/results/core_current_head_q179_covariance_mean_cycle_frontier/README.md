# Q179: complete initial covariance means and finite-cycle frontier

Parent: Q178, `d992833f7984553de28d1b794f1a8be4e9df95e4`.
The current remote head, Q178, Q150--Q153, Q159, their primitive frozen
evidence and current CI were read first. Historical evidence, the physical
primary family, roots, profile, pressure gauge, supports, accepted scales and
the `C^-20` width remain unchanged.

Q178 isolated the raw primary covariance power `Q^(-1-2h*kappa)` and asked
for the complete mean cancellation. The later frozen source chain shows that
this is not the first current obstruction. The initial cancellation and the
first same-physical-family feedback are already available. This checkpoint
binds them to Q178 and records the first quantifier that is still open.

## 1. Complete initial cancellation

For the Q149 physical primary `w`, retain the full symmetric covariance

\[
 W_{ab}=\langle w_aw_b\rangle_\theta,
 \qquad a,b\in\{r,\theta,z\}.                         \tag{1}
\]

All six entries are present, including `W_rz`, `W_theta z`, `W_zz`, axial
velocity, radial velocity and the complete curl square. The pressure memory
is also complete:

\[
 q^{2A}p_B=\Pi_0+\chi_1q^{2h}
 \left(K_1-{\mu^2\over4X}\right)+\chi_2q^{4h}K_2.    \tag{2}
\]

Both exceptional tails cancel only by subtracting the same exact base
equation. The compact pressure and axial moment retain the entire moving
product `c_rho P`, including `(partial_Z c_rho)P`.

At the initial zero-old-mean state, the zero-Haar temporal inverse and
compact curl give exactly zero averaged fast velocity. After reconstructing
pressure, applying the slow five-map and reconstructing again,

\[
 P^{new},J_\theta^{new},J_z^{new}\in S_2,             \tag{3}
\]

and the angular-mean annular residual has order

\[
 2-3\kappa,qquad
 Q^{h(2-3\kappa)}.                                    \tag{4}
\]

The exponent in (4) is positive. Thus Q178's negative raw covariance power
is canceled by the exact six-entry identities at this initial stage. It is
not absorbed term by term. The limiting loss from order two is the retained
fast viscosity term.

## 2. Same-family principal stress and first feedback

The canonical and inherited physical primaries remain different objects.
Polynomial frame/root closeness cannot transfer a positive epsilon power
between them. The inherited physical family instead has its own exact root
and cover identities. For every active physical label `beta`,

\[
 H_\beta y_\beta=T_\beta,qquad
 a_\beta T_\beta=\Sigma_0,qquad
 \sum_\beta\eta_\beta^2=1.                            \tag{5}
\]

The band exponent in the second identity is
`-2A+h+A+1/2=h-A+1/2=0`. Consequently

\[
 \overline{(W_{0,principal})}_{re}
 =\sum_\beta\eta_\beta^2a_\beta H_\beta y_\beta
 =\Sigma_{0,e},\qquad e\in\{\theta,z\}.              \tag{6}
\]

No partition is divided out and neither sign is counted twice. The complete
curl remainder has order `2-2*kappa`; after the `epsilon^-1` stress scaling
it has order `1-2*kappa`. The required order is
`D-1=1/2-2*kappa`, leaving the exact margin `1/2`.

With

\[
 B=1-\kappa,\qquad C=B+\tfrac12,
 \qquad D=C-\kappa=\tfrac32-2\kappa,                 \tag{7}
\]

the exact scalar and row subtractions give postfast scalar order `C`, fast
row order `D`, slow increment order `C`, and both orientations of the first
fresh mean/current-wave feedback at order `D`. These identities retain all
old/new crosses, new squares, cylindrical terms, viscosity, compact pressure,
the nonzero `A1` remainder, axial velocity and (2).

## 3. Precise all-finite-cycle obstruction

The base and first feedback do not establish the target identities for every
later recomputed finite cycle. For each finite cycle `j` and each fixed mixed
jet `I`, the actual signed/full-covariance target must satisfy

\[
 |D^I\Xi_j|\le C_{j,I}\epsilon^{\alpha_j}\zeta
 y_a^{-p_{j,I}}y_b^{-q_{j,I}},                        \tag{8}
\]

together with both exact complete moments

\[
 \int R^2\overline{E_{\theta,j}}
 =\epsilon\,\partial_ZJ_{\theta,j},                  \tag{9}
\]
\[
 \int R\overline{E_{z,j}}
 =\epsilon\,\partial_Z(J_{z,j}+c_\rho P_{mean,j}).   \tag{10}
\]

Equations (8)--(10) must retain current axial velocity, all six covariance
entries, `Pi1/Pi2`, the whole moving pressure product, both velocity moments,
zero-Haar data and the reserved hull. The weight `zeta` vanishes at both flat
shell edges, while the reciprocal signed root costs half a weight. Finite
ordinary compact seminorms do not imply divisibility by `zeta`, and support
alone does not imply the zero totals in (9)--(10).

Q153's generated-source, tail and complete-potential diagonal result is a
conditional downstream implication. Its finite-stage recurrence assumes the
all-cycle target contract; it does not prove (8)--(10) for every `j`.
Therefore the complete covariance/mean hierarchy, complete residual, force
admissibility, exclusion and blow-up remain unproved.

The exact tests verify rational exponents, the same-family assembly, moment
fixtures, provider flags and immutable provenance. They do not replace the
continuum induction required by (8)--(10).

Next: derive (8)--(10) inductively from the actual physical cycle
expressions, retaining the complete pressure and support memory, before
reusing Q153's downstream source and residual diagonal.

Reproduce with `python3 test_audit.py -v` and
`python3 audit.py --verify-record evidence.json --output fresh.json` here.
