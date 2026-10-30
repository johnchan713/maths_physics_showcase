# Q17: projected mean supports and full-pressure aggregation

Parent: `b5b4383ec413a9d7d9272c7242ce0d6ead950887` (Q16).
This closes Q16's separate mean/pressure **overlap-count** obligation on
its retained geometry. It does not certify its normalized coefficient inputs.
The new join uses the frozen slow mesh and compact maps; it does not choose
new bands, profiles, roots, carriers, cutoffs or pressure gauges. Historical
evidence is unchanged. All statements are on q>0 and the accepted later tail.

## 1. Count after radial and auxiliary support have filled

The frozen `auxiliary_geometry_audit`, section 6, proves at most three
active closed dyadic bands, three mesh supports per coordinate and two
signs. The common covering difference is at most three, with at most
14^3=2744 lifts. It also proves that a fixed radial line in one native
chart meets at most

    18 (ceil(L_R S_i^3)+3)

slow labels, where L_R is an upper bound for the length of the entire
fixed radial chart hull. Use a maximum native hull length over the fixed
chart family, including the reserved profiles; L_R is independent of band.
The coarsest common chart of `core_primary_aggregate_atlas` has
S_i <= (25/16) S, S=ell_*^2 >=1. Thus, even counting every common lift,
the number of original labels whose slow support projects to a fixed
(Z,T) is bounded by

\[
 N(S)\le54\cdot2744\{\lceil L_R(25S/16)^3\rceil+3\}
 \le C_{\rm proj}S^3,\quad
 C_{\rm proj}=54\cdot2744\{L_R(25/16)^3+4\}.       \tag{1}
\]

The second inequality uses ceil(x)<=x+1 and S>=1. This intentionally
coarse count is for projected supports, not the sharper count at one
radial/auxiliary point. The latter's single-wave separation cannot be
assigned to means. Derivatives of smooth zero extensions do not enlarge
closed support, so the same count works at every fixed jet order.

Why does the projection suffice? The shifted T_e and A_e integrate along
radius and translate the torus while keeping physical z,t fixed. Haar
averaging and the zero-Haar temporal inverse act only on the torus. A
radial moment multiplied by a fixed q-based reserved profile keeps z,t
fixed too. These operations preserve the projection onto z,t although
they can fill radial gaps and the entire torus. The scalar q(z,t) is
independent of radius, hence no additional dyadic band or older covering
index is introduced. Physical slow derivatives preserve that closed
projection. The original global physical definitions supply chart gluing;
independent fixed-R bumps in different charts would not do so.

Equation (1) applies to an original label sum after any finite composition
of these linear maps. It is NOT a count of nonlinear monomials at all
stages. For nonlinear products use the aggregate norm algebra below;
there is no claim that only N(S) product words exist.

## 2. Aggregation costs logs, with no power of q

Suppose fixed-order native estimates, converted to one common physical
normalization and made unweighted by Q15, give

\[
 \|D^I F_a\|\le K_I q^E P_I(S),\quad |I|\le m,          \tag{2}
\]

uniformly over the labels a. Here P_I has nonnegative coefficients and
finite degree. Every normalization factor, covering derivative, fixed
width cost and finite operator constant is inside K_I P_I. Then (1) and
the triangle inequality give

\[
 \left\|D^I\sum_a F_a\right\|
 \le K_I C_{\rm proj}q^E S^3P_I(S).                    \tag{3}
\]

One may instead bound the aggregate before a linear map, using its
whole-hull operator estimate. This often removes the counting loss, as
in the primary atlas proof. Equation (3) remains valid without such a
sharpening. Normalized Haar averaging has norm one; a covering does not
multiply its measure by 2744. The deck factor in (1) is only a deliberately
redundant count of local terms, never an averaging normalization.

For multiindices, the exact Leibniz bound is

\[
 \|D^I(FG)\|\le\sum_{J\le I}{I\choose J}
       \|D^JF\|\|D^{I-J}G\|.                         \tag{4}
\]

Thus sums, products, all covariance entries, radial moments, and the
retained bounded linear mean maps preserve a finite polynomial S bound
at each fixed node of the finite construction. Derivative losses request
higher input jets; they do not justify those input jets. Inverse root and
frame operations still require their separate quantitative lower bounds.
A finite graph is not by itself a proof of those bounds or of its assigned
epsilon orders. No constant uniform in construction stage is asserted.

The Q16 comparison S<=16|log q|^2 turns S^3 into at most
4096(1+|log q|)^6. Hence each use of (3) adds at most six to the log degree
and spends no q reserve. Products add log degrees. The smallest Q16 margin
remains h(3/2-kappa)+m(1-5h/2)>0 at the exact root h; Q12/Q13's symbolic
cutoff selection can absorb any resulting fixed finite log degree using
its existing strict reserve. This does not evaluate the actual constants.

## 3. Full pressure is an aggregate compact increment

At every current mean state use the same linear map

\[
 P=\int\overline{g_r}\,dR,\qquad
 p_m=T_0(g_r-\rho P).                                  \tag{5}
\]

Its input has zero *radial Haar moment*, not necessarily zero shifted
moment pointwise. The exact identities are

\[
 D_rp_m-g_r=-\rho P-A_0(g_r-\rho P),\qquad
 \overline{A_0(g_r-\rho P)}=0.                          \tag{6}
\]

The A0 term remains in the pointwise residual; no support count deletes
it. The compact map and moving rho preserve the same slow projection, so
(3) applies to its increments. Because rho is fixed across states,
Delta p_m=T0(Delta g_r-rho int bar(Delta g_r)). The *nonlinear* change
Delta g_r must first include all old-new and new-old products and all
six complete-wave covariance entries, including axial terms. A pressure
map's linearity does not make g_r linear in velocity.

Radial integration by parts on the compact Haar pressure gives

\[
 \int R\bar p_m=-\tfrac12\int R^2\bar g_r+c_\rho P,
 \quad c_\rho=\tfrac12\int R^2\rho.                   \tag{7}
\]

Consequently the axial moment uses the entire moving derivative
\(\epsilon\partial_Z(J_z+c_\rho P)\), including
\(\epsilon c_{\rho,Z}P\). For the fixed Q4 profile,
\(c_\rho=s\int x^2\widehat\rho/2\), so this term need not vanish.
The auxiliary-independent A1 zero-moment case is exactly compact, while
the general axial velocity remains gamma_d-A1 gamma_d. Both enter the
next recomputation.

Total pressure is the same base plus all compact/full wave changes.
The base Pi1/Pi2 memory is not compactified by this argument. It cancels
only when subtracting two states with the identical complete base and
gauge; its axial gradient is retained in each state. Q16 concerns those
complete fixed-gauge differences. Equation (7) is not permission to
replace them with a leading-pressure approximation.

## 4. Scope, checks and next obligation

The support aggregation needed by Q16 is now a consequence of the frozen
geometry and the exact mean/pressure maps, conditional on their retained
smooth weighted inputs. It is no longer an independent unknown label
multiplicity. Q16's normalized coefficient inequalities remain open here:
for each representative component c and each finite corrected Q14 request,
prove for the *actual* freshly recomputed field

\[
 \sup_{\ell\ge\ell_*}(1+|\log\epsilon_\ell|)^{-P_{j,R}}
 \epsilon_\ell^{-\alpha_c}
 \max_{|I|\le R}\sup_{\Omega_{j,\ell}}|D^I F_{c,j,\ell}|<\infty,
                                                               \tag{8}
\]

with the appropriate weighted edge/pole bounds before Q15, original
whole-path source bounds, root gaps and full pressure moments. In (8),
F denotes the normalized coefficient envelope before physical fast-phase
differentiation. Q5's source induction statement must be connected to
these component orders and Q14's fully expanded domains, rather than
promoted from a norm-container name or this count. That is next.

Exact rational tests check integer mesh counts against (1), aggregation
and Leibniz bounds, gap filling, the compact pressure identity, and the
moving-pressure commutator. Polynomial fixtures have finite endpoint
regularity and are not actual pulses. Parent evidence is regenerated and
its entire inherited hash graph checked. Neither these tests nor this
support theorem prove complete representative bounds, endpoint closure,
unconditional residual improvement, forcing regularity or blow-up.

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q17_mean_support_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q17_mean_support_join/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q17_mean_support_join/evidence.json --output FRESH.json
```
