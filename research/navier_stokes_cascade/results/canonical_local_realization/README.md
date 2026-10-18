# One canonical local realization, with complete representatives

Parent: `20bb5a62a4770cc93fb7372705b1991d03905a4e`.
Status: `conditional-canonical-single-local-realization-and-all-jet-flat-residual`.
On 18 October 2026 the parent's new canonical finite-cycle CI step passed;
the complete hosted workflow was still running. Historical files and
frozen records are unchanged.

The [canonical finite-cycle bridge](../canonical_finite_cycle_join/README.md)
now supplies the inputs for ONE new canonical local field. Reapply the
reviewed summation mechanism to its OWN complete potentials and pressure
differences, not the older realized field. Conditional on the retained
exact-core/background, defining/inverse/gap and common-chart analytical
premises, the result is real, divergence-free, smooth for t<1 and q>0, with

\[
\forall m,N\in\mathbb N_0,\qquad
|\mathcal R(u_c^*,p_c^*)|_m=O(q^N),\quad q\downarrow0.             \tag{1}
\]

Here R=u_t-Delta u+(u.grad)u+grad p, viscosity is one, and the norm includes
all Cartesian space-time derivatives through order m. This is a conditional
analytic implication, not an unconditional complete-residual improvement,
construction of an admissible global force or blow-up proof. No actual
infinite coefficient or cutoff sequence is numerically evaluated.

The manuscript Section 9, (9.15)--(9.21), is being audited, not used as a
conclusion. Its PDF SHA256 is
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The generic mechanism in [realization foundation](../realization_foundation/README.md)
is reusable under its displayed hypotheses. Its historical fields are NOT
identified with the canonical family, even though their order budgets agree.

## 1. Supply the hypotheses from the canonical maps

Keep Md64, unit viscosity, original h,N, fixed positive lambda, radius,
width C_setup^-20, supports, carriers and the same accepted qstar<1.
Put A=1/2+h, D=1/2-h and

\[
q=1-t+z^2q^{2h},\quad \eta=z/q^D,\quad X=r^2/(2q),\quad
L=1-2h\eta^2\ge1-2h>0,\quad0<h<1/100.                            \tag{2}
\]

The parent constructs the canonical seed and all finite canonical states
U_c^[j], on this SAME domain, with B0=3/2-2kappa,
delta=1/2-kappa, rho=1-kappa, kappa=10^-5, Bj=B0+j delta, and

\[
|\mathcal R(U_c^{[j]})|_m\le C_{jm}q^{hB_j-K_m},
\qquad K_m=2A+3/2+2m.                                           \tag{3}
\]

This includes the full means, radial pressure reconstruction, exterior/base
errors, viscosity, interactions and additive Gaussian cutoff ledger, at
each FIXED finite stage. Stage constants need not be uniform.

For the cycle j-1 to j, collect the newly computed physical sourced and
signed potentials and the actual fast/slow azimuthal mean potentials:

\[
A_j=A_j^{wave}+\Psi_j e_\theta,\quad B_j^{sw}=b_je_\theta,
\quad\partial_\theta b_j=0,\quad p_j=p_c^{[j]}-p_c^{[j-1]},
\qquad \Delta u_j=\operatorname{curl}A_j+B_j^{sw}.                 \tag{4}
\]

The superscript sw distinguishes the vector swirl from the scalar order
Bj. This decomposition is supplied by the actual operators: it is not
inferred from a velocity bound or a claim that every compact solenoidal
velocity has a compact full vector potential. The mean torus map is
theta-independent. Its radial/axial component is the exact curl of Psi,
including actual compact-inverse A1 corrections; direct swirl stays separate.

With B=B_(j-1), C=B+1/2, wave potentials and pressures have order C, fast
means order C-kappa and slow means C. Current full pressure reconstruction
loses at most kappa. A sufficient common weakened ledger is:

| Actual canonical representative | Normalized order |
|---|---|
| Complete wave velocity | B-kappa |
| Wave potential | B+1/2-kappa |
| Mean azimuthal potential, direct swirl, full pressure difference | B+1/2-2kappa |
| Mean radial velocity | B+3/2-2kappa |

Since B_(j-1)=delta*j+rho, every row has slope delta and positive margin
over delta*j; the weakest is 1-2kappa. Flat pieces fit each fixed-stage,
fixed-jet order. Own canonical Hc/ac and current sources/means are used at
every map. No historical root, source, pressure target or potential value
is imported by polynomial closeness.

All representatives have smooth zero extensions supported in the same
closed annulus Xa<=X<=Xb, Xa>0, and agree under the prescribed physical
chart identifications. Both exceptional base pressure tails remain in p0;
only SAME-base changing parts cancel in (4). All six covariance entries,
axial velocity, moving pressure moment and intermediate updates telescope.
These support/gauge properties follow conditionally from the canonical
source grammar and compact mean maps, not from the algebra fixtures.

Physical normalizations Q^-A,Q^-2A,Q^(1/2-A) are bounded by q^-2A on Q~q.
Each fixed derivative costs at most q^-2, with fixed logarithmic powers.
Finite-stage harmonic ranges affect constants, not these exponents. Thus,
for Zj=(Aj,Bj^sw,pj), gj=h delta j,

\[
|Z_j|_m\le C_{jm}(1+|\log q|)^{P_{jm}}q^{g_j-\ell_m},\quad
\ell_m=2A+2m,\quad L_m=\ell_{m+1},\quad H_m=L_m+1.                \tag{5}
\]

Lm reserves one curl and its cutoff derivative. The fixed base and each
finite total state satisfy |Uc^[J]|m<=C_Jm q^-Hm after absorbing fixed logs.
This does not assert a stage-independent constant.

## 2. Cut complete potentials, swirl AND pressure

Fix smooth chi=1 for s<=1/2, zero for s>=1. For chi_a=chi(aq), use

\[
V_{j,a}=\operatorname{curl}(\chi_a A_j)+\chi_a b_je_\theta
 =\chi_a\Delta u_j+\nabla\chi_a\times A_j,\qquad
P_{j,a}=\chi_a p_j.                                             \tag{6}
\]

Each V is exactly divergence-free. Direct swirl is divergence-safe because
q and b are theta-independent. The mean radial commutator is
-a chi'(aq) q_z Psi_j, and pressure differentiation keeps p_j grad chi_a.
All wave curl terms remain. Cutting only the velocity need not preserve
divergence. Multiplication and curl do not enlarge closed supports. For
each t<1 the annulus is away from the axis, and zero extension is smooth.

The exact weighted operators on q^w F(eta,sigma), sigma=aq, are

\[
\partial_z:q^{w-D}\left[\frac{2w\eta}{L}F+
 \frac{1-\eta^2}{L}F_\eta+\frac{2\eta\sigma}{L}F_\sigma\right],
\quad
\partial_t:q^{w-1}\left[-\frac wL F+
 \frac{D\eta}{L}F_\eta-\frac\sigma L F_\sigma\right].              \tag{7}
\]

Starting with F=chi(sigma), each positive jet contains sigma^k chi^(k),
k>=1, and bounded rational eta factors with denominators powers of L.
Positive derivatives are supported on 1/2<=sigma<=1. Induction gives
|partial_z^a partial_t^b chi_a|<=C_ab q^(-Da-b), uniformly in the scale.
Leibniz and (5), with the extra curl, imply

\[
|(V_{j,a},P_{j,a})|_m+|\chi_a Z_j|_m
 \le\widetilde C_{jm}(1+|\log q|)^{\widetilde P_{jm}}
 q^{g_j-L_m}.                                                   \tag{8}
\]

### Gauge audit: equal uncut states do not mean equal cut states

There is an additional representation constraint. Although curl grad phi=0,

\[
\operatorname{curl}(\chi(A+\nabla\phi))-
 \operatorname{curl}(\chi A)=\nabla\chi\times\nabla\phi,\qquad
\nabla[\chi(p+c(t))]-\nabla(\chi p)=c(t)\nabla\chi.               \tag{9}
\]

Both can be nonzero. Divergence-free alone does not bound the first
artifact. An arbitrary stage-dependent gauge can violate (5) or support,
despite identical uncut velocities/pressure gradients. Use the actual
operator representatives and compact reconstruction gauge BEFORE choosing
scales. No unsupported change of potential gauge is allowed. Admissible
changes satisfying the same increasing-order ledger can instead be
included in the finite constraint lists and handled by the same argument.
A COMMON c(t) added to every full state cancels from every pressure
increment; it only changes p0 by a spatial constant. Independent c_j(t)
do not cancel through different chi(a_jq). Tests distinguish these cases.

There is a continuum counterexample to using finite residual bounds alone.
Take uncut u^[j]=0, p^[j]=j: EVERY finite residual is exactly zero, but
the pressure increments are 1, not increasing-order compact representatives.
Cut them as in (12). Choose sigma in (1/2,1) with chi'(sigma)!=0,
eta=1/2 and q_j=sigma/a_j. Strict doubling makes all OTHER positive cutoff
derivatives zero there. Since z=eta*q_j^D and 1-t=(1-eta^2)*q_j>0,
these are actual points in the implicit-q chart, and

\[
|\mathcal R(0,\sum_i\chi(a_iq))|\ \ge
\frac{\sigma|\chi'(\sigma)|}{1-h/2}\ q_j^{-D}\longrightarrow\infty.
                                                                    \tag{9a}
\]

This is a smooth locally finite pressure on q>0, not the canonical
pressure. It violates the representative hypotheses, not the theorem.
It proves why equal uncut equations plus good finite residuals cannot
replace the fixed-gauge, increasing-order representative ledger.

## 3. An exact sufficient scale recipe on the WHOLE interval

For every j collect (8) through m=j. Pad C>=1,P>=0. Write gamma=gj/2>0
and c>=ln C. For s=-ln q>=1, ln(1+s)<=sqrt(s), so

\[
c+P\ln(1+s)-\gamma s+j\ln2
\le c+P^2/(2\gamma)-\gamma s/2+j.                               \tag{10}
\]

Consequently any sj=ln a_j strictly larger than

\[
\max\{s_{j-1}+1,\ \ln(qstar^{-1})+1,\ 1,
       \max_{(c,P)}[2(c+j)/\gamma+P^2/\gamma^2]\}                 \tag{11}
\]

satisfies a_j>2a_(j-1), a_j>2/qstar, and
C(1+|ln q|)^P q^(gj/2)<=2^-j for EVERY 0<q<=1/a_j. The same inequality
holds for all larger s because the right side of (10) decreases.
The implementation uses rational upper bounds and adds one to (11),
avoiding roundoff in its sufficient majorant. Synthetic extreme constants
test arithmetic only; actual finite-stage constants need not be supplied
numerically for this existence implication. Rational upper bounds always
exist for the finitely many real constants at a fixed stage.

Choose this ONE sequence and define

\[
u_c^*=u_c^{[0]}+\sum_{j\ge1}V_{j,a_j},\qquad
p_c^*=p_c^{[0]}+\sum_{j\ge1}P_{j,a_j}.                           \tag{12}
\]

On compact subsets of q>0, only finitely many stages and physical labels
are active. Thus this is one smooth real divergence-free field with the
same protected inner velocity and exterior. No raw series convergence or
uniform stage constants is required. For FIXED J>=max(1,m), q<=1/(2a_J),
the first J cutoffs equal one with zero positive derivatives, and telescope.
For j>J use (8)--(11), or exact zero if inactive, to obtain

\[
|U_c^*-U_c^{[J]}|_m
 \le2^{-J}q^{h\delta(J+1)/2-L_m}.                               \tag{13}
\]

This holds across transition intervals. The cutoff sequence is chosen
once for all fixed jets and targets; the comparison stage is NOT j(q).

## 4. Complete residual and remaining scope

For (e,pi)=Uc*-Uc^[J], the EXACT residual difference is

\[
e_t-\Delta e+\nabla\pi+(u_c^{[J]}\cdot\nabla)e+
 (e\cdot\nabla)u_c^{[J]}+(e\cdot\nabla)e.                        \tag{14}
\]

With E=|(e,pi)|_(m+2), its m-jet is bounded by
C_m(E+C_Jm q^-H_(m+2) E+E^2). Fix m,N FIRST; choose a finite J>=m+2 so

\[
h\delta(J+1)/2\ge N+1+L_{m+2}+H_{m+2},\qquad
hB_J-K_m\ge N+1.                                               \tag{15}
\]

Both slopes are positive. Equations (3),(13)--(15) prove (1) for the
SAME canonical field. Time, viscosity, pressure, both nonlinear crosses
and new squares are all charged. The complete additive Gaussian ledger
belongs to this ONE fixed finite comparison state; it is not summed
separately over infinitely many stages or fed into a forward inverse.

At diagnostic h=1/200,m=N=0 the sufficient budget is J=12816. This is
neither evaluation of that many fields nor an optimized count, progress
percentage or time-to-proof estimate. The parent's arbitrary-diagonal
counterexample remains valid; (11)--(15) supply the missing controlled
construction rather than contradicting it.

Exact controls reuse the independently checked implicit-q operators through
order six and the complete Cartesian residual/curl polynomial fixtures.
New local polynomial gauge fixtures detect both artifacts (9), common vs
independent pressure gauges and premature pressure truncation. They are
not actual supported canonical pulses or a continuum proof. All-order
claims follow from the analytical estimates and their retained hypotheses.

This closes only the canonical finite-state-to-one-local-field implication.
Next audit the closed-time auxiliary atlas and ALL complete one-sided
potential, velocity and pressure jets away from the singular corner for
THIS field. Then join localization/initialization, admissible smooth
compact forcing, whole-space energy and global comparison. Older downstream
fields are not automatically identical. Retained exact-core/background
and analytical defining/existence/inverse/gap premises still need independent
discharge before unconditional pulse, residual or blow-up conclusions.

```sh
python3 research/navier_stokes_cascade/results/canonical_local_realization/test_audit.py -v
python3 research/navier_stokes_cascade/results/canonical_local_realization/audit.py \
 --verify-record research/navier_stokes_cascade/results/canonical_local_realization/evidence.json \
 --output /tmp/canonical_local_realization_fresh.json
```
