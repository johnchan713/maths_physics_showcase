# One current local field from the finite-stage foundations

Parent: `07b30739332662b1bca55018fa735bb0391ae5d5`.
Status: `current-source-single-local-realization-relative-to-retained-core-and-finite-analytic-lemmas`.
Parent hosted CI was still running, without a failed step, when inspected
on 2026-10-17. Historical sources and frozen evidence are unchanged.

The [current finite induction](../finite_induction_foundation/README.md)
supplies one complete state at each finite stage. This checkpoint joins its
ACTUAL complete representatives to ONE real divergence-free field on the
inherited local domain, smooth for t<1 and q>0, with

\[
\forall m,N\in\mathbb N_0,\qquad
|\mathcal R(u_*,p_*)|_m=O(q^N),\qquad q\downarrow0. \tag{1}
\]

Here R=u_t-Delta u+(u.grad)u+grad p, viscosity is one, and |.|_m is the
Cartesian space-time jet norm through total order m. This is relative to
the retained exact-core and finite defining/existence/inverse/gap lemmas
in the reviewed foundation chain. It does NOT discharge those dependencies,
join the closed-time endpoint atlas, supply an admissible global force,
or prove blow-up. Finite algebra diagnostics below are distinct from the
analytic argument. No actual infinite coefficient or cutoff sequence is
numerically evaluated.

The manuscript source is Section 9, especially (9.15)--(9.21), Lemma 9.8
and Proposition 9.9, in the supplied PDF with SHA256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its final conclusion is never an input.

The earlier [conditional realization](../current_potential_realization/README.md)
used B0=3/2-6kappa and delta=1/2-3kappa. Its summation mechanism is sound
under its stated inputs. We now independently match that mechanism to the
reviewed current source grammar, complete physical pulses and compact mean
maps, with B0=3/2-2kappa and delta=1/2-kappa. This closes the single-field
quantifier gap left by the preceding current finite foundation; it does
not merely rename its per-target finite states.

## 1. Complete representatives and stage-independent physical losses

Keep Md64, h,N, positive lambda, width C_setup^-20, original roots, supports,
dyadic labels, chart map and operation order unchanged. Write

\[
A=\tfrac12+h,\quad D=\tfrac12-h,\quad
q=1-t+z^2q^{2h},\quad \eta=z/q^D,\quad X=r^2/(2q),
\quad L=1-2h\eta^2\ge1-2h>0. \tag{2}
\]

Use the same qstar<1, 4Q<qstar and enlarged accepted charts. For j>=1,
the COMPLETE j-1 to j cycle has the actual fixed representatives

\[
\Delta u_j=\operatorname{curl}A_j+b_je_\theta,
\quad A_j=A_j^{wave}+\Psi_je_\theta,\quad\partial_\theta b_j=0,
\quad p_j=p^{[j]}-p^{[j-1]}. \tag{3}
\]

The sourced and signed full potentials/pressures have order C=B+1/2,
where B=B_(j-1). Actual fast means have H=C-kappa, slow means C; their
radial velocities have one extra order. The compact mean foundation
supplies the radial/axial azimuthal potential separately from direct swirl.
Thus a sufficient common (weakened) ledger is:

| Actual increment | Sufficient normalized order |
|---|---|
| Complete wave velocity | B-kappa |
| Wave potential | B+1/2-kappa |
| Mean azimuthal potential, direct swirl, COMPLETE pressure difference | B+1/2-2kappa |
| Mean radial velocity | B+3/2-2kappa |

These are bounds on representatives from the actual finite maps, not a
deduction of compact potential existence from a velocity norm. In
particular no compact FULL vector potential is claimed for direct swirl.
The source grammar, finite harmonic induction, full mean recomputation
and same-band smooth zero extensions propagate these bounds at each finite
stage. Stage constants and harmonic counts may grow arbitrarily.

The exact current indexing identity is

\[
B_{j-1}=\delta j+\rho,\qquad
\delta=\tfrac12-\kappa,\quad\rho=1-\kappa,
\quad B_0-\delta=\rho. \tag{4}
\]

Every row exceeds delta*j; the weakest margin is 1-2kappa>0. Physical
velocity, pressure and potential factors are Q^-A, Q^-2A and Q^(1/2-A).
Every fixed Cartesian derivative costs at most q^-2; the reviewed graph,
frame and phase costs are strictly smaller than two. Logs and finite
harmonics affect the constants/polynomial degrees. Q~q on the SAME charts.
Consequently Z_j=(A_j,b_je_theta,p_j) obeys

\[
|Z_j|_m\le C_{jm}(1+|\log q|)^{P_{jm}}q^{g_j-\ell_m},
\quad g_j=h\delta j,\quad\ell_m=2A+2m. \tag{5}
\]

This includes actual compact flat velocities/potentials. The designated
Gaussian residual tails are NOT new inverse inputs. Complete pressure is
in the fixed compact reconstruction gauge and telescopes after ALL four
operations. Both Pi1/Pi2 background tails, axial pressure and the moving
moment remain in each full state; subtraction cancels only the SAME fixed
base. Adding a time-only pressure gauge BEFORE a cutoff would generate
c(t) grad chi and is not an allowed silent replacement of (3).

## 2. Cut the potentials and keep all pressure and curl terms

Fix smooth chi=1 on s<=1/2 and chi=0 on s>=1. Set chi_a=chi(aq).
The cut increment is exactly

\[
V_{j,a}=\operatorname{curl}(\chi_a A_j)+\chi_a b_je_\theta
 =\chi_a\Delta u_j+\nabla\chi_a\times A_j,
\qquad P_{j,a}=\chi_a p_j. \tag{6}
\]

Every V is divergence-free: the first term in its defining expression is
a curl; the second is axisymmetric swirl. The mean radial commutator is
-a chi'(aq) q_z Psi_j. Pressure differentiation retains p_j grad chi_a.
The supported representatives vanish smoothly outside the same closed
annulus Xa<=X<=Xb, Xa>0. Products/curls do not enlarge support. At fixed
t<1 the annulus stays away from the Cartesian axis, so its zero extension
is smooth there. The protected inner field and inherited exterior agree
exactly with the finite states. Mean supports may fill the auxiliary torus;
no wave-label separation is falsely imposed on them.

For an explicit all-order cutoff argument, put sigma=aq. For any coefficient
F(eta,sigma), weighted differentiation of q^w F gives

\[
\partial_z:\quad q^{w-D}\left[
\frac{2w\eta}{L}F+\frac{1-\eta^2}{L}F_\eta
 +\frac{2\eta\sigma}{L}F_\sigma\right],
\quad
\partial_t:\quad q^{w-1}\left[
-\frac wL F+\frac{D\eta}{L}F_\eta-\frac\sigma L F_\sigma\right]. \tag{7}
\]

Start with F=chi(sigma), w=0. Every positive derivative is a finite sum
of sigma^k chi^(k)(sigma), k>=1, with rational eta coefficients whose
only denominators are finite powers of L. On |eta|<=1 these are bounded;
positive chi derivatives are supported on 1/2<=sigma<=1. Induction proves
|partial_z^alpha partial_t^beta chi_a|<=C_alpha,beta q^(-D*alpha-beta),
uniformly in the cutoff scale a (the derivative indices are fixed).
Radial/angular derivatives
vanish. This is one fixed domain for every fixed jet, not a sequence of
jet-dependent starting bands. Leibniz and the extra curl give

\[
|(V_{j,a},P_{j,a})|_m+|\chi_a Z_j|_m
 \le\widetilde C_{jm}(1+|\log q|)^{\widetilde P_{jm}}
 q^{g_j-L_m},\quad L_m=\ell_{m+1}, \tag{8}
\]

uniformly in a>=1. The fixed base and finite sums obey
|U^[J]|_m<=C_Jm q^-H_m with H_m=L_m+1, absorbing fixed logarithms by
q^-1. No uniformity in J or m is asserted.

## 3. ONE diagonal sequence, then a fixed comparison stage

At step j collect (8) for 0<=m<=j; pad constants C>=1,P>=0. Choose
a_j>=2a_(j-1), a_1>2/qstar, so that for EVERY requirement

\[
C(1+|\log q|)^P q^{g_j/2}\le2^{-j}
\quad\hbox{for ALL }0<q\le1/a_j. \tag{9}
\]

This whole-interval condition exists for any finite constants. Put
s=-log q, gamma=g_j/2. For s>=1, log(1+s)<=sqrt(s), so

\[
\log C+P\log(1+s)-\gamma s
 \le\log C+P^2/(2\gamma)-\gamma s/2. \tag{10}
\]

It suffices to choose log a_j larger than log a_(j-1)+log2,
log(2/qstar), 1, and every
2(log C+j log2)/gamma+P^2/gamma^2. The actual finite constants need not
be evaluated or satisfy any growth cap. Endpoint sampling alone is not
the argument. Define

\[
u_*=u^{[0]}+\sum_{j\ge1}V_{j,a_j},\qquad
p_*=p^{[0]}+\sum_{j\ge1}P_{j,a_j}. \tag{11}
\]

On compact subsets of q>0, only finitely many j and physical dyadic labels
are active. Thus (11) is a real smooth field, exactly divergence-free.
For FIXED J>=max(1,m) and q<=1/(2a_J), the first J cutoffs and derivatives
are exactly their constant-one values, so full state/pressure differences
telescope. Equation (9) sums the later cut contributions, including their
transition commutators, to

\[
|U_*-U^{[J]}|_m\le2^{-J}q^{h\delta(J+1)/2-L_m}. \tag{12}
\]

This proves an asymptotic expansion by finite STATES, not convergence of
the uncut series: coefficients exp(j^2) can defeat that series at every
fixed q. Neither fixed-q contraction nor a substitution J=J(q) is used.

For e=u_*-u^[J], pi=p_*-p^[J], the COMPLETE difference is

\[
\mathcal R(U_*)-\mathcal R(U^{[J]})
=e_t-\Delta e+\nabla\pi+(u^{[J]}\cdot\nabla)e
 +(e\cdot\nabla)u^{[J]}+(e\cdot\nabla)e. \tag{13}
\]

With E=|(e,pi)|_(m+2), its m-jet is bounded by
C_m(E+C_Jm q^-H_(m+2) E+E^2). Fix m,N FIRST and choose a finite J>=m+2
with

\[
\frac{h\delta(J+1)}2\ge N+1+L_{m+2}+H_{m+2},\qquad
hB_J-K_m\ge N+1,\quad K_m=2A+\tfrac32+2m. \tag{14}
\]

Both slopes are positive. The current finite foundation bounds the full
finite residual by C_Jm q^(hB_J-K_m), including viscosity, full means,
pressure and exact additive tails. Equations (12)--(14) prove (1) for this
SAME field for EVERY fixed m,N. The Gaussian ledger is used at this ONE
fixed finite stage; it is never separately summed over infinitely many
stages or forwarded through an inverse. All cutoff time, pressure, viscous
and old/new/new-square terms are included via (13).

At DIAGNOSTIC h=1/200,m=N=0 a sufficient coarse budget is J=12816. This
is an unevaluated existence budget, not a numerical construction, progress
percentage, optimized stage count or estimate of distance to a proof.

## 4. Evidence, limits and the next obligation

`realized_checks.py` verifies 226 exact rational mixed-cutoff coefficient
comparisons through total order six at two allowed h values, using the
ACTUAL implicit-q weighted operators (7). The all-order result comes from
the induction and denominator gap, not that finite sample. Separate exact
Cartesian polynomial fixtures check full curls, divergence, pressure
products and all six groups of (13), with two different cutoffs. Every
individual omitted group, the naive velocity cutoff, pressure-gauge
cutoff and viscous commutators produce nonzero defects. Those polynomial
fixtures are LOCAL algebra controls, neither the actual implicit q nor
C-infinity compact pulses; actual-q checks are separately identified.
Exact rational budgets check current indexing and (14); extreme SYNTHETIC
constants check the recipe for (9), not actual infinite scales.

The parent frozen chain is reproduced, preserving 961 inherited provenance
bindings and 18 new direct source bindings. Use fresh output paths:

```sh
python3 research/navier_stokes_cascade/results/realization_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/realization_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/realization_foundation/evidence.json \
  --output /tmp/realization_foundation_fresh.json
```

Next join the closed-time auxiliary atlas and complete one-sided potential,
velocity and pressure jets away from the singular corner to THIS field.
Local finiteness for q>0 and flatness as q->0 do not by themselves supply
that endpoint atlas. Then join complete localization/initialization,
admissible smooth global forcing, whole-space energy and the global solution
comparison. Older conditional downstream results are not automatically
promoted. Retained exact-core and finite analytic lemmas must still be
discharged before unconditional residual or blow-up claims.
