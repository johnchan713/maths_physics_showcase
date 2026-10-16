# Current-state complete-potential diagonal realization

Status: `conditional-current-state-complete-potential-diagonal-realization`.

This closes a bounded stream: the [current finite-stage hierarchy](../finite_stage_induction/README.md)
now implies ONE real divergence-free local field with a complete momentum
residual flat in every fixed Cartesian space--time derivative. It does so
under the explicitly retained leading-profile, summed-background, inverse,
support and auxiliary-chart analytical premises. Neither tests nor this
conditional implication independently verify those premises. There is no
unconditional solution, admissible global force or blow-up claim.

The older [realization argument](../full_potential_realization/README.md)
used a different hypothetical hierarchy. Its frozen evidence is preserved
and reproduced, not silently treated as a realization of the new states.
Here the actual current cycle, complete representatives, stronger stage
gain and Gaussian-tail exclusion are joined explicitly. The actual infinite
stage coefficients and cutoff constants are not numerically evaluated.

## 1. Current input and complete representatives

Use the same inherited domain `t<1`, `0<q<qstar<1`, with

\[
q=1-t+z^2q^{2h},\quad 0<h<1/100,\quad A=1/2+h,\quad D=1/2-h,
\quad X=r^2/(2q),\quad\eta=z/q^D.
\]

The implicit derivative denominator is `1-2h eta^2>=1-2h>0`.
The latest actual audited finite state is `U^[0]`. The conditional
four-operation induction supplies `U^[j]`, for every finite integer j,
on the SAME domain, with

\[
\kappa=10^{-5},\quad\delta=1/2-3\kappa,\quad
B_j=3/2-6\kappa+j\delta,\quad\rho=1-3\kappa,
\]
\[
|\mathcal R(U^{[j]})|_m\le C_{j,m}q^{hB_j-K_m},
\qquad K_m=2A+3/2+2m.\tag{1}
\]

This is the COMPLETE residual, not merely its regular nonzero angular
part. At every fixed j,m its exact designated Gaussian cutoff tails and
other flat remainders have been absorbed by choosing their finite powers.
They are never forward-inverted. Constants can grow arbitrarily in j,m;
harmonics are finite at each stage but not uniformly bounded in j.

For the cycle `j-1 -> j`, collect the actual particular and signed wave
potentials, plus the actual fast and slow azimuthal mean potentials:

\[
\mathcal A_j=\mathcal A_j^{\rm wave}+\Psi_j e_\theta,
\quad\mathcal B_j=b_j e_\theta,\quad\partial_\theta b_j=0,
\quad p_j=p^{[j]}-p^{[j-1]},
\]
\[
u^{[j]}-u^{[j-1]}=\operatorname{curl}\mathcal A_j+\mathcal B_j.\tag{2}
\]

The mean torus evaluation depends on r,t, not theta, so the direct swirl
is axisymmetric even when it has auxiliary oscillations. Its divergence
is zero. No compact full vector potential for this direct swirl is assumed
or asserted; (2), with the separately divergence-safe swirl, is the
complete representation used here. Radial/axial mean velocity is its EXACT
azimuthal-potential curl, including the flat compact-inverse correction.
The full pressure retains the SAME p_B, including both exceptional base
tails, in p^[0]. Only changing parts telescope in p_j: these include wave
pressure, full six-covariance mean pressure, the complete axial-memory
contribution and its current moment bump, after ALL four operations.
A fixed pressure gauge is used throughout; no intermediate pressure
change or unchanged base tail is discarded.

Every representative has a smooth zero extension on the same closed
annulus `Xa<=X<=Xb`, `Xa>0`, with chart agreement. Distinct wave labels have
the inherited closed-support separation; means may fill their torus.
These are inputs established conditionally by the finite induction and
its operator lemmas, not inferred from manufactured examples.

Put `B=B_(j-1)=delta*j+rho`, `C=B+1/2`. The current sourced and signed
complete wave increments have order B and potentials/pressures C. Fast
means and their pressure/potential changes have order `H=C-kappa`; slow
ones have order C; their radial velocities gain one. Recomputing mean
pressure after each wave change costs at most kappa, as in the induction.
The following deliberately weaker common ledger reserves extra kappa:

| Complete entry | Sufficient normalized order | Margin over delta*j |
|---|---|---|
| Wave velocity | B-kappa | rho-kappa |
| Wave potential | B+1/2-kappa | rho+1/2-kappa |
| Mean azimuthal potential, direct swirl, complete pressure difference | B+1/2-2kappa | rho+1/2-2kappa |
| Mean radial velocity | B+3/2-2kappa | rho+3/2-2kappa |

Every row has EXACT slope delta, with smallest margin `1-4kappa>0`.
Thus this is an all-stage argument, not a large finite scan. Actual flat
pieces remain in the representatives and products; for each fixed stage
and jet order they fit this ledger by their all-orders flat estimates.

The inherited physical normalizations are `Q^-A`, `Q^-2A`,
`Q^(1/2-A)` for velocity, pressure, potential. On accepted charts `Q~q`.
The physical phase/graph jet estimates, flat edge envelopes and local
label separation bound each fixed derivative by a cost smaller than
`q^-2`; finite harmonic ranges change constants only. For
`Z_j=(A_j,B_j,p_j)`, the ledger therefore implies

\[
|Z_j|_m\le C_{j,m}(1+|\log q|)^{P_{j,m}}
q^{g_j-\ell_m},\quad g_j=h\delta j,\quad\ell_m=2A+2m.\tag{3}
\]

One extra derivative for the curl gives `L_m=ell_(m+1)=2A+2m+2`.
The fixed base and each finite total state obey
`|U^[j]|_m<=C_jm q^-H_m`, `H_m=L_m+1`, absorbing fixed logarithmic powers
by `q^-1`. None of these constants is assumed uniform in stage.

## 2. Whole-interval diagonal cutoffs, not a raw series

Fix smooth chi equal to one on `s<=1/2` and zero on `s>=1`.
For `chi_a=chi(aq)`, implicit differentiation and the chain rule give

\[
|\partial_z^a\partial_t^b\chi(\lambda q)|
\le C_{a,b}q^{-Da-b},\tag{4}
\]

uniformly in the cutoff scale lambda: derivatives are supported where
`1/2<=lambda*q<=1`, and cancel their scale factors there. Radial/angular
derivatives vanish. The exact cutoff increment is

\[
V_j=\operatorname{curl}(\chi(a_jq)\mathcal A_j)
       +\chi(a_jq)\mathcal B_j,
\quad P_j=\chi(a_jq)p_j.\tag{5}
\]

In particular `V_j=chi*Delta u_j+grad chi cross A_j`. For the mean potential
the extra radial velocity is `-a_j chi'(a_jq) q_z Psi_j`. Dropping it breaks
divergence. Direct-swirl cutting is divergence-safe since both b_j and q
are theta independent. Pressure is cut too: `p_j grad chi` is retained.
Using (3) at m+1 and (4), a common bound for the cut representatives,
velocity and pressure is `C_jm(1+|log q|)^P_jm q^(g_j-L_m)`, uniformly
in a_j. Include every such seminorm through m=j in a finite constraint list.

Choose ONE sequence `a_j>=2a_(j-1)`, `a_1>2/qstar`, so

\[
C(1+|\log q|)^Pq^{g_j/2}\le2^{-j}
\quad\hbox{for ALL }0<q\le1/a_j\tag{6}
\]

for its entire finite list. Such scales exist with arbitrary finite C,P.
Pad `C>=1,P>=0`. For `s=-log q>=1` and `gamma=g_j/2>0`,

\[
\log C+P\log(1+s)-\gamma s
\le\log C+P^2/(2\gamma)-\gamma s/2.
\]

Consequently any log(a_j) exceeding
`2(log C+j log 2)/gamma+P^2/gamma^2` for every listed pair, plus
`log(a_(j-1))+log 2`, `log(2/qstar)` and 1, suffices on the WHOLE interval.
The implementation tests this majorant with extreme synthetic constants
and the CURRENT g_j, not evaluated flow constants. Existence is analytic;
there is no uniform-constant or fixed-q norm-contraction premise.

Define `u_loc=u^[0]+sum V_j`, `p_loc=p^[0]+sum P_j`.
On each compact subset of q>0 only finitely many j are active. Combined
with the locally finite physical-label atlas, this gives ONE smooth real
field, exactly divergence-free, with unchanged core and exterior.
For fixed `J>=max(1,m)` and `q<=1/(2a_J)`, the first J cutoffs equal one
with all derivatives zero. Telescoping (2), then using (6) for j>J, gives

\[
|(u_{\rm loc},p_{\rm loc})-U^{[J]}|_m
\le2^{-J}q^{h\delta(J+1)/2-L_m}.\tag{7}
\]

For inactive j the summand and all its derivatives are exactly zero.
For active j use (6) and `q<1` to sum `2^-j`. Thus (7) holds throughout
the neighborhood, including every transition interval. No convergence
of the uncut series is used; constants like `exp(j^2)` can defeat that.

## 3. Complete residual and the quantifier change

Fix m,N FIRST. Set `(e,pi)=(u_loc,p_loc)-U^[J]`. Exactly,

\[
\mathcal R(U_{\rm loc})-\mathcal R(U^{[J]})
=e_t-\Delta e+\nabla\pi+(u^{[J]}\cdot\nabla)e
 +(e\cdot\nabla)u^{[J]}+(e\cdot\nabla)e.\tag{8}
\]

With `E=|(e,pi)|_(m+2)`, Leibniz bounds its m-jet by
`C_m(E+C_Jm q^-H_(m+2) E+E^2)`. This retains EVERY cutoff curl,
time derivative, pressure gradient, viscosity and nonlinear interaction.
Choose ONE finite J satisfying

\[
J\ge m+2,\quad h\delta(J+1)/2\ge N+1+L_{m+2}+H_{m+2},
\quad hB_J-K_m\ge N+1.\tag{9}
\]

Both left sides grow linearly, so J exists. Equations (1),(7)--(9)
give `|R(U_loc)|_m=O(q^N)` for every fixed m,N, for the SAME field.
J is fixed before the q-limit, and only the comparison neighborhood and
constant depend on m,N,J. This is the justified quantifier change absent
from finite-stage induction alone. Exact flat Gaussian tails occur in
the residual of this one fixed finite comparison state; they are neither
independently summed across stages nor fed into a forward inverse.

At `h=1/200,m=N=0` the rational sufficient-stage calculation is stored
in the evidence. It is an existence budget, not a computation of that
many actual fields or proof that those conservative losses are sharp.

## 4. What ends here, and what does not

The current finite-state-to-local-realization stream reaches its CONDITIONAL
stopping milestone: complete representatives and one all-jet flat-residual
field, with exact protected inner velocity. If the inherited positive
core value is accepted, annular cutoffs cannot erase it. This does not
independently establish that value, endpoint regularity or global energy.

All assertions here concern t<1. Flatness as q->0 is not a proof of smooth
one-sided jets on t=1 away from the singular point, nor an extension to
t>1. The next bounded obligation is to audit that closed-time auxiliary
atlas and all complete potential/pressure jets; then join physical
localization and initialization, every induced commutator, admissible
smooth compact forcing and whole-space energy to THIS realization.
Independent foundation and inverse/geometry review remains mandatory
before any unconditional residual, solution-exclusion or blow-up claim.

## 5. Reproduce

```sh
python3 research/navier_stokes_cascade/results/current_potential_realization/test_audit.py -v
python3 research/navier_stokes_cascade/results/current_potential_realization/audit.py \
  --verify-record research/navier_stokes_cascade/results/current_potential_realization/evidence.json \
  --output /tmp/current_potential_realization_recheck.json
```

The audit replays the current induction, historical realization and frozen
provenance; exact rational checks cover all-stage margins and fixed-J
budgets. Reused manufactured actual-q Taylor-jet controls test the complete
curl and residual identities with missing-term controls. They are algebraic
and numerical diagnostics, NOT a continuum proof or an evaluated global field.
