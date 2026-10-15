# Conditional full-potential realization and flat local residual

This checkpoint proves an asymptotic summation implication for the COMPLETE
wave and mean corrections in the preceding
[finite-cycle result](../residual_cycle_closure/README.md). Conditional on
that result and its inherited leading-profile, support and auxiliary-chart
hypotheses, there is a real divergence-free local field, smooth for q>0,
whose complete Navier--Stokes momentum residual satisfies

\[
\forall m,N\ge0,\qquad
|\mathcal R(u_{\rm loc},p_{\rm loc})|_m=O(q^N),\qquad q\downarrow0.\tag{1}
\]

Here derivatives are Cartesian space--time derivatives. All cutoff curl
commutators, pressure gradients, viscosity and nonlinear interactions are
included. The annular corrections leave the protected inner velocity
exactly unchanged. A conditional local kinetic-energy bound is also derived.
The proof requires no uniform constants in correction stage and no
convergence of the raw infinite series.

This is a conditional analytic implication, with separate manufactured
diagnostics. The actual infinite sequence of stage coefficients or cutoff
constants is NOT numerically materialized. Independent foundation review,
one-sided endpoint geometry away from the singular point, localization,
initialization and admissible global smooth forcing are still open. No
unconditional solution or blow-up theorem is asserted. In particular, the
conclusions of the proposed manuscript are never used as hypotheses.

## 1. Inputs and their scope

Fix the inherited 0<h<1/100 and put

\[
A=\tfrac12+h,\quad D=\tfrac12-h,\quad \kappa=10^{-5},\quad
\tau=1-t,\quad q=\tau+z^2q^{2h},\quad
\eta=z/q^D,\quad X=r^2/(2q).\tag{2}
\]

For tau>0 this is the positive branch with |eta|<1 and
L=1-2h eta^2>=1-2h>0. All statements concern one fixed local domain
0<q<qstar<1; constants may depend on its inherited fixed parameters.

The preceding checkpoint supplies actual finite states U^[j]=(u^[j],p^[j]),
j>=0, with these CONDITIONAL properties:

1. Each state is real and divergence-free. All increments have smooth
   Cartesian zero extensions supported in the same closed annulus
   Xa<=X<=Xb with Xa>0. They agree under the prescribed chart changes,
   and closed supports of distinct wave labels are separated as required
   by the finite-cycle proof. Mean support is allowed to fill its auxiliary
   torus. The fixed background is the previously summed background.
2. The cycle from j-1 to j is a finite sum of complete wave curls, mean
   azimuthal-potential curls, direct axisymmetric azimuthal velocities and
   pressure increments. Their coefficient/phase bounds hold uniformly
   over all their physical dyadic labels, at every fixed derivative order.
   The same qstar works for every fixed stage; constants may grow arbitrarily.
3. With B_j=7/10+j/10 and K_m=2A+3/2+2m,
   \[
   |\mathcal R(U^{[j]})|_m\le C_{j,m}q^{hB_j-K_m}.\tag{3}
   \]
   This is the COMPLETE pointwise residual, after evaluating the auxiliary
   map and absorbing its fixed-stage flat terms. It includes the background
   stress cancellation, radial pressure equation and remaining axial terms.
   Outside the active supports only the inherited flat base error remains;
   in the heat exterior the error is exactly zero.
4. The leading background is bounded after velocity normalization q^A on
   each fixed compact X range. Its realized inner swirl is
   q^-A(e0+O(q^(2h))), e0>0, at one fixed Xin<Xa, eta=0. Its exterior
   swirl obeys |K(r,tau)|<=C r^(-1-2h). These are inherited conditional
   background/heat conclusions, not new numerical inputs.

The new proof below uses these explicit finite-state estimates. It does
not replace the missing independent review of their construction. No
endpoint extension of the full auxiliary atlas is silently added to this
list: Section 6 separates that remaining question from flatness as q->0.

## 2. Representatives and physical derivative loss

Let U0=U^[0] include the finite initialization. For j>=1 define physical
representatives

\[
\mathcal A_j=\mathcal A_j^w+\Psi_j e_\theta,\quad
\mathcal B_j=b_j e_\theta,\quad \partial_\theta b_j=0,\quad
p_j=p^{[j]}-p^{[j-1]},
\]
\[
\Delta u_j=\operatorname{curl}\mathcal A_j+\mathcal B_j
             =u^{[j]}-u^{[j-1]}.\tag{4}
\]

Use the ACTUAL full potentials, including every wave particular solution,
signed covariance update and both mean operations in that cycle. Their
azimuthal-potential coefficients have the same normalization as the wave
potentials. Pressure differences telescope after all four operations;
intermediate pressure updates are not discarded.

Care with indexing is useful: the cycle j-1 -> j uses
B_(j-1)=3/5+j/10. Its normalized lower orders are

| Entry | Lower order |
|---|---|
| Complete wave velocity increment | B_(j-1)-kappa |
| Wave potential | B_(j-1)+1/2-kappa |
| Mean potential and direct mean swirl | B_(j-1)+1/2-2kappa |
| Pressure increment | B_(j-1)+1/2-2kappa |
| Mean radial velocity increment | B_(j-1)+3/2-2kappa |

Every row exceeds j/10 by a fixed positive margin. Flat realization
remainders fit every specified finite order by choosing their power after
j and the derivative order have been fixed. In a Q chart, q is comparable
to Q. Physical rescalings for velocity, pressure and potentials are
Q^-A, Q^-2A and Q^(1/2-A), respectively, all bounded by a constant times
q^-2A. The preceding physical phase/graph calculation permits a coarse
q^-2 cost per fixed Cartesian space--time derivative; its worst displayed
phase cost 1+5h/2 is strictly below two. Finite harmonic ranges and powers
of the stage number affect constants only.

Set g_j=hj/10. Leibniz' rule therefore gives, for the tuple
Z_j=(A_j,B_j,p_j),

\[
|Z_j|_m\le C_{j,m}(1+|\log q|)^{P_{j,m}}
                q^{g_j-\ell_m},\qquad \ell_m=2A+2m.\tag{5}
\]

One additional spatial derivative reconstructs the velocity. The same
loss must cover the derivative of its summation cutoff. A sufficient
common loss for representatives AND realized velocity/pressure is

\[
L_m=\ell_{m+1}=2A+2m+2.\tag{6}
\]

These losses do not depend on j. For the finite TOTAL states (including
the fixed base), absorbing each fixed logarithmic polynomial by q^-1
gives the deliberately coarse bound

\[
|U^{[j]}|_m\le C_{j,m}q^{-H_m},\qquad H_m=L_m+1.\tag{7}
\]

For example, use (5) at order m+1 for potentials, take their curls, and
use g_j>=0. The fixed base obeys a better derivative bound and fits (7).
This argument does not assume the constants are uniform in j or m.

## 3. Cut potentials and retain their curls

Choose one smooth chi with chi(s)=1 for s<=1/2 and chi(s)=0 for s>=1.
For a>=1 write chi_a=chi(aq). Since q is independent of r and theta,

\[
\operatorname{curl}(\chi_a\mathcal A_j)
 =\chi_a\Delta u_j-\chi_a\mathcal B_j
   +\nabla\chi_a\times\mathcal A_j.\tag{8}
\]

Thus the cut velocity is chi_a Delta u_j + grad chi_a cross A_j.
The additional MEAN radial component is specifically

\[
(\nabla\chi_a\times\Psi_j e_\theta)_r
       =-a\chi'(aq)q_z\Psi_j.\tag{9}
\]

All components of the wave commutator in (8) are retained as well. Cutting
the velocity directly would give div(chi_a Delta u_j)=grad chi_a dot
Delta u_j, which need not vanish. Cutting b_j e_theta directly is valid:
div(chi_a b_j e_theta)=r^-1 partial_theta(chi_a b_j)=0.

Differentiating (2) gives smooth bounded eta coefficients and

\[
|\partial_z^a\partial_t^b q|\le C_{a,b}q^{1-Da-b}
\quad(a+b\ge1).
\]

This follows by implicit differentiation and induction with L bounded
below; equivalently the weighted homogeneity of (2) gives the powers and
the closed eta interval bounds the coefficients. In a positive derivative
of chi(aq), every term is chi^(k)(aq) a^k times k differentiated factors
of q. On its support 1/2<=aq<=1, the q^k from those factors cancels a^k.
Consequently

\[
|\partial_z^a\partial_t^b\chi(aq)|
       \le C_{a,b}q^{-Da-b}\le C_{a,b}q^{-a-b},\tag{10}
\]

with constants INDEPENDENT of a. Radial and angular derivatives vanish.
Combining (5), (10) and the one curl proves

\[
|\chi_a Z_j|_m+
|\operatorname{curl}(\chi_a A_j)+\chi_a B_j,\chi_a p_j|_m
\le \widetilde C_{j,m}(1+|\log q|)^{\widetilde P_{j,m}}
                      q^{g_j-L_m},\tag{11}
\]

uniformly in the cutoff scale. One may use (5) at order m+1 and the coarse
two-per-derivative bound throughout; this includes every split of the
derivatives in the product rule. Smooth zero extensions at support edges
are part of the inherited representatives. Multiplication and a local
curl cannot enlarge their closed supports. For each t<1 the annulus stays
away from r=0, so its zero extension is smooth at the Cartesian axis.

## 4. A diagonal scale choice on whole intervals

At step j collect the finitely many constants in (11) for 0<=m<=j. Pad
C>=1, P>=0. Choose a_j>=2a_(j-1) large enough that, for EVERY requirement,

\[
C(1+|\log q|)^Pq^{g_j/2}\le2^{-j}
                    \quad(0<q\le a_j^{-1}).\tag{12}
\]

Such a scale exists even for arbitrarily large C and P. For an explicit
sufficient recipe put gamma=g_j/2 and s=-log q. The elementary inequalities
log(1+s)<=sqrt(s) and P sqrt(s)<=gamma*s/2+P^2/(2gamma) give

\[
\log C+P\log(1+s)-\gamma s
 \le \log C+P^2/(2\gamma)-\gamma s/2.
\]

It suffices that log a_j exceed

\[
\max\left\{\log a_{j-1}+\log2,\log(q_\star^{-1})+\log2,1,
 \max_{\rm requirements}
 \left[\frac{2(\log C+j\log2)}\gamma+\frac{P^2}{\gamma^2}\right]\right\}.
\tag{13}
\]

The implementation reuses the independently tested background recipe in
log coordinates and adds strict padding. It tests SYNTHETIC constants;
actual Navier--Stokes constants are not supplied or evaluated. The
inequality above proves the requirement for all smaller q, not merely
the finite diagnostic samples. No effective bound on their growth is
needed for the existence implication.

Define

\[
u_{\rm loc}=u^{[0]}+\sum_{j\ge1}
 \{\operatorname{curl}(\chi(a_jq)A_j)+\chi(a_jq)B_j\},\quad
p_{\rm loc}=p^{[0]}+\sum_{j\ge1}\chi(a_jq)p_j.\tag{14}
\]

On any compact set with q>=c>0 only finitely many stage cutoffs are
nonzero. The inherited dyadic labels are locally finite there too.
Therefore (14) is smooth for q>0, real, and EXACTLY divergence-free.
It has the same protected core and heat exterior. The fixed base has
its own already realized potentials, so (14) also supplies full local
potential representatives. There is no raw infinite-series limit being
interchanged with differentiation.

Fix m and J>=max(1,m). For q<=1/(2a_J), the first J cutoffs equal one,
with all their positive derivatives zero. For j>J apply (12) to (11).
Since q<1, summing the geometric factors gives

\[
|U_{\rm loc}-U^{[J]}|_m
 \le 2^{-J}q^{h(J+1)/20-L_m}.\tag{15}
\]

This bound covers pressure, potentials and actual velocity, including all
commutators. J is fixed BEFORE taking q small. The raw constants C_j,m
may grow too fast to permit substituting j(q) into any uncut bound.
For example coefficients exp(j^2) defeat a raw sum q^(hj/10) at every
fixed q>0; they do not obstruct (12).

## 5. The complete equation error is flat

Fix a Cartesian derivative order m and a desired power N>=0. Set
e=u_loc-u^[J] and pi=p_loc-p^[J]. For viscosity one the EXACT identity is

\[
\begin{aligned}
\mathcal R(U_{\rm loc})-\mathcal R(U^{[J]})
={}&\partial_te-\Delta e+\nabla\pi
 +(u^{[J]}\cdot\nabla)e+(e\cdot\nabla)u^{[J]}
 +(e\cdot\nabla)e.\tag{16}
\end{aligned}
\]

Use ordinary Cartesian product estimates on jets through m+2. With
E=|(e,pi)|_(m+2), (7) gives

\[
|\mathcal R(U_{\rm loc})-\mathcal R(U^{[J]})|_m
 \le C_m\{E+C_{J,m}q^{-H_{m+2}}E+E^2\}.\tag{17}
\]

Choose ONE integer J>=m+2 so that

\[
\frac{h(J+1)}{20}\ge N+1+L_{m+2}+H_{m+2},\qquad
h(7/10+J/10)-K_m\ge N+1.\tag{18}
\]

Both left sides grow linearly in J, with positive slopes independent of
m. Equation (15) makes all three terms of (17) O(q^(N+1)); equation
(3) does the same for the finite-state residual. Restricting q below
1/(2a_J) is allowed after this fixed choice. This proves (1), including
ALL derivatives, viscosity, cutoff terms, radial/axial pressure and every
old/new and new/new product.

`comparison_budget` gives exact rational sufficient choices. For the
diagnostic h=1/200, m=N=0 it returns J=64079; this is a deliberately
coarse EXISTENCE budget, not the evaluation of 64079 stages. It does not
measure how many stages a numerical approximation would require. The
finite-stage flat remainders are handled in (3) at this one fixed stage.
They are never independently summed across j.

## 6. Singular corner, protected growth and local energy

All fixed residual derivatives tend to zero as q->0. Moreover (2) implies

\[
q\le2\tau+(2z^2)^{1/(1-2h)}
 \le C(\tau+|z|)\le C(|x|+|t-1|)\quad\hbox{near }(0,1),\ t\le1.\tag{19}
\]

Indeed either q<=2tau or the second term in (2) is at least q/2.
Thus every residual jet is flat also in ordinary space--time distance
when approaching the singular corner from below. In the exterior the
residual is already zero. If the inherited finite atlas has compatible
smooth one-sided limits on q>=c>0, local finiteness supplies those limits
for (14), and (1), (19) attach the zero residual jet at (0,1). This gives
a one-sided smooth residual through the corner. The endpoint atlas and
the subsequent extension to t>1 have NOT been independently verified in
this checkpoint. Flatness at the corner alone is not a global smooth force.

At Xin<Xa all new potentials, direct means and pressures vanish in a
neighborhood. The cutoffs cannot change this. Hence at
x_tau=(sqrt(2 Xin tau),0,0), t=1-tau,

\[
(u_{\rm loc})_\theta(x_\tau,1-\tau)
       =\tau^{-A}\{e_0+O(\tau^{2h})\},\qquad e_0>0.\tag{20}
\]

For sufficiently small tau its magnitude is at least (e0/2)tau^-A.
This is conditional on the inherited positive leading value and
background sum, but now includes the ENTIRE realized annular correction.

There is also useful size control away from that path. Every fixed
finite-stage annular correction has normalized velocity of order at
least 1/2, so
|u^[J]-u_base|<=C_J q^(-A+h/2)(1+|log q|)^P.
Absorb the logarithm by q^(-h/4). In (15), choose a fixed J large
enough that its velocity tail exponent exceeds -A+h/4. It follows that

\[
|u_{\rm loc}-u_{\rm base}|\le Cq^{-A+h/4}
\quad\hbox{on the correction annulus},\qquad
|u_{\rm loc}|\le Cq^{-A}\quad(0\le X\le X_{\rm ext}).\tag{21}
\]

Outside this range only the inherited exterior swirl remains. For a
fixed small cylinder r<R0, |z|<Z0, split the radial energy integral at
r=C sqrt(q). Equation (21) and |K|<=C r^(-1-2h) give

\[
\int_0^{R_0}|u_{\rm loc}|^2 r\,dr\le Cq^{1-2A}=Cq^{-2h}.
\]

If the splitting radius exceeds R0, integrate (21) only to R0; the same
upper bound still holds. Since |eta|<=1 implies q>=|z|^(1/D),

\[
\sup_{0<\tau<\tau_0}\int_{|z|<Z_0,r<R_0}|u_{\rm loc}|^2\,dx
 \le C\int_{-Z_0}^{Z_0}|z|^{-2h/D}\,dz<\infty,
\quad 2h<D.\tag{22}
\]

The last strict inequality holds for every allowed h since D-2h=1/2-3h
>0. The energy in an axial slab |z|<delta^D tends to zero uniformly like
delta^(D-2h). This is LOCAL kinetic energy, not a proved whole-space
admissible solution, energy identity, dissipation bound or uniqueness
argument. Those require the later physical localization and initialization.

## 7. Reproducible checks and the next gate

The diagnostics use the actual implicit q in (2) but manufactured compact
wave/mean potentials. Third-order Taylor arithmetic checks exact curls
and their velocity derivatives through order two, with independent
Cartesian/cylindrical equation evaluation. Controls detect a missing or
reversed curl commutator, cutoff time terms, pressure cutoff gradients,
viscous commutators and each nonlinear product. Exact rational budgets
check the stage indexing, all-stage margins, fixed-J comparison and energy
integrability. Synthetic extreme constants test the all-interval scale
recipe. These checks audit the algebra and scope; they do not simulate
the actual global pulses or independently establish the inherited premises.

```sh
python3 research/navier_stokes_cascade/results/full_potential_realization/test_audit.py -v
python3 research/navier_stokes_cascade/results/full_potential_realization/audit.py \
  --verify-record research/navier_stokes_cascade/results/full_potential_realization/evidence.json \
  --output /tmp/full_potential_realization_recheck.json
```

Next: check the closed-time auxiliary atlas and one-sided potential jets
away from (0,1); localize and initialize the COMPLETE potential field;
prove the resulting residual extends to an admissible smooth compact force
through t=1. This must include every localization curl, pressure, viscosity
and nonlinear term. Then check whole-space energy and the exact argument
excluding a global smooth solution with the same data. Independent review
of the leading profile, background and global support geometry remains a
parallel requirement. Earlier frozen evidence is unchanged.
