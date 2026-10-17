# Local core foundation: all orders, reciprocals and physical scales

This checkpoint starts the independent foundation pass requested after
[current_solution_comparison](../current_solution_comparison/README.md).
It quantifies and checks the LOCAL analytic-core step directly from the
infinite coefficient norm and defining similarity map. It supplies an
explicit joint all-order majorant, one zero-free reciprocal neighborhood,
and Cartesian axis/physical-scale estimates for the stated pressure-family
core. The finite Taylor iterates are NOT used as the infinite solution.

Parent: `19e184739d019e14c3e151a9f3a49dbe7a2d9f65`. The defining inputs are
the infinite-space construction in
[core_contraction_audit](../core_contraction_audit/README.md), the reference
pressure and frozen parameters in
[core_growth_activation](../core_growth_activation/README.md), and the
actual Cartesian ansatz in [paper_profile_audit](../paper_profile_audit/README.md).
Their frozen records, source bytes, constants, pressure memory and width
C^-20 are unchanged. This audit checks the norm-to-field and physical
derivative steps; it is not a new formal or peer-review certificate for the
entire nonlinear/global construction.

The pressure-family hypothesis remains explicit: P>=16, 0<h<=.01,
Pi_axis=-P^2 int_(0,1)(1+eta^2)^(-2theta) dmu(theta), mu positive with
mass<5. The original reference has the additional positive atom used for
growth. No sampled pressure fit replaces this datum. Identification with
the completed global backward pressure remains a separate defining-map
obligation, not a consequence of local analyticity.

## 1. Recheck the infinite-space input, then sum its majorant

For H=sum_n H_n(eta)Y^n the retained norm gives, for EVERY n,k,

\[
|\partial_\eta^kH_n|\le B\,
\frac{20^{-n}\rho^{-k} k!\binom{n+k}{k}}{(n+1)^2(k+1)^2}.\tag{1}
\]

Here rho=sigma/2048 and B is the norm, not a finite truncation error.
The parent proof's product inequality follows by counting k-subsets in two
disjoint blocks, its two squared-denominator convolutions have product
constant 256, and the shifted derivative ratio is at most
80(i+1)/rho. In a J-integrated mixed monomial the remaining factor is
(i+1)j/[(i+j+1)(i+j+nu)]<=1. Thus no bounded bare eta-derivative operator
is assumed. The factorial inverse series for the angular operator and all
three pressure-source terms remain included. The parent contraction,
original-equation and failure controls are regenerated with their original
comparator, in addition to the current downstream record.

Expand the coefficient functions in Taylor disks centered at real
eta_0 in [-1,1]. Put a=|Y|/20 and b=|eta-eta_0|/rho. Dropping only positive
denominators in (1), the exact generating function is

\[
\sum_{n,k\ge0}\binom{n+k}{k}a^n b^k=\frac1{1-a-b},\qquad a+b<1.
\]

The disks glue on overlaps by their common real analytic coefficients,
including disks centered at the endpoints. Absolute convergence on every
smaller joint domain permits termwise mixed differentiation. Differentiating
the positive majorant gives the explicit all-order estimate

\[
|\partial_Y^r\partial_\eta^k H|
\le B\,(r+k)!20^{-r}\rho^{-k}
 (1-|Y|/20-|eta-eta_0|/\rho)^{-r-k-1}.\tag{2}
\]

The mixed factorial is (r+k)!, NOT r!k!. Nor may the two convergence
conditions a<1 and b<1 replace a+b<1. The code has exact counterchecks
for both shortcuts. Formula (2) proves finiteness for arbitrary fixed
orders; tests of selected orders do not prove its quantifiers.

## 2. An explicit common zero-free neighborhood

For the actual infinite core, ||Phi||<=S, S>=2, and the parent real bound
on 0<=Y<=4.1 is Phi>33/125. On |Y|<=5 and complex eta displacement
<=rho/8 the joint margin is at least 5/8. In particular (2) gives

\[
|\Phi_Y|\le(16/125)S,\qquad
|\Phi_\eta|\le(64/25)S/\rho.
\]

Around EVERY real point of [0,4.1] x [-1,1], choose the SAME polydisc radii
a0=1/(100S), b0=rho/(100S). These disks stay strictly inside the above
majorant domain. Integrating the two complex derivatives along coordinate
segments changes Phi by at most

\[
\frac{16/125+64/25}{100}=84/3125.
\]

Therefore Re Phi>741/3125>1/5 throughout their union. In particular
|1/Phi|<5 there. Cauchy's formula on each polydisc supplies

\[
|\partial_Y^r\partial_\eta^k(1/\Phi)|
\le5r!k!(100S)^r(100S/\rho)^k.\tag{3}
\]

The principal log is holomorphic on the same union, since Re Phi>0.
Its supremum, hence its fixed-order Cauchy constants, is finite and may
depend on S. Real positivity by itself would NOT give a common complex
strip: 1+n^2 eta^2 is >=1 on the real line but has zeros at +/-i/n.
No jet constant uniform in derivative order or as C grows is asserted.

## 3. Full pressure, radial averages and regular quotients

Keep g=phi_star/C with the original nonvanishing holomorphic exponential.
The selected C gives ||g||<=1, while |1/g|<=C exp(2 Lambda M) on the
retained complex rectangle. From the SAME Banach algebra and radial
integration estimate, the actual fixed core satisfies

| Field | Coefficient norm upper bound |
|---|---:|
| F=g Phi | 256 S |
| U=U_star+u/Lambda | 9+S/Lambda |
| p, with p_Y=g^2 Phi^2 | 80*256^3 S^2 |
| Pi=Pi_axis+p/Lambda | 10P^2+80*256^3 S^2/Lambda |

Thus (2) bounds every fixed mixed derivative, retaining the datum AND
pressure increment, including pressure eta derivatives. The reciprocal
of F is bounded on the smaller neighborhood by 5C exp(2 Lambda M), and
all its fixed jets follow by Cauchy. These enormous finite constants are
not numerically materialized or declared order-independent.

The average Ubar=int_0^1 U(sY,eta)ds has the same neighborhood; fixed
radial derivatives pass under the integral with a factor s^r<=1. Put
A=1/2+h, D=1/2-h, d=1-eta^2 and L=1-2h eta^2. The regular quantity

\[
v_0=\frac{2eta U-2Deta\,\overline U-d\partial_\eta\overline U}{L}
\]

has every fixed jet on a smaller neighborhood by (2), with one additional
parameter derivative. This does not make bare eta differentiation bounded
in the original Banach norm. On the real closed interval L>=49/50; its
complex reciprocal is controlled after a fixed further shrink if needed.
E=sqrt(2X)F vanishes on the axis: E^-1 is NOT an axis-regular quantity.
F^-1 and the displayed v0 are the regular objects.

## 4. Exact physical derivative operators, including the axis

Define q by q-z^2 q^(2h)=1-t, and xi=x/sqrt(q), upsilon=y/sqrt(q),
eta=z/q^D. At positive q the implicit denominator is L>=49/50, even at
eta=+/-1. For any scalar profile G and exponent beta, direct differentiation
gives

\[
\begin{aligned}
\partial_x(q^{-\beta}G)&=q^{-\beta-1/2}G_\xi,\\
\partial_y(q^{-\beta}G)&=q^{-\beta-1/2}G_\upsilon,\\
\partial_z(q^{-\beta}G)&=q^{-\beta-D}Z_\beta G,\\
\partial_t(q^{-\beta}G)&=q^{-\beta-1}T_\beta G,
\end{aligned}
\]

where
Z_beta=[d partial_eta-eta(xi partial_xi+upsilon partial_upsilon)-2beta eta]/L,
and T_beta=[beta+(xi partial_xi+upsilon partial_upsilon)/2+Deta partial_eta]/L.
The exponent beta MUST be updated after EACH derivative. Freezing it
already fails at the second time derivative, as an exact control detects.

Use Y=Lambda(xi^2+upsilon^2)/2 and the actual Cartesian ansatz

u1=q^-1/2 xi v0/2-q^-A upsilon F,
u2=q^-1/2 upsilon v0/2+q^-A xi F,
u3=q^-A U, p_phys=q^-2A Pi.

There is no division by sqrt(x^2+y^2). All profile coefficients are smooth
functions of xi,upsilon and eta, including xi=upsilon=0. Applying each
operator preserves bounded fixed profile jets on the compact normalized
core. An induction, with shifted beta, proves for every alpha,j and
0<q<=1,

\[
|D_x^\alpha\partial_t^j u|\le C_{\alpha j}q^{-A-w},\quad
|D_x^\alpha\partial_t^j p_{phys}|\le C_{\alpha j}q^{-2A-w},\quad
w=(\alpha_x+\alpha_y)/2+D\alpha_z+j.\tag{4}
\]

The radial term initially has exponent 1/2; q<=1 and A>1/2 let it be
dominated by the stated velocity bound. Constants depend on the frozen
core parameters and requested derivatives, NOT on q. This proves local
all-physical-scale size/derivative bounds, not a bound uniform in C,N or
all derivative orders at once. It does not extend velocity to q=0.

## 5. Scope boundary and the precise next inequality

The unmodified analytic SEED exists through Y=4.1, which supplies room for
the complex estimates. The FINAL leading profile is guaranteed unchanged
only through Y=4, X<=Xa=4/Lambda. Its actual flat activation starts to the
right. Formula (4) is claimed for that unchanged leading core, NOT silently
for the whole finalized profile through Y=4.1. It does not imply vanishing
complete residual: leading stress cancellation omits later viscous and
interaction obligations.

The next common-collar certificate must be derived for the ACTUALLY
activated leading functions, not their analytic-core reference. With
b^2=Xa exp(tc/4), for EVERY fixed r it requires some delta_r,B_r>0 such that

\[
\sup_{0\le X\le b^2,\;dist(eta,[-1,1])<\delta_r}
\sum_{H\in\{F_0,U_0,Pi_0,M_0/X,V_0/X\}}|\partial_X^r H|
\le B_r,\tag{5}
\]

with holomorphy in eta and the regular axis values. This audit supplies
the left/core part of (5); it does NOT yet supply its part in Xa<X<=b^2.
The actual flat shear prescription and forward moments/pressure must be
used to prove that part. Fixed-width radial derivatives may cost powers
of C^20; no width-independent constants are allowed. Finite C2/C6
closeness or replacing the activated functions by their reference cannot
prove (5) at arbitrary order.

After (5), independently audit the modulation/repair, global pressure
moments and weighted edges, then the background and actual pulse inverse
estimates. The global conditional exclusion chain remains conditional.
No unconditional pulse bound, complete residual improvement or blow-up
claim is made here.

## Reproduce and interpret the diagnostics

The exact rational diagnostics solve the implicit q equation as a
four-variable Taylor jet, and independently apply normalized derivative
words. All 70 multiindices of total order <=4 for all four fields give
280 exact identities, including mixed orders. Separate axis and off-axis
divergence, nonzero full residual, pressure and axial-viscosity controls
prevent interpreting a manufactured profile as the actual core/solution.
Their h=1/200, Lambda=3 are fixtures, never replacements for the selected
parameters. The analytical proofs are (1)--(4), not those finite checks.

```sh
python3 research/navier_stokes_cascade/results/core_all_order_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_all_order_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_all_order_foundation/evidence.json \
  --output /tmp/core_all_order_foundation_fresh.json
```

Evidence binds 13 current source/record files and 758 inherited source
checks. No historical source, evidence or tolerance is rewritten.
