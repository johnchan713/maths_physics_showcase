# Conditional admissible-stress realization

Status: **conditional-stress-realization-bounded-global-frequency-uninstantiated**.

This checkpoint checks the mechanism in Appendix C of the supplied manuscript:
a prescribed-mean loop of admissible shears, a small radial modulation of the
profile values, and an exact five-moment repair. It supplies an analytic
argument, outward scalar bounds, stable numerical diagnostics, and an explicit
acceptance condition for the modulation frequency. It does **not** instantiate
that condition for the actual joined axis/exterior profile.

In particular, `full_admissible_stress_realized`,
`actual_joined_profile_frequency_selected`, `full_PDE_corrections_verified`,
`smooth_force_verified`, and `blowup_verified` remain **false**. This is neither
a proof-assistant certificate nor independent mathematical peer review.
The preceding [analytic core attachment](../axis_core_attachment/README.md)
and its historical evidence are unchanged. The old `Md=4` outer-cone failure
is retained, not replaced by these manufactured fixtures.

Source: the user-supplied
[Navier-Stokes manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA-256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Relevant printed pages: 29-31, 38, 127-129, 158-165. The PDF page containing
the loop and variance formulas was visually checked as well as text-extracted.
Parent local checkpoint: `6aedb89b2ba2174abbd3237100778a49a5583d80`.

## 1. What must be preserved

For positive shear coordinate `a`, write

\[
t=-b/a,\qquad v=a(1+t^2),\qquad
P_c=p_1+p_2t,\qquad J_c=p_2-p_1t.
\]

The four gaps from (4.35), with the manuscript's sign convention for `b`, are

\[
\Psi=(a,\ v-2,\ P_c-v,\ 2(P_c-v)^2-(v-2)J_c^2).
\]

All four must be positive. Checking only `Pc>2`, or only the fourth polynomial,
does not imply admissibility. The joined profile currently supplies the
relaxed condition, with `v>2` on preserved boundary collars but not everywhere.

The integrated coefficient `p=ps` must be recomputed from the modified
profiles and their cumulative moments, using the same axis pressure datum.
It is not a free parameter in the actual profile. Our scalar loop fixtures
hold it fixed solely to test Lemma C.1; their `p` is not asserted to come from
the original moment formulas.

## 2. Variance without a removable-singularity failure

With normalized auxiliary-angle mean, put

\[
M(z)=\langle e^{z\sin\theta}\rangle=I_0(z),\qquad
t(\theta)=t_s+d_0\frac{e^{\mu p_2\sin\theta}/M(\mu p_2)-1}{p_2}.
\]

Choose `0<d0<min(Pc(ts)-2)/2`. The quotient is smooth at `p2=0`; there it is
`d0 mu sin(theta)`. Exactly,

\[
\langle t\rangle=t_s,\qquad P_c(t)\ge P_c(t_s)-d_0>2.
\]

Directly subtracting `M(2z)/M(z)^2-1` can erase the variance near `p2=0`.
Multiplication of the positive series for `I0` instead gives the identity

\[
V(\mu,p)=d_0^2\mu^2
\frac{\displaystyle\sum_{k\ge1}
 \left(1-\frac{\binom{2k}{k}}{4^k}\right)
 \frac{z^{2k-2}}{(k!)^2}}
 {M(z)^2},\qquad z=\mu p.
\]

Every displayed numerator coefficient is positive. Its value at zero is
`1/2`, so `V(mu,0)=d0^2 mu^2/2` is kept exactly. For the numerator's tail,
bound the coefficient in parentheses by one; the factorial-series ratio is
at most `z^2/(k+2)^2` after the first omitted term. Once this is below one,
the corresponding geometric tail bounds the infinite remainder. The same
argument, with ratio `z^2/[4(k+2)^2]`, bounds `M`. `variance_iv` uses these
outward enclosures, not a finite Taylor polynomial asserted to be exact.

The original monotonicity argument is valid: for `g=log M`, `g''` is the
strictly positive variance of `sin(theta)` under a positive tilted density.
For `p!=0, mu>0`,

\[
\partial_\mu[g(2\mu p)-2g(\mu p)]
=2p[g'(2\mu p)-g'(\mu p)]>0.
\]

At `p=0` strict monotonicity follows from the exact quadratic formula. This
analytic argument, together with the verified signs of a root bracket,
establishes a unique exact variance root between its endpoints. A grid of
monotone samples alone would not establish that conclusion.

### An explicit uniform cap exists

Suppose `a>=amin>0` and `|p2|<=P`. These are input bounds, not values inferred
from seven samples. For `z>=1`, the inequalities
`cos(theta)<=1-2 theta^2/pi^2` on `[-pi,pi]` and
`cos(theta)>=1-theta^2/2` imply

\[
M(z)\le e^z\sqrt{\pi/(8z)},\qquad
M(2z)\ge e^{2z-1}/(\pi\sqrt z),\qquad
\frac{M(2z)}{M(z)^2}\ge\frac{8\sqrt z}{e\pi^2}>\frac{\sqrt z}{4}.
\]

The last constant is checked outward using `pi<22/7`. The positive variance
series also gives `V>=d0^2 mu^2 exp(-2|z|)/2` for all `z`. Consequently set

\[
Z=16\left(1+\frac{8P^2}{a_{\min}d_0^2}\right)^2,\qquad
\mu_{\max}=\frac{4e^Z}{d_0\sqrt{a_{\min}}}.
\]

For `|mu_max p|<=Z`, the latter bound gives `V>=8/amin`. On the complement,
the Bessel-ratio bound gives the same lower bound, using `|p|<=P`.
Thus this one finite cap reaches strictly more than `3/amin` uniformly,
including `p=0`. The implementation records its logarithm. The cap is
deliberately coarse and is not an efficient numerical search prescription.

## 3. All-phase cone margins and the correct mean

An excursion bound avoiding division by `p` is

\[
|t-t_s|\le16d_0\mu\qquad(\mu\ge0,\ p\in\mathbb R).
\]

For `|mu p|<=1`, differentiate `e^{mu p sin}/M(mu p)` with respect to `p`
and use `|g'|<=1`; the quotient is bounded by `2mu e^2<16mu`.
For `|mu p|>=1`, integrating over `|theta|<=1/sqrt(|mu p|)` gives
`e^{|z|}/M(z)<=pi e^(1/2) sqrt(|z|)` and hence a bound
`(pi e^(1/2)+1)mu<16mu`. Both numerical constants are enclosed outward.

Given a verified cap, let

\[
T=t_{\max}+16d_0\mu_{\max},\quad
J=P_2+P_1T,\quad
\gamma=\min(P_c(t_s)-2)-d_0>0.
\]

Choose a positive `delta` no larger than
`min(1/2, gamma/2, gamma^2/[4(1+J^2)])` and smaller than the positive
`vs-2` margin on the preserved radial boundary collars. This order matters:
the variance cap is fixed **before** delta.

Use the manuscript's smooth cutoff `zetaL(vs)` (one below `2+delta/8`, zero
above `2+delta/4`) and put

\[
v_*=2+\delta/2,\quad
\rho=\zeta_L(v_s)^2(v_*-v_s),\quad v=v_s+\rho,\quad
V(\mu,p_2)=\rho/a.
\]

On the cutoff support, `sqrt(rho)=zetaL sqrt(v*-vs)` is smooth because
`v*-vs>=delta/4`; extend it by zero outside. Near `mu=0`, the signed square
root of `V` is smooth and odd with derivative `d0/sqrt(2)`. The implicit
function theorem therefore gives a smooth `mu`, also through zero variance
and through `p2=0`. Away from zero use the strict derivative above.

Where the cutoff is nonzero, the four gaps have the explicit lower bounds

\[
\frac{2}{1+T^2},\qquad \frac\delta8,\qquad
\gamma-\delta/2,\qquad
2(\gamma-\delta/2)^2-\delta J^2/2>0.
\]

In particular the second bound is `delta/8`, **not** `delta/2`, in the
cutoff transition. Where the cutoff is zero, the original shear is unchanged
and `vs>=2+delta/4`. Its strict relaxed condition then implies admissibility.
Compactness on that closed unmodified set gives the remaining positive
minimum, but an actual quantitative certificate must supply that minimum.

The auxiliary-angle vectors have the form `v(1,t)/(1+t^2)`. Their
unweighted mean generally is **not** the requested shear. Define the lift

\[
\frac{d\phi}{d\theta}=\frac{a(1+t^2)}{2\pi v},\qquad
(a_L,-b_L)=\frac{v(1,t)}{1+t^2}.
\]

Since `a<1+t^2>=vs+aV=v`, this positive derivative integrates to one.
Changing variables gives exactly

\[
\int_0^1a_L\,d\phi=a,\qquad
\int_0^1(-b_L)\,d\phi=a\langle t\rangle=-b_s.
\]

The inverse circle map is smooth for fixed compact input bounds. The
preserved boundary collars have zero variance and unchanged shear.

For a separate manufactured family `a=.8, ts=.2, p1=10, p2 in [-1,1]`,
the smaller cap `mu=64` is certified over the **entire** pressure interval:
128 dyadic interval cells on `[0,1]`, together with exact evenness, give
`V>=4.9273>3/.8=3.75`. The outward calculation is repeated at 80 and 110
digits. Here `d0=1, delta=1e-8` yield `aL>1.9066e-6`,
`v-2>=1.25e-9`, `Pc-v>6.79999999`, and fourth gap `>91.95`.
These are manufactured input bounds, not measurements or
certificates of the actual joined profile. Seven numerical representatives
add independent mean and removable-limit diagnostics; they do not replace
the full interval cover.

## 4. Small fields do not mean small radial derivatives

Let the unique zero-mean periodic antiderivatives satisfy

\[
\partial_\phi\mathcal A=-\tfrac12(a_L-a),\qquad
\partial_\phi\mathcal B=\tfrac12E(b_L-b_s).
\]

They vanish on the preserved collars. With phase `phi=N log X`, define

\[
E_N=E e^{\mathcal A/N},\qquad U_N=U+\mathcal B/N.
\]

The exact differentiated shears are

\[
a_N=a_L-\frac{2D_X\mathcal A}{N},\qquad
b_N=e^{-\mathcal A/N}\left(b_L+\frac{2D_X\mathcal B}{NE}\right),
\]

where `DX` on the right holds phase fixed. The axial exponential denominator
cannot be dropped. A closed-form `p2=ts=0` fixture independently differentiates
both modulated fields and detects that omission.

For each fixed angular order `m`, values, moments, and shear-versus-loop
differences are `O_m(1/N)` on the fixed compact interval. Formula (4.17)
controls `ps_N-ps` at order `m` using **one extra eta derivative** of values
and moments. Its denominators include positive lower bounds for `X,E,H,L`.
It uses no comparison of radial derivatives. For radial order `r>=1`,
profile differences may instead be `O_{r,m}(N^(r-1))`.

The phase must be independent of eta for this angular estimate. For example,
`sin(2 pi N(y+eta))/N` has an order-one eta derivative. Tests retain this
countercontrol and the order-one first/order-N second radial derivatives.
The smallness constants still need to be bounded from the actual loop and
joined profile; they are not extracted from the manufactured examples.

## 5. Stable five-moment repair, with the lambda loss retained

Use the **first** reserved interval `(Tw-25,Tw-20)` from (A.9), before the
heat-compensation and two later correction patches. Let its left radius be
`rhop`, put `x=X/rhop`, `y=log x`, and write the unchanged fields as
`U0=0`, `E0=e(eta) x^alpha`, `alpha=-1/2-lambda`, `0<lambda<=.01`.

Choose unit-mass bumps
`beta_c(y)=sigma'((y-c)/.1+.5)/.1`: two U centers `1,2`, three E centers
`.5,1.5,2.5`. All five supports are disjoint and lie inside `(0,5)`.
Set `deltaU=e sum u_j beta_j`, `deltaE=e sum z_j beta_j`.

Normalize the ordinary moment changes as follows. These factors are
different; combining M and J is done only **after** normalization.

| Moment | Dividing factor |
|---|---|
| M | `rhop e` |
| J | `sqrt(2) rhop^(3/2) e^2` |
| I | `sqrt(2) rhop^(3/2) e` |
| S | `rhop e^2` |
| Cp | `e^2` |

Denote normalized rows by `m,j,i,s,c`. Use transformed rows
`(m,(m-j)/lambda,i,-s,c)`. In `dy`, their two linear blocks have weights

\[
B_U:\quad e^y,\ e^y g_\lambda(y),\qquad
g_\lambda(y)=\frac{1-e^{-\lambda y}}\lambda;
\]

\[
B_E:\quad e^{3y/2},\ e^{(1/2-\lambda)y},\ e^{(-1/2-\lambda)y}.
\]

Evaluate `g` with `-expm1(-lambda*y)/lambda`. It has the smooth limit `y`,
and `g'=exp(-lambda*y)>0`. A naive subtraction returns zero for the selected
`lambda=exp[-4(exp(64)+10)]` at the diagnostic precision. The stable formula
does not set this positive lambda to zero.

### Bounds on the exact matrices, without bump quadrature

For the U block, the positive first-row masses and separation between
weighted averages of `g` give

\[
|\det B_U|\ge .9e^{.95+1.95-.0205},\qquad
\|B_U^{-1}\|_\infty\le
\frac{3.05e^{2.05}}{.9e^{.95+1.95-.0205}}<2.
\]

For the E block, at ordered sample points put `x_j=exp(y_j)` and factor
`x_j^(-1/2-lambda)` from column j. The remaining determinant has rows
`(x_j^(2+lambda),x_j,1)`. Its absolute value is the Vandermonde product
times the second divided difference of `x^(2+lambda)`. On `x>=1`, the second
derivative is at least two, so that divided difference is at least one.
Multilinearity and the bumps' exact unit masses give

\[
|\det B_E|\ge e^{-.51(4.65)}
(e^{1.45}-e^{.55})(e^{2.45}-e^{.55})(e^{2.45}-e^{1.55})=:d_E.
\]

Every matrix entry is at most `exp(1.5*2.55)`. Bounding each cofactor by
twice its square gives `||B_E^-1||inf<=6 exp(7.65)/dE<788<1000`, verified
outward at 80 and 110 digits. Thus use **beta=1000** for the transformed
five-row inverse, uniformly for `0<lambda<=.01`. It is independent of eta.

This does **not** remove the lambda loss from arbitrary incoming moments:

\[
\|(\Delta m-\Delta j)/\lambda\|_{C^1}
\le\frac{\|\Delta m\|_{C^1}+\|\Delta j\|_{C^1}}\lambda.
\]

The original U-block determinant is `-lambda det(B_U)`. At lambda zero its
two original rows coincide, and general independent M/J discrepancies cannot
be corrected. The solver explicitly rejects zero lambda, even though the
transformed limiting matrix is perfectly well conditioned. Normalization
and eta derivatives of `e` must also be included in the input error bound.

Stabilizing the matrix cannot recover a small `M-J` already erased by rounding
the incoming moments separately. For the actual profile, that discrepancy
must be carried in a stable integral ledger or bounded before subtraction.
An explicit rounded-input failure is retained; the helper transforming
ordinary numerical targets does not claim to solve this data-loss problem.

### Exact quadratic map and small root

The transformed map is `F(c)=Bc+Q(c,c)`. Its only nonlinear rows are

\[
Q_{-S}=-\sum_{U}u_j^2\int e^y\beta_j^2dy
+\tfrac12\sum_E z_j^2\int e^y\beta_j^2dy,\qquad
Q_{C_p}=\tfrac12\sum_E z_j^2\int\beta_j^2dy.
\]

The J cross term is exactly zero because U/E supports are disjoint; the
U-square term in S is **not** discarded. Since `sigma'<=64`, each unit-mass
bump satisfies `beta^2<=640 beta`. This step bound follows from
`sigma(1-sigma)<=1/4` on `[1/4,3/4]` and the flat exponential estimate on
the reflected endpoint intervals; the resulting constants are checked
outward. Using the maximum of value and first-derivative norms, products
cost a factor two. Hence a valid C1 bilinear norm bound is

\[
2(640)(2e^{2.05}+1.5e^{2.55})<10^5=:q_*.
\]

If the transformed target has `C1` norm at most `d`, then
`8 beta^2 q_* d<=1` gives an exact smooth small root of norm at most
`2 beta d` by contraction. This conclusion concerns exact bump integrals;
the independently refined quadrature roots are only diagnostics of that map.
Higher fixed derivatives are finite by implicit differentiation. There is
no requirement to make infinitely many derivative norms small at once.

At the right of the patch, all five moments agree as functions of eta,
and the profile values are unchanged. Lemma 4.4 then gives exact agreement
of pressure, `Qs,Ns,ps`, and stress at every later radius. Four matching
moments without Cp would leave a pressure offset; that failure is tested.
The original axis pressure datum is never changed to hide such an offset.

## 6. A finite-frequency acceptance condition

Suppose the exact loop and the unmodified compact portion through the repair
patch have gap minimum `kappa>0`, first shear coordinate at least `amin`,
and all four state coordinates `(a,b,p1,p2)` of absolute value at most `B`.
On the segment box for a coordinate error at most `min(1,amin/2)`, set
`R=max(2,B+1,2/amin)`. Then

\[
\|\nabla v\|_1,\ \|\nabla P_c\|_1,\ \|\nabla J_c\|_1\le4R^4,
\quad |v-2|,|P_c-v|\le4R^3,\quad |J_c|\le2R^3.
\]

Differentiating the fourth polynomial bounds its gradient by
`128 R^7+80 R^10<=208 R^10<256 R^10`; the other three gradients obey
the same bound. Thus a sufficient state error tolerance is

\[
\epsilon=\min\left(1,\frac{a_{\min}}2,
\frac{\kappa}{512R^{10}}\right).
\]

Let certified constants `Cstate,D,Ccorr,r>0` have the following meanings:

- modulation state/pressure error is at most `Cstate/N`;
- transformed moment discrepancy has C1 norm at most `D/N`, including
  the dimensional, angular, and **lambda^-1** losses;
- correction state/pressure error is at most `Ccorr ||c||C1`;
- coefficient norm below `r` preserves the required positive-field and
  comparison bounds on the correction patch.

Choose any finite integer N strictly larger than

\[
\max\left(1,\frac{2(C_{\rm state}+2C_{\rm corr}\beta D)}\epsilon,
16\beta^2q_*D,\frac{4\beta D}{r}\right).
\]

The moment contraction has Lipschitz constant at most one half, coefficients
are smaller than `r/2`, and combined state error is smaller than epsilon.
The cone therefore persists. This choice is made **after** all profile,
collar, bump, and lambda choices and **before** physical scale q, later bands,
and PDE correction stages. Radial derivative constants may depend on this
fixed N; they are not claimed to be small.

The preserved inner collar is before both modifications. On the outer side,
exact moment restoration preserves the heat exterior, pressure normalization,
and both still-reserved patches. Appendix C.3 would then retain the original
edge factorizations and their common flat weight. We have not independently
reproved every earlier edge estimate here, nor constructed later waves.

## 7. What this closes, and what it does not

The scalar mechanism has an analytic mean-preserving construction; the exact
repair has a quantitative, nondegenerate transformed inverse and an explicit
small-root condition. Numerical failure at tiny lambda is avoidable, but
the increased frequency requirement is genuine in this general estimate.

The **next actual-profile obligation** is to supply interval endpoints,
derivative bounds, positive compact cone margins, and `Cstate,D,Ccorr,r` for
the selected joined profile, and use them to select N. The universal cap and
frequency formulas establish conditional sufficiency, not those missing
input bounds. Positivity in manufactured data must not promote the global
status flag. After this come the full Navier-Stokes corrections, convergence,
and admissible smooth forcing; none is verified by this checkpoint.

## Reproduction

With the pinned mpmath dependency installed, from this directory:

```bash
python test_audit.py -v
python audit.py --verify-record evidence.json --output reproduced.json
```

`evidence.json` distinguishes outward scalar certificates, exact identities,
numerical diagnostics, inherited source hashes, and unresolved global inputs.
All scientific status flags are checked independently of numerical success.
