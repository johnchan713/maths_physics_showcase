# Axial moments and the complete reference strength root

Starting from `a41b941add8dc9c0368e3cea7c4f0aed979a1f78`, this checkpoint
independently closes the retained reference `M,J` solve and scalar
`S(infinity)=0` strength solve. It proves

\[
 \lambda\|B^{-1}\|_\infty<4,\qquad
 1.0100502<A(\eta)<1.0100504,\qquad
 F_A(A(\eta),\eta)>.49,\qquad |A_\eta|<T\lambda/30.
\]

The inverse improves the previous `20/lambda` bound and the last estimate
improves `3000 T lambda`. The SAME accepted roots and parameters are used.
The complete energy equation retains its nonzero affine end-bump terms.
All fixed one-sided eta jets exist and have the explicit recursion below.
This does not discharge the annular/compact-repair maps, complete-core
backward-pressure identification or remaining finite-background/pulse lemmas.
There is no unconditional complete residual or blow-up conclusion.

The defining maps are
[`pulse_moment_audit/pulse.py`](../pulse_moment_audit/pulse.py) and
[`outer_pressure_pilot/schedule.py`](../outer_pressure_pilot/schedule.py).
The previous [moment audit](../pulse_moment_audit/README.md),
[angular/heat root foundation](../root_foundation/README.md) and all historical
frozen records remain unchanged. The new checkpoint distinguishes a full
analytic derivation, outward enclosures and finite numerical diagnostics.

## 1. Parameters, moments and the actual supports

Keep `Md=64`, `T=exp(64)+10>=64`, `lambda=exp(-4T)`, `h=lambda^2`,
`Tw=240T`, `Tf=1000`, `co=1/1000`, unit viscosity and width `C_setup^-20`.
Let `f=(1+eta^2)^-1`, `Ep=eb f`, and use pulse logarithmic coordinate y.
The prescribed fields have

\[
 E=E_pe^{-(1/2+\lambda)y},\qquad
 U=E\left[A(\eta)R_0(\lambda y)+c_1\beta_1(y)+c_2\beta_2(y)\right],
 \quad 0\le y\le13/\lambda.
\]

Here `R0(xi)=phi(xi)[1-sigma(xi-10)]`,
`phi(xi)=integral_0^xi sigma(v/.02)dv`; thus `0<=R0<=11`, supported up to
`11/lambda`. The width `.3` bumps are centered at `Y=13/lambda-3` and `Y+2`.
Their supports are disjoint from the main pulse and each other. Their
positive masses are `.3`; `sup sigma'<=64` is supplied by the parent proof.
Both supports lie at positive y. Consequently

\[
 0<w_j:=\int e^{-2\lambda y}\beta_j(y)^2\,dy<.3\cdot64<20.
\]

No bump is removed because its coefficient is very small.
The moment conventions are exactly

\[
 M=\int U\,dX,\quad J=\int UH\,dX,\quad
 S=\int(U^2-E^2/2)\,dX,\qquad H=\sqrt{2X}E.
\]

The factors of `sqrt(2)` cancel between J and its incoming normalization. Write
`m0=M/(Xp Ep)`, `j0=J/(Xp Hp Ep)`, `s0=S/(Xp Ep^2)`.
For coefficient functions use
`max_j(sup|cj|+sup|cj_eta|)`, including one-sided endpoints.

## 2. An independent prefix bound

The exact amplitude scale is

\[
 \log e_b=-119.5T+.3-\lambda(240T+.5).
\]

The incoming positive integrals use `0<=k(t)<=4` on `0<=t<=e^Md-1`:
`Ka=4e^-T+integral e^(t-T)k(t)dt` and the analogous K2a with initial 16
and k squared. Hence `0<Ka<=4`, `0<K2a<=16`.
For `r=I/(XH)`, the first ramp has `r'=1-(1+l)r`, `r(0)=5/8` and
`0<=l<=.6`, so `0<r1<=1`. The exact
`Ja=Ka+4(r1-1)e^-T` satisfies `0<Ja<=4`, because `Ka>=4e^-T`.
These are analytic bounds at Md64; the executable `Prefix` class's Md4
fixtures are not being promoted to an Md64 numerical profile.

Putting `g=lambda(Tw+.5)`, the exact normalized incoming moments are

\[
 m_0=\eta(1+\eta^2)K_ae^{-1-T_w}/e_b,\quad
 j_0=\eta(1+\eta^2)J_ae^{-1-T_w+g}/e_b.
\]

The polynomial has full C1 norm six. Since `g<=1`, both full C1 norms are
at most `24 exp(-120.5T+.7)<1`. The positive U-squared part of s0 has
value at most `64 exp(-T+.4)<1` and eta derivative at most
`256 exp(-T+.4)<1`; its normalized E-squared part is eta-independent.

For an independent bound on that E memory, its first-ramp normalized mass
G obeys `G'=1-2lG<=1`, `G(0)=5/6`, so `G1<=11/6`.
This coarser bound avoids needing the earlier `G1<1.5` estimate.
The following stages, including the whole reserved interval, give

\[
 G_p\le e^{2g}(G_1+T+1+T_w)\le2(241T+17/6).
\]

Here `g<1/4` and `e^(1/2)<2`. Therefore
`|s0|<=241T+23/6<500T`, `|s0_eta|<1`. No earlier energy is discarded.

## 3. Both axial moments: keep the small determinant factored

Let `s1=1/2-lambda`, `s2=1/2-2lambda`, and
`b(s)=integral exp(s t) beta(t)dt` for the centered width `.3` bump.
Normalizing both rows at Y gives the actual system

\[
 Bc=d,\quad B=\begin{pmatrix}b_1&b_1n_1\\b_2&b_2n_2\end{pmatrix},
 \quad n_i=e^{2s_i},
\]
\[
 d_i=-m_{i,0}e^{-s_iY}
       -A\int_0^{11/\lambda}e^{s_i(y-Y)}R_0(\lambda y)dy,
 \quad(m_{1,0},m_{2,0})=(m_0,j_0).
\]

For `0<lambda<=.01`, `.48<=si<=.5` and
`b_i>=b_min=.3 exp(-.075)`. More precisely,

\[
 n_1-n_2=e^{1-4\lambda}\operatorname{expm1}(2\lambda)
             \ge2\lambda e^{.96}>0.
\]

The inverse row sums satisfy

\[
 \lambda\|(B^{-1})_{0,:}\|_1
 \le\frac{1+e^{.02}}{2b_{min}}<3.629239,\qquad
 \lambda\|(B^{-1})_{1,:}\|_1
 \le\frac{e^{-.96}}{b_{min}}<1.375714.
\]

Thus the scaled inverse norm is below four uniformly. The UNscaled matrix
at `lambda=0` is singular; no inverse at that endpoint is asserted.
The continuum bounds use rational exponential enclosures from positive
Taylor sums, reciprocal bounds and outward rounding, not sampled inverses.

The main pulse lies at least `2/lambda-3` before Y. Using `R0<=11`,
`si>=.48`, and the full prefix C1 bounds proves, for `0<=A<=1.2`,

\[
 \|d\|_{C^1},\ |d_A|\le
 e^{1.5}\left(1+1.2\cdot11/.48\right)e^{-.96/\lambda}
 <200e^{-.8/\lambda}.
\]

The same bound controls `sup_(A,eta)|d|+sup_(A,eta)|d_eta|` over this
whole trial domain: the main part is independent of eta and bounded using
`A<=1.2`, while the incoming part already has its separate value and
derivative bounds. This stronger uniform form also passes through the
constant matrix inverse. Therefore

\[
 \|c(A,\cdot)\|_{C^1},\ |c_A|
 \le800\lambda^{-1}e^{-.8/\lambda}\le\lambda^{20}.
\]

The last comparison is equivalent to
`lambda(log 800+21 log(1/lambda))<=.8`. Its left side increases with
lambda in the selected range, so its maximum occurs at `e^-256` and is
less than `(7+21*256) 2^-352`. The binary bound follows from
`e^8>2^11`, proved by a positive Taylor sum, and `log 800<7`.
This evaluates neither `exp(-.8/lambda)` nor the huge support coordinate.
The matrix is constant in eta, so the SAME solve gives every fixed jet.
In particular `c_j=a_j A+b_j(eta)`, with constant a_j and
`b_j(eta)=v_j eta(1+eta^2)`. M and J close exactly for every trial strength.

## 4. The complete post-pulse E energy

Set `R=integral_(Xend)^infinity E^2 dX/(Xend Eend(eta)^2)` for the angularly
corrected POWER reference, before the later Gamma heat attachment.
Interpolation contributes at most Tf=1000. The uniform interval contributes
at most `2Tu`, `Tu=120T`, including both signed angular energy increments.
This uses the parent's actual angular root and relative edit bound `.01`.
The first release transition contributes at most one, the steep `l=-1`
hold at most one-half, and that hold reduces `XE^2` by exactly h^8.
All remaining release, waiting, terminal-cutoff and infinite-tail stages have
`l<=-3h/4`: the nonzero terminal bump contributes at most
`32 co h/(1-co h)<h/4` to l. The entire remaining mass is therefore
at most `(2/3)h^7`. Consequently

\[
 0<R\le1002+240T\le256T.
\]

All endpoints and the reference wait are eta-independent. The normalized
squared interpolation/post shape has eta logarithmic derivative at most
two, plus `2*.01/(1-.01)` on an edited patch; four is a common bound.
Differentiation of the complete integrals gives `|R_eta|<=4R<=1024T`.
The energy of a very long waiting stage is bounded by its decay integral,
not by discarding the stage or replacing it with a finite-tail truncation.
The parent proves the scalar reference wait and the actual angular roots.

The Gamma heat attachment changes the pointwise tail. Its three-bump
compensation has exact net Cp, S and renormalized-I change zero as eta
functions, with their fixed jets, by the parent foundation. Thus the COMPLETE
S equation below remains the same after that attachment. Both edited regions
have U=0, so M and J are unchanged as well. Its pointwise
power-tail formula is not assigned to the final heat-modified field.

## 5. Full scalar quadratic, signs and its positive root

Disjoint supports give the exact normalized equation

\[
 F(A,\eta)=\frac{\lambda S(\infty)}{X_pE_p^2}
 =A^2K-C+\lambda s_0-\frac\lambda2e^{-26}R
                          +\lambda\sum_jw_jc_j(A,\eta)^2=0,
\]

where `K=integral_0^11 exp(-2xi) R0(xi)^2 dxi`, `C=(1-exp(-26))/4`.
The affine corrections make this exactly `k A^2+b(eta) A+d(eta)`, with

\[
 k=K+\lambda\sum w_ja_j^2,\quad b=2\lambda\sum w_ja_jb_j,\quad
 d=-C+\lambda s_0-\tfrac\lambda2e^{-26}R+\lambda\sum w_jb_j^2.
\]

All three bump-energy increments are retained. The cross increment is
generically odd in eta, so evenness of A must not be assumed.
The trial domain includes A=0; this also justifies the bound on d itself.
For `0<=A<=1.2`, the preceding estimates give

\[
 |\mathcal E|,|\mathcal E_\eta|\le1000T\lambda<10^{-100},\quad
 |\mathcal E_A|\le80\lambda^{41},\qquad F=A^2K-C+\mathcal E.
\]

For example the eta error has the more useful bound
`lambda+512 T lambda exp(-26)+80 lambda^41`.
The scalar product `T exp(-4T)` decreases for `T>=64`; hence
`1000 T lambda<64000*2^-352<10^-100`. These are exact rational comparisons.

The actual positive monotone enclosure for K is independently recomputed
with 1024 boxes and 80 interval digits using the existing defining integral:

\[
 .2450496182718<K<.2450496221363.
\]

For the initial interval `[0,.02]` use left/right sums for the nondecreasing
step primitive; the middle `[.02,10]` is integrated exactly using
`phi(xi)=xi-.01`; the final `[10,11]` uses monotonicity of the cutoff and
exact positive weights. The stored rational endpoints strictly enclose the
outward interval endpoints, converted to exact dyadic rationals. Quadrature
errors do not enter this continuum enclosure.
These bounds imply

\[
 F(1.0100502,\eta)<-3\cdot10^{-8},\qquad
 F(1.0100504,\eta)>6\cdot10^{-8},\qquad
 F_A>.44\ \text{on }[.9,1.2].
\]

The unique root in the retained bracket is consequently the SAME root.
In fact `k>49/200` and `d<-249/1000`; its discriminant
`D=b^2-4kd>4*(49/200)*(249/1000)>.49^2`. A positive quadratic with negative
constant has exactly one positive root and one negative root. At its positive
root,

\[
 A=-\frac{2d}{b+\sqrt D},\qquad F_A=2kA+b=\sqrt D>.49.
\]

Thus its inverse Jacobian is below `100/49`. The stable displayed formula
also retains the tiny linear term, rather than setting it to zero.
Finally

\[
 |A_\eta|\le\frac{T\lambda}{.49}
 \left[\frac1{64}+512e^{-26}+\frac{80\lambda^{40}}{64}\right]
 <\frac{T\lambda}{30},
\]

because the bracket is below `49/3000`. `axial_bounds.py` checks this using
rational exponential enclosures (128 Taylor terms for exp(-26)).
The actual composed coefficients satisfy the chain rule
`dc/deta=c_eta+c_A A_eta`, so their full C1 norm is at most
`(1+Tlambda/30)lambda^20<=2lambda^20`, as required by the retained budget.

## 6. Fixed jets and remaining obligations

All target eta derivatives are smooth on the compact interval, one-sided
at the endpoints. For R, the eta-dependent angular edits have compact support;
the later infinite tail has a common positive integrable radial majorant at
each fixed derivative order. Each derivative can therefore pass under the
complete integral. No common complex heat strip is required.
With normalized jets `u_n=A^(n)/n!`, `b_n=b^(n)/n!`, `d_n=d^(n)/n!`,
the SAME positive root has the exact recursion for n>=1

\[
 (2kA+b)u_n=-k\sum_{j=1}^{n-1}u_ju_{n-j}
                  -\sum_{j=1}^{n}b_ju_{n-j}-d_n.
\]

The reciprocal is bounded by `100/49`. This gives every fixed jet with
finite constants; it is not a simultaneous smallness assertion in all orders.
The axial coefficients then obey
`c_j^(n)=a_j A^(n)+v_j [eta(1+eta^2)]^(n)`.

The next gate is the ACTUAL annular and compact-repair targets, inverses and
gaps, retaining derivatives of their eta-dependent bilinear terms. Complete
the current core's global backward-pressure identification with all matching
roots and the fixed gauge, then the remaining finite-background/pulse lemmas.
Closing these reference roots alone does not finish those obligations.

The 36 tests include independent exact polynomial moment/energy identities,
omitted-cross-term and parity failures, the composed-derivative chain term,
three finite actual-bump matrix/stable-solver checks, scalar discriminant and
implicit derivative controls, rational continuum inequalities and frozen-record
scope/provenance mutation controls. Polynomial and moderate-parameter fixtures
are diagnostics, not evaluations of the selected Md64 profile. Frozen evidence
binds 18 sources and preserves 1051 inherited source checks; it is never
overwritten. Passing tests is not a formal proof of the remaining construction.

From this directory, with the listed requirements installed:

```sh
python -m unittest discover -s . -p 'test_*.py' -v
python audit.py --verify-record evidence.json --output /a/fresh/axial_replay.json
```
