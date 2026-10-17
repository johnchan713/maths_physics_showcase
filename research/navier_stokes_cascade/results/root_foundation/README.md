# Reference angular, heat and waiting roots

This checkpoint starts at `bf4794ad6e0dc309acab0cc6e84ad74a8ea558a1`.
It supplies existence, uniqueness in the retained small balls, inverse bounds
and strict gaps for **two reference correction maps**, and the scalar terminal
waiting equation. It does not discharge every construction lemma. In particular,
the axial two-moment/strength solve, annular and compact-repair maps, complete
core backward-pressure identification and later finite-pulse lemmas remain open.
No unconditional residual improvement or blow-up follows.

The defining formulas are the existing
[`pulse_moment_audit/pulse.py`](../pulse_moment_audit/pulse.py),
[`outer_pressure_pilot/schedule.py`](../outer_pressure_pilot/schedule.py), and
[`heat_exterior_audit/construction.py`](../heat_exterior_audit/construction.py).
The previous moment and heat checkpoints remain unchanged. This is an
independent foundation for their retained root bounds, not a new root choice.
The accepted parameters, pressure gauge, labels, supports, unit viscosity and
width `C_setup^-20` are unchanged.

## 1. Exact bump and norm

Write

\[
 T=e^{64}+10\ge64,\quad 0<\lambda=e^{-4T}\le e^{-256},\quad
 h=\lambda^2,\quad X_*\ge100,\quad w=3/10,\quad c_o=1/1000.
\]

The step is extended by 0 and 1 outside `(0,1)` and inside is

\[
 \sigma(x)=[1+\exp(x^{-2}-(1-x)^{-2})]^{-1}.
\]

It has nonnegative derivative of integral one. We use the deliberately coarse
continuum bounds `sup sigma' <= 64`, `sup |sigma''| <= 20000`.
Indeed on `[1/4,3/4]`, put `g=2/x^3+2/(1-x)^3`; then
`sigma'=sigma(1-sigma)g`, `g<=256`, `|g'|<=3072`. Thus the first bound holds
and `|sigma''| <= (256^2+3072)/4 < 20000`.
For `0<x<=1/4`, `sigma(1-sigma)<=exp(4-x^-2)`; derivative majorants are

\[
 (2x^{-3}+16)e^{4-x^{-2}},\qquad
 (4x^{-6}+6x^{-4}+64x^{-3}+352)e^{4-x^{-2}}.
\]

Each summand increases up to `1/4`, where these sums are below one.
Reflection handles the other endpoint; flatness supplies the joins.
For `b_j(y)=sigma'((y-c_j)/w+1/2)` the mass is exactly `w`, the support
has width `w`, and `b_j^2<=64 b_j`. All bumps within each map are disjoint.

The coefficient space is the Banach algebra with norm

\[
 \|c\|_{C^1}=\max_j\{\sup_{[-1,1]}|c_j|+
                         \sup_{[-1,1]}|\partial_\eta c_j|\}.
\]

Endpoint derivatives are one-sided. Both maps below have constant matrices and
bilinear maps in this space. The same assertion is NOT being made about the
eta-dependent annular map. A values-only estimate would not suffice here.

## 2. Actual angular matrix and quadratic term

Use local coordinates with first bump center zero (global center `Tu-3`),
second center two, `Tu=30 log(1/lambda)`, and `beta=1-lambda`.
For relative edit `epsilon=sum c_j b_j`, the exact map is

\[
 B_Ac+Q_A(c,c)=z_A,\qquad
 B_A=\begin{pmatrix}
  \int e^{\beta y}b_0&\int e^{\beta y}b_1\\
  2\int e^{-(1+2\lambda)y}b_0&2\int e^{-(1+2\lambda)y}b_1
 \end{pmatrix},
\]
\[
 Q_A(c,d)=\left(0,\int e^{-(1+2\lambda)y}
                       \sum_j c_jd_j b_j^2\,dy\right),
 \qquad z_A=(-\delta I\,e^{-\beta(Tu-3)},0).
\]

The pressure row has linear factor two and quadratic factor one; it is not
the heat normalization's half-factor. It sets the incoming pressure-memory
change exactly to zero, not to a sampled or truncated approximation.

Let `i,p` be the two first-column masses, and
`n_i=exp(2 beta)`, `n_p=exp(-2(1+2lambda))`. Then
`B_A=diag(i,p) [[1,n_i],[1,n_p]]`. Positive mass bounds give, uniformly on
the larger CLOSED interval `0<=lambda<=1/100`, inverse row-sum bounds

\[
 (\|B_A^{-1}\|\text{ row bounds})
 \le(2.825741,\ 1.108826),\qquad \|B_A^{-1}\|<3.
\]

For clarity the bounding numerators are `p_hi n_p_hi+i_hi n_i_hi` and
`p_hi+i_hi`, divided by `i_lo p_lo(n_i_lo-n_p_hi)`.
Also

\[
 \|Q_A\|\le64w e^{0.153}(1+e^{-2})<25.402264<32.
\]

All mass and exponential bounds used here are rational enclosures, not
numerical quadrature; their construction is given in section 5.

## 3. Actual heat matrix and quadratic terms

Normalize the pressure, minus-energy and angular rows by
`e_*^2 f^2`, `X_* e_*^2 f^2`, `sqrt(2) X_*^(3/2) e_* f`, respectively,
where `f=(1+eta^2)^-1`. The three centers are `1/2,5/2,9/2` and

\[
 a_0=-1-2\lambda,\quad a_1=-2\lambda,\quad a_2=1-\lambda,
 \qquad (B_H)_{kj}=\int e^{a_k y}b_j(y)\,dy.
\]

For rows zero and one,
`(Q_H(c,d))_k=(1/2) integral exp(a_k y) sum c_j d_j b_j^2 dy`;
row two is zero. The minus-energy sign in the row definition is essential.
The half-factor follows from the relative-square expansion `epsilon+epsilon^2/2`.

Writing `m_k=integral exp(a_k y)b_0(y)dy`, `n_k=exp(2a_k)`, gives
`B_H=diag(m_k) V`, `V_kj=n_k^j`. The nodes stay separated:

\[
 n_0\in[e^{-2.04},e^{-2}],\quad n_1\in[e^{-0.04},1],\quad
 n_2\in[e^{1.98},e^2].
\]

Column `k` of the inverse consists of the coefficients of
`prod_(i!=k)(x-n_i)/(m_k prod_(i!=k)(n_k-n_i))`.
Bounding their absolute values with separated interval endpoints gives

\[
 (\|B_H^{-1}\|\text{ row bounds})
 \le(8.809572,\ 14.240873,\ 1.811035),\qquad \|B_H^{-1}\|<15.
\]

The quadratic norm is at most `(64/2)*3w=144/5<29`, since the two nonzero
rows have nonpositive rates on their nonnegative supports.
This closes the actual positive-bump inverse problem, including the orientation
of the inverse: a transposed Lagrange matrix is not a valid substitute.

## 4. Complete targets, including the infinite heat tail

For the angular target let `r=I/(XH)`. Before interpolation,
`r'=1-(1+l)r`, `r(initial)=5/8`, `l>=-lambda`; hence `0<=r<=1/beta`.
The eta-independent prefix discrepancy `delta0` has absolute value at most
`1/beta`. Put `a=log(2/(1+eta^2))`, so `0<=a<=log 2` and `|a'|<=1`.
The fixed interpolation duration is `Tf=1000`. The existing formula reads

\[
 \delta I=e^{-\beta Tf+a}\delta_0+
 \frac a\beta\int_0^1e^{-\beta Tf(1-z)+a\sigma(1-z)}\sigma'(z)\,dz.
\]

Its value is bounded by `4/beta` and its eta derivative by `6/beta`, so
`||delta I||_C1<=10/beta`. Consequently

\[
 \|z_A\|_{C^1}\le\frac{10}{\beta}\lambda^{30}
       \exp(3\beta+30\lambda\log(1/\lambda))\le700\lambda^{30}.
\]

For the last inequality use `beta>=.99`, `e<11/4`, and monotonicity of
`lambda log(1/lambda)` below `e^-256`, so its multiple by 30 is less than one;
`(1000/99)(11/4)^4<700`. The matrix's larger parameter interval does NOT
justify extending this target bound beyond the stated selected range.

The heat multiplier has its COMPLETE Gamma representation

\[
 H(Z)=\Gamma(1+h)^{-1}\int_0^\infty e^{-v}v^h(1+Zv)^{-h}\,dv,
 \qquad Z=2(1-\eta^2)/X.
\]

With `alpha=chi_heat(1-H)` and `E_new=E_ref(1-alpha)`, positivity and
differentiation of this integral give
`0<=alpha<=3h/X`, `|alpha_eta|<=5h/X`, `||alpha||_C1<=8h/X`.
Thus `||(1-alpha)^2-1||_C1<=18h/X` when `X>=100,h<=.01`.
Let `X_K` be the start of the heat edit and `e_K=E_ref(X_K)`.
The exact terminal reference tail, including its flat end bump, gives

\[
 E_{ref}(X_Kx)\le2e_K x^{-1/2-h}\quad(x\ge1).
\]

Indeed its extra factor is the ratio of
`F_o(y)=1-rho sigma((3-y)/2)`, `rho=h/1000`, at two points, bounded by
`1/(1-rho)<2`. Direct integration of the complete tails yields

\[
 \|\Delta C_p\|_{C^1}\le18he_K^2/X_K,\quad
 \|\Delta S\|_{C^1}\le36he_K^2,\quad
 \|\Delta I\|_{C^1}\le16\sqrt2 e_K\sqrt{X_K}.
\]

The last estimate retains `h integral_1^infty x^(-1-h) dx=1`.
Setting `h=0` before that integral would lose a genuine contribution.
The steep hold of slope `l=-1`, length `4 log(1/h)`, attenuates
`H=E sqrt(X)` by `h^4`. Subsequent release has nonpositive slope, and the
angular correction has already returned to the reference at both endpoints.
Therefore `r_H=e_K sqrt(X_K)/(e_* sqrt(X_*))<=h^4` and
`r_E=e_K/e_*<=h^4`. The bump perturbation does not undo this endpoint
statement; its slope gaps are checked below.

Normalizing requires DERIVATIVES of `f^-1` and `f^-2`:
their full C1 norms are 4 and 12. Thus the normalized target row bounds are

\[
 216h^9/X_*,\quad432h^9/X_*,\quad64h^4/X_*.
\]

All are bounded by `delta_H=64h^4/X_*`. Signs are negative for the pressure
and I losses and positive for the minus-S row. All integrations and identities
are eta FUNCTION identities, not only finite sets of point evaluations.

## 5. Contraction, same accepted roots, jets and gaps

For `Bc+Q(c,c)=z` with bounds `beta,q,delta`, the map
`B^-1(z-Q(c,c))` on radius `r=2 beta delta` is a strict contraction into
that ball if `8 beta^2 q delta<1`. Its Lipschitz constant is at most
`4 beta^2 q delta<1/2`; the image bound is `beta delta+beta q r^2<r`.
Use `(beta,q)=(3,32)` and `(15,29)` to obtain

\[
 \boxed{\|c_A\|_{C^1}\le4200\lambda^{30}},\qquad
 \boxed{\|c_H\|_{C^1}\le1920h^4/X_*}.
\]

These sharpen the previous `7000 lambda^30` and `6400 h^4/X_*` radii.
The previous inverse constants 5 and 50 also give contractions on those
previous balls at the selected parameters. The new balls are contained in
them, so uniqueness identifies the SAME accepted roots. No profile, support,
parameter or cutoff sequence has been replaced.

The pointwise Jacobian `J=B+2Q(c,.)` has inverse norm at most `2 beta`,
namely 6 and 30. Differentiating the SAME roots gives, for normalized fixed
jets `d_k=||partial_eta^k c||/k!`, `z_k=||partial_eta^k z||/k!`,

\[
 d_k\le2\beta\left(z_k+q\sum_{j=1}^{k-1}d_jd_{k-j}\right).
\]

This proves each fixed smooth one-sided jet from the corresponding smooth
target. It is not a simultaneous contraction in all derivative orders, and
does not assert a common complex strip for heat targets at the endpoints.

The relative value and eta derivative of each edit are separately at most
`64r`; its radial derivative is at most `(20000/w)r`. Hence positivity
holds and the log-slope loss is at most `(20000/w)r/(1-64r)<lambda/2`.
The reproducible binary upper exponents (use `lambda<2^-256`) are:

| Quantity | Angular | Heat |
| --- | ---: | ---: |
| contraction condition left side | -7659 | -2032 |
| coefficient radius | -7667 | -2043 |
| relative edit, or each eta derivative | -7661 | -2037 |
| radial log-slope change divided by lambda | -7395 | -1771 |

These are analytic upper bounds, not evaluations of the enormous selected
coordinates or tiny accepted coefficients.

`root_bounds.py` uses rational Taylor enclosures: for nonnegative `x<n+2`,
the sum through order `n` is a lower bound for `e^x`, and the tail is at most
`[x^(n+1)/(n+1)!]/[1-x/(n+2)]`. Negative arguments use reciprocals.
Order 64, rational outward rounding to denominator `10^9`, and final upward
rounding to `10^6` give the displayed inverse bounds. All enclosing arguments
here have absolute value at most 2.04. The separate finite moderate-parameter
quadratures in `root_checks.py` only diagnose transcription errors.

## 6. Positive unique waiting time for the power reference

This scalar solve belongs to the power reference BEFORE adding the heat
attachment and its compensating bumps. Distinguish
`q_w=(1-h)I/(XH)-1` from the similarity-coordinate variable q.
After angular correction, `q_w=(lambda-h)/(1-lambda)>0`, independent of eta.
The exact ODE is

\[
 q_w'+(1+l)q_w=-l-h.
\]

On both release transitions the source is nonnegative. The intervening
`l=-1` hold of length `4 log(1/h)` gives `q_w>=4(1-h)log(1/h)`.
The last transition loses at most a factor `e^-1`, so
`q_before>=4*.99*512/(11/4)=18432/25>1`.
For the exact nonzero terminal flat bump put

\[
 q_p=\frac\rho{1-\rho}\int_0^1
       e^{(1-h)(1+2v)}\sigma'(v)\,dv,\qquad \rho=h/1000.
\]

It is strictly positive for every selected `h>0` and at most
`1331/6399936<1`. Therefore

\[
 \omega=\frac{\log(q_{before}/q_p)}{1-h}>0
\]

is finite, unique and eta-independent. On the terminal segment
`l=-h+F_o'/F_o`, the integrating factor is
`exp((1-h)y) F_o(y)/(1-rho)`. Integrating its source `-F_o'/F_o`
and substituting `y=1+2v` gives exactly the displayed `q_p`, hence `q_w(3)=0`.
The reference power tail has `I=XH_pow/(1-h)`: its renormalized I constant
vanishes. The later Gamma heat attachment changes the pointwise tail; the
heat map restores its COMPLETE renormalized I moment, not this pointwise
power formula or the reference q_w profile. No reference waiting root is
being re-solved on the final heat-modified profile.
The raw total I integral is not being asserted finite. There is no claim
of a uniform scalar-root condition number as `h` tends to zero.

## 7. Pressure memory, evidence and next gate

The angular Cp row is exactly zero. The heat correction cancels the complete
Cp, S and renormalized-I changes of the Gamma heat edit, including all fixed
eta jets. Neither map depends on the axial-velocity strength. The waiting
condition depends on E and I, not that strength, so there is no circular
strength choice in these three roots. Heat targets use the unedited future
tail on disjoint supports; heat correction restores the future total S.

This does NOT independently establish the complete-core identity
`Pi_axis+Cp_infinity=0`: other exact roots and the global joining/repair
formulas are still required. Inherited conclusions keep those dependencies.
The next justified gate is the actual axial two-moment solve and scalar
strength root, with all later energy and nonzero end-bump terms retained.
Then audit the annular/compact-repair targets, inverses and gaps and finish
complete-core global backward-pressure identification. Finite-background and
finite-pulse existence/inverse/gap obligations remain thereafter.

`evidence.json` binds 18 sources and preserves 1033 inherited source checks.
It is immutable once frozen. The 40 tests include rational enclosures,
44 two-sided exact inverse-product components, transposed-inverse rejection,
disjoint relative-square coefficient controls, 12 Catalan jet orders,
factorial rejection, terminal integrating-factor controls, six finite bump
integral diagnostics and scope/provenance mutation controls. Polynomial bump
proxies in algebraic tests are NOT the actual flat bump; the continuum matrix
proof above uses the actual bump's positivity, exact mass and derivative bound.
Passing these checks is not a formal certificate or a continuum proof of
the remaining Navier--Stokes construction.

Reproduce with the requirements in this directory:

```sh
python -m unittest discover -s . -p 'test_*.py' -v
python audit.py --verify-record evidence.json --output /a/fresh/root_replay.json
```

The output path must be fresh; the audit refuses to overwrite evidence.
