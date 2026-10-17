# Current coefficient jets, supports and reserved mean foundation

Parent: `135c9da9a41f070ce0d583ba2a70d9a8fe21d490`.
Status: `current-linear-coefficient-jets-reserved-velocity-and-exceptional-pressure-supports-joined-relative-to-leading-inputs`.

This checkpoint supplies an explicit all-fixed-angular-derivative estimate for
the **actual positive-order linear moment solve**, sharpens its normalized
inverse bound from 1000 to 800, and joins the current global moment identities
to the common coefficient supports and reserved mean **velocity**. It also
derives the general-order conservative stress totals with the previous-order
axial viscosity present. The argument is relative to the retained finalized
leading/core input and its regular collar equations. It does not establish
unconditional all-scale physical pulses, complete covariance means, a full
residual estimate, or blow-up. No historical evidence is edited.

The parameters, exact leading roots, frequency, width `C^-20`, labels, gauges
and supports remain unchanged. This advances the first remaining obligation
of the [current pressure/moment join](../global_pressure_moment_foundation/README.md),
not an earlier completed checkpoint. The historical
[support audit](../coefficient_support_audit/README.md) supplied a qualitative
induction with an inverse bound 1000; the new work exposes the actual
normalization jet costs and independently checks the full general-order
conservative equations with exact rational radial integration.

## 1. The actual linear map, not another leading quadratic repair

On the unedited I3 patch let `X=rho exp(y)`, `rho=rho3`, `0<y<5`.
The current leading inputs, after the same-root global join, are

\[
 U_0=0,\quad E_0=e(\eta)e^{-(1/2+\lambda)y},\quad
 e(\eta)=\frac{e_*}{1+\eta^2}>0,\quad M_0=\mu\eta,\quad V_0=-\mu.
 \tag{1}
\]

Here `e_*>0` is the unchanged eta-independent power amplitude; the leading
radial velocity is not zero. At order `n>=1`, all lower orders have already
been globally finalized. Define the five ordinary total discrepancies by

\[
 m_1=\int U_n,\quad m_2=\int\sqrt{2X}E_n,\quad
 m_3=\int(\Pi_n)_X,\quad
 m_4=\int\sqrt{2X}\sum_{i+j=n}U_iE_j,\quad
 m_5=\int\left(\sum_{i+j=n}U_iU_j-X(\Pi_n)_X\right).
 \tag{2}
\]

All integrals include the complete lower pressure source through its true
support, not just the I3 patch. At the current order the only terms involving
`n` are `(0,n)` and `(n,0)`; a current-order square belongs to `2n>n`.
For additive unit-mass bumps `delta U_n=sum c_j beta_j` and
`delta E_n=sum c_j beta_j`, the exact changes in (2) are

\[
 \int\delta U_n,\quad\int\sqrt{2X}\delta E_n,\quad
 \int E_0\delta E_n/X,\quad
 \int\sqrt{2X}E_0\delta U_n,\quad-\int E_0\delta E_n.
 \tag{3}
\]

The two U centers are 1,2; the three E centers are .5,1.5,2.5; all widths
are .1. Divide the ordinary rows by
`(rho,sqrt(2)rho^(3/2),e,sqrt(2)rho^(3/2)e,rho e)` to obtain
`(M,I,Cp,J,S)`, then use `(M,(M-J)/lambda,I,-S,Cp)`.
The matrix is exactly the **linear part B** of the reviewed compact repair:

\[
 B_U:\ (e^y,e^y(1-e^{-\lambda y})/\lambda),\qquad
 B_E:\ (e^{3y/2},e^{(1/2-\lambda)y},e^{(-1/2-\lambda)y}).
 \tag{4}
\]

The same positive unit-mass determinant bounds, recomputed with rational
exponential enclosures, give `||B_U^-1||<2`, `||B_E^-1||<787.562` and
therefore `||B^-1||<800`, for the retained `0<lambda<=.01`.
Unlike the leading repair, the coefficient solve is exactly

\[
                       c_n=-B^{-1}z_n.                    \tag{5}
\]

It needs no small-target condition. Reusing `Bc+Q(c,c)` here introduces
spurious terms from later expansion orders. The diagnostic deliberately
uses non-small targets and detects that mistake.

## 2. Explicit angular jets, with amplitude and slope costs

Use the product algebra norm `|f|_k=sum_(j<=k) sup|partial_eta^j f|/j!`.
Writing `K=1/(sqrt(2)rho^(3/2))`, the exact reciprocal-amplitude norm is

\[
 W_k=|e^{-1}|_k=\frac1{e_*}
 \begin{cases}2&k=0,\\4&k=1,\\5&k\ge2.\end{cases}       \tag{6}
\]

All higher reciprocal jets vanish because `e^-1=(1+eta^2)/e_*`; the
first and second ones do not. If `mu_i=|m_i|_k`, (5) gives the concrete bound

\[
 |c_n|_k\le800\max\left\{
 \frac{\mu_1}{\rho},\quad
 \frac{\mu_1/\rho+KW_k\mu_4}{\lambda},\quad
 K\mu_2,\quad\frac{W_k\mu_5}{\rho},\quad W_k\mu_3\right\}.
 \tag{7}
\]

This is a bound for every fixed order, not a uniform bound as `n`, `k`,
`lambda^-1` or `e_*^-1` grow. Since B is eta-independent, the same inverse
acts on every differentiated normalized target. No fresh root or high-order
smallness hypothesis is used. Freezing `e(eta)` while differentiating the
ordinary moments is incorrect; independent re-solved second jets detect it.

For each fixed radial order `r`, disjoint supports give
`|D_y^r delta U_n|_k+|D_y^r delta E_n|_k<=2 B_r |c_n|_k`,
where `B_r=max_j sup|partial_y^r beta_j|<infinity`.
The conversion `partial_X^r=X^-r sum_j s(r,j)D_y^j` retains its fixed
radius and bump-width costs. Thus all fixed mixed correction seminorms are
finite on the actual patch. This is not a physical all-scale pulse estimate.

## 3. Finite induction and the signed totals needed for support

For each finite n the inner six-variable Volterra system has finite
coefficient/source bounds on the **same radial collar**. Its analytic eta
strip can shrink with n. Every derivative block maps the first four
components into the last two and kills the last two, including at different
radii with an intervening diagonal kernel. Hence at most `ceil(k/2)`
derivative blocks survive a Picard word of length k. With finite coefficient
bound C, strip loss Delta and `z=2Ca`, the sum is bounded by
`(z+beta)exp(beta)`, `beta=e z^2/Delta`. This constructs each finite inner
coefficient from finalized lower sources, not a uniform analytic series in n.
The regular axis primitives and eta-independent inner cutoff preserve the
required collar data. Equation (5) then completes that coefficient globally.
There is no future-support feedback into its zero axis pressure datum.

To justify backward **stress**, local support alone is insufficient. Put
`A=1/2+h`, `D=1/2-h`, `L=1-2h eta^2`, `d=1-eta^2`, `lambda_n=2nh` and

\[
 a_n=1-h+\lambda_n,\quad b_n=3/2-2A+\lambda_n,\quad
 g_n=1-A+\lambda_n,\quad s_n=1-2A+\lambda_n.
\]

After removing the physical similarity powers, the conservative residual
totals are exactly

\[
 \int R^2\mathcal R_{\theta,n}\,dR
 =\frac{-a_n I_n+D\eta(I_n)_\eta+2b_n\eta J_n+d(J_n)_\eta}{L}
 -\nu Z_{a_{n-1}-D}Z_{a_{n-1}}I_{n-1},                 \tag{8}
\]
\[
 \int R\mathcal R_{z,n}\,dR
 =\frac{-g_n M_n+D\eta(M_n)_\eta+2s_n\eta S_n+d(S_n)_\eta}{L}
 -\nu Z_{g_{n-1}-D}Z_{g_{n-1}}M_{n-1}.                 \tag{9}
\]

On a scalar moment, `Z_b f=(2b eta f+d f_eta)/L`.
For `n=1`, the lower angular moment in (8) is the **renormalized** leading
moment, not divergent raw I. The stationary power tail has zero physical
axial derivative, and the fixed positive-h heat majorants justify the
subtraction and differentiation. It is zero by the current global join.
The lower axial total is zero by leading M. At later orders both previous
moments vanish by induction. All current moments vanish by (5).
Therefore (8)--(9) vanish identically and the negative forward stress
integrals equal the positive backward ones with no unknown constants.
The axial viscosity belongs to `n-1`, not n; dropping it gives the wrong
identity for imperfect or uncorrected lower moments.

`conservative_checks.py` uses a new Fraction-only angular-jet/Laurent-radial
calculus, not the old floating-point residual routine. It directly integrates
48 nonzero-component polynomial residual identities at n=1,2,4, including
eta=+-1. The tests also exercise n=7. Lower moments are deliberately nonzero
so omission of viscosity is detectable. These algebraic fixtures do not
solve the actual coefficient equations or prove continuum existence.

## 4. Actual supports and what remains on the mean patch

All I3 bumps end at `y=2.55`, before `Xvel=rho3 exp(3)`. The reserved I4
has `6<log(X/rho3)<11`. Outside Xvel, `F_n=U_n=0`; the exact first total
makes `M_n=0` and hence `V_n=0`, including all fixed eta jets.
The coefficient pressure is reconstructed, not cut off:

\[
 (\Pi_n)_X=\sum_{i+j=n}F_iF_j-\Omega_{n-1}/(2X).
 \tag{10}
\]

After Xvel its only possible sources are the leading radial source in
order one and the shifted leading radial axial viscosity in order two.
The third zero total in (2) removes their final pressure constants.
Together with (8)--(9), the exact right edges are:

| Order | Velocity | Pressure | Stress |
|---|---|---|---|
| 1 | Xvel | Xv | Xb |
| 2 | Xvel | Xv | Xv |
| n>=3 | Xvel | Xvel | Xvel |

The inner equations give `T_n=0` through `Xinner=Xa exp(tc/16)`.
Thus higher stress supports share a compact interior interval and T1 alone
can reach the outer heat edge, consistently with the weighted-stress audit.

On I4, the leading constant `V0=-mu` gives `Omega0=-mu^2/(2X)` and
`Z_-D Z_0 V0=0`. Thus

\[
 \Pi_1=K_1(\eta)-\mu^2/(4X),\qquad
 \Pi_2=K_2(\eta),\qquad\Pi_n=0\quad(n\ge3).             \tag{11}
\]

K1 and K2 contain complete later-source memory and are not assigned zero.
In a manufactured example `M0=eta B(X)`, exact rational differentiation at
eta=0 gives `Z_-D Z_0 V0=6X B'`; a later smooth drop from `B=mu` to zero
leaves backward `Pi2=3nu mu`, explicitly nonzero on the plateau. That
constant is a control, not the value of actual K2.

Cutting the streamfunction by `chi(c_n q)` adds
`-2eta(c_nq)chi'(c_nq)M_n/L` to the normalized radial flux. It vanishes on
I4 because **M_n and all its jets** vanish there. Hence any locally finite
potential-cutoff sum is exactly leading in velocity on the reserved mean
patch, but its normalized pressure still includes

\[
 q^{2A}p_B=\Pi_0+\chi(c_1q)q^{2h}[K_1-\mu^2/(4X)]
                       +\chi(c_2q)q^{4h}K_2.             \tag{12}
\]

The remaining mean/covariance work must retain (12), leading radial velocity
and full-wave axial velocity. Neither zero positive-order velocity nor a
zero order-two local pressure derivative removes those terms.

## 5. Reproduction and remaining obligation

Run `python -m unittest -v test_audit.py` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
The CLI refuses to overwrite evidence. New tests include exact rational
jets/integrals, independent ordinary-integrand repairs of non-small inputs,
re-solved second eta jets, nonlinear-repair and frozen-amplitude failures,
pressure memory controls, immutable sources and scientific scope guards.
The new diagnostic restores the complete Python import search path after
loading legacy helpers; an isolated-process regression covers the module
selection failure discovered during this checkpoint's first counted run.
Quadrature is labelled diagnostic; the continuum inverse is the rational
positive-bump determinant argument. Historical records remain unchanged.

Next independently discharge the remaining actual-pulse defining/existence,
inverse, strict-gap and uniform value/derivative/support estimates at every
physical scale. Then join complete covariance means, including (12) and
axial velocity, and the full viscosity/cutoff/interaction residual. Retained
leading-core dependencies and these obligations still prevent any
unconditional pulse, residual or blow-up conclusion.
