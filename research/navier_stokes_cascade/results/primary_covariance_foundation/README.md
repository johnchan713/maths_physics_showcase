# Primary covariance roots: explicit cone gaps, Haar inverse and edge weights

Parent: `b17351e44a98b0ca8ef02a936ea7a0d0e94a4fc4`.
Status: `explicit-primary-covariance-root-gap-inverse-and-weighted-jets-from-reviewed-leading-frame-inputs`.

The sourced inverse is already explicit. This checkpoint supplies the NEXT
finite positive-root/gap estimate for the ORIGINAL primary covariance
matrix. It gives a sharp entry-error positivity threshold, a safe inverse
bound, explicit target-freezing and pulse-length budgets, and all fixed
root and reciprocal-root jets with the actual Haar area and shell weights.
No primary root, frequency, support, width `C^-20`, pressure gauge or
historical evidence is changed.

The result is relative to the retained leading/core, background and
moving-frame estimates. It does not independently prove every frame floor,
every generated nonlinear source/mean target, a full residual theorem,
or blow-up. Numerical homogeneous pulses below are manufactured controls,
not evaluations of the enormous actual selected background. The continuum
argument is the set of inequalities below, not a finite corner scan.

The original defining maps are (7.21)--(7.32) of the supplied manuscript,
SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its final theorem is not assumed. The current same-root global moment/cone
join and coefficient/reserved-pressure corrections are retained through
the [sourced-inverse parent](../sourced_inverse_foundation/README.md).

## 1. Construct the homogeneous columns without assuming their roots

On the SAME accepted chart freeze the representatives, carriers, signs,
radius and covers. Let L=Ls, r=r0, and let lambda(v)>0 denote the pulse
growth rate (not the small matching parameter also called lambda).
The reviewed frame equation is

\[
 z'=[\operatorname{diag}(\lambda,-\lambda)+E-dI]z,
 \quad z(0)=(P(0),0),\quad P'/P=\lambda-d_{ref},
 \quad |E_{ij}|\le C_E/L,\quad|d-d_{ref}|\le C_d/L.       \tag{1}
\]

The parent's finite-interval Picard argument applies to these nonzero
initial data as well as zero data. Hence existence and uniqueness do not
require a separate finite homogeneous-ODE lemma or source-size condition.
For mu=min lambda>0, choose L>=8C_E/mu. The Riccati barrier for
`r_z=zminus/zplus` has radius rho=C_E/(mu L)<=1/8 and outward derivative
at either boundary at most `-47 C_E/(64 L)<0`. Integration of the plus
equation then gives

\[
 |r_z|\le\rho,\quad |\log(z_+/P)|\le C_d+C_E(1+\rho),
 \quad x_-P\le x=z_++z_-\le x_+P,
 \quad x_\pm=(1\pm\rho)e^{\pm[C_d+C_E(1+\rho)]}.         \tag{2}
\]

The positive scalar exponential cannot vanish; this closes the barrier's
continuation argument. Fixed differentiated ODEs and the m=1 propagator
bound give polynomial S losses for every coefficient jet relative to P.
There is NO positive lower bound claimed for the oscillating cutoff curl.

The reviewed Gaussian bounds have positive cG,CG:
`exp[-CG(v-L/2)^2/(2L)]<=P<=exp[-cG(v-L/2)^2/(2L)]`.
Let J be the actual Haar rectangle Jacobian, with positive fixed bounds
Jmin,Jmax; `ci=2r/L`. With the original angular average 1/2, transverse
cutoff chi and temporal cutoff psi,

\[
 h_\sigma={J\over2}\left(\int\chi^2d\xi\right)c_i
             \int_0^L\psi^2x_\sigma^2dv,
 \quad s_0={r^2\over\sqrt L},\quad h_\sigma=s_0j_\sigma.  \tag{3}
\]

For L>=25 the central interval `|v-L/2|<=sqrt(L)` lies where psi=1.
Using `(r/2)e^-2/3<=integral chi^2<=r` gives explicit floors

\[
 j_*=J_{min}e^{-2/3}x_-^2e^{-C_G}\le j_\sigma\le
 J_{max}x_+^2\sqrt{\pi/c_G}=j^*.                         \tag{4}
\]

The Gaussian first absolute moment divided by L is bounded, while its mass
is bounded below by a constant times sqrt(L). Thus the normalized column
direction error is at most `Cdir/sqrt(L)`. Differentiating the integrals
uses the same Gaussian majorants and polynomial jet bounds; j_sigma and
the normalized column matrix have all fixed polynomial S jet bounds.
No inverse of a root is used in deriving these column estimates.

## 2. Sharp normalized cone threshold and a safe actual inverse

Let N,K be the frozen orthonormal directions, A=|c0|sqrt(1+u^2)>0 and
u>0 the retained tilt. In the basis (-N,K), the actual two columns are
`h_sigma[(A,-sigma*u)+error]`. Normalize the two rows by A and u and
let every normalized entry error have magnitude at most eta<1:

\[
 C=\begin{pmatrix}1+a&1+b\\-1+c&1+d\end{pmatrix},
 \quad |a|,|b|,|c|,|d|\le\eta,
 \quad H=O\operatorname{diag}(A,u)C\operatorname{diag}(h_+,h_-).
                                                               \tag{5}
\]

Use the ORIGINAL target `T=(Q/q)^(A_similarity+1/2)T0`, not an altered
stress. Put V=(-N.T)/A>0, W=(K.T)/u, and suppose
`|W|<=(1-gamma)V`, `0<gamma<=1`. For q_sigma=h_sigma*y_sigma,
`Cq=(V,W)`. Exact determinant and Cramer formulas give

\[
 2(1-\eta)^2\le\det C\le2(1+\eta)^2,\quad
 \|C^{-1}\|_1,\|C^{-1}\|_\infty,\|C^{-1}\|_2
 \le{1+\eta\over(1-\eta)^2},                             \tag{6}
\]
\[
\min\!\left\{{n\over2(1+\eta)^2},{n\over2(1-\eta)^2}\right\}V
 \le q_\sigma\le{(1+\eta)(2-\gamma)\over2(1-\eta)^2}V. \tag{7}
\]

Here `n=gamma-(2-gamma)*eta`. When n>=0 the lower bound uses the
larger determinant; when n<0 it must use the smaller determinant.
This sign distinction was corrected during pre-publication review and is
covered by a negative-numerator regression, not inferred from the positive
corner controls.

Thus the uniform positivity condition is precisely
`eta<gamma/(2-gamma)`. For 0<gamma<1 it is SHARP in this entrywise error class:
take `C=[[1,1+eta],[-1,1-eta]]`, `W=(1-gamma)V`. Its determinant is 2,
but `qplus=[gamma-(2-gamma)eta]V/2` vanishes at the threshold and becomes
negative above it. Matrix invertibility alone does not give positive roots.
For gamma=1, the threshold eta=1 is the boundary of the assumed error class.

The convenient safe budget eta<=gamma/4 gives

\[
 \det C\ge9/8>1,\quad\|C^{-1}\|_2\le20/9<3,
 \quad {\gamma V\over8}\le q_\sigma\le3V,\quad
 {\gamma V\over8s_0j^*}\le y_\sigma\le{3V\over s_0j_*},
 \quad\|H^{-1}\|_2\le{3\over s_0j_*\min(A,u)}.          \tag{8}
\]

The 2-norm follows from the displayed 1/inf bounds; O is orthogonal.
The physical Haar inverse therefore costs `sqrt(L)/r^2`, equivalently
`r^-3/2 S^1/2` with the retained L comparable to rS. This is not a
radius-uniform bound. The square roots of (8) are the unique positive
primary roots on the OPEN active shell, not at its zero-stress endpoints.

## 3. Join the actual global cone, freezing and the existing length budget

The current global leading cone supplies an unfrozen margin gamma0>0.
Suppose freezing changes V,W by at most tau*V0, where V0>0 is the unfrozen
first coordinate. Then

\[
 V\ge(1-\tau)V_0,\quad |W|\le(1-\gamma_0+\tau)V_0,
 \quad\gamma_{frozen}\ge1-{1-\gamma_0+\tau\over1-\tau}. \tag{9}
\]

The exact half-margin condition is tau<=gamma0/(4-gamma0).
The padded tau<=gamma0/8 suffices; use gamma=gamma0/2 in (8).
Smooth limiting stress DIRECTION, its compact nonzero component floors,
and the O(S^-3) box diameter give tau<=Cfreeze S^-3. The possibly unbounded
weighted stress magnitude cancels in this RELATIVE comparison.

If A>=Amin, the column bound following (4) gives
`eta<=Cdir*max(Amin^-1,u^-1)/sqrt(L)`. It suffices that

\[
 S^3\ge8C_{freeze}/\gamma_0,\qquad
 L\ge[4C_{dir}\max(A_{min}^{-1},u^{-1})/\gamma]^2.       \tag{10}
\]

These are finite value/gap conditions, independent of derivative order.
The retained padded length from geometric pulse replay,
`[Cdir*(Amin^-1+u^-1)*(4+gamma)/gamma]^2`, dominates the new L condition.
The freezing condition was already among the preliminary leading-value
acceptance conditions; here it is made explicit, not numerically evaluated.
No actual first integer band is materialized, and no frequency, width,
support, root branch or derivative-dependent starting band is introduced.

## 4. Weighted jets: do not assume a bounded normalized outer target

Keep `zeta=exp(-t1^2/ya^2-4/yb^2)`. The actual weighted first-stress audit
has outer asymptotic

\[
 T_0=e^{-4/yb^2}yb^{-3}(b_0,yb^6b_z),\quad b_0>0.       \tag{11}
\]

Consequently T0/zeta is generally UNBOUNDED, and a primary root divided
by sqrt(zeta) can grow like yb^-3/2. Smooth zero extension of the full
root does not imply boundedness or smoothness of this quotient. The
historical weighted estimates allow polynomial inverse-distance losses;
they are preserved, not replaced by a bounded normalized-target premise.

Write V=zeta*Vhat, where Vhat>=vstar>0 and every fixed jet has polynomial
S and inverse-edge-distance bounds, not necessarily a bounded edge value.
Let `b_sigma=s0*y_sigma/zeta`. Equations (8) give

\[
 b_\sigma\ge\gamma v_*/(8j^*)>0,\quad
 b_\sigma\le3\widehat V/j_*,\qquad
 a_\sigma=\sqrt{y_\sigma}
 ={L^{1/4}\over r}\sqrt\zeta\sqrt{b_\sigma}.            \tag{12}
\]

All normalized H jets are polynomial S; the usual differentiated inverse
identity gives all fixed weighted jets of b. Its positive lower bound
controls the negative powers in the square-root and reciprocal-root
derivatives. Thus (12) has every fixed jet bounded by sqrt(zeta) times
finite polynomial S/edge losses, with retained scale `r^-3/4 S^1/4`.
The EXACT half-exponential constants are t1^2/2 and 2. Every fixed
polynomial edge loss is absorbed by them, giving C-infinity zero extension.
This argument does not assert b is smooth on the closed shell.

For a real, auxiliary-independent subsequent signed target Sigma with the
stated full zeta-weighted jets, use the original map `d=H^-1(Sigma/epsilon)`,
`delta a_sigma=d_sigma/(2a_sigma)`, the input's FULL zeta weight and the
reciprocal half weight leave sqrt(zeta). No positivity of an updated square
or smallness of Sigma is required. Exactly the selected principal cross
covariance is Sigma; the signed square is nonzero and must remain in the
COMPLETE later covariance tensor. Other covariance entries, axial velocity,
leading radial flow and Pi1/Pi2 pressure tails are not erased.

## 5. Checks and the precise stopping point

Fraction-only corner controls check (6)--(8), including the invertible
nonpositive boundary witness. Independent Taylor-polynomial Cramer solves,
positive square-root recurrences and reciprocal-root signed updates through
order six retain all matrix and denominator derivatives. New homogeneous
noncommuting numerical models integrate both sign columns at two pulse
lengths, checking radial barriers, actual Haar masses and direction budgets.
These are manufactured controls, not continuum certificates or the actual
gigantic selected background. The analytical inequalities, not these tests,
justify the conditional estimates.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
Use fresh outputs; historical evidence is never overwritten.
The new helper modules use distinct covariance-prefixed names: an initial
replay exposed a collision with the historical root-bounds module. No
historical file was changed; a fresh-process import regression covers it.

Next independently close the actual frame/background input floors and
generated COMPLETE covariance/mean targets against the current same-root
pressure/support construction, then every viscous/cutoff/interaction residual
and the single-field realization. The explicit finite primary root solve is
now supplied from those stated inputs; retained exact-core/frame hypotheses
and remaining nonlinear obligations still prevent an unconditional pulse,
full residual or blow-up theorem.
