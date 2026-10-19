# Same heat-root C3 bound and the remaining global swirl inputs

Status: `same-heat-root-third-jets-and-global-swirl-input-envelope`.
Parent: `ca0e6a3aec5b491bfd0069418261a092d8fdd03e`.

The preceding checkpoint supplied the leading axial and pointwise radial
pressure-source inputs to the first-stress estimate. This checkpoint supplies
the remaining swirl inputs through the later reference, compensation and
heat regions. A new third-angular target estimate identifies the SAME heat
compensation root, rather than assuming that its old C1 bound is a C3 bound.
No root, profile, cutoff, width, frequency, gauge, support or mean correction
is changed. Historical evidence is not overwritten.

All conclusions are analytical implications of the retained exact core,
regional input, matching and positive-tail lemmas. They are not independent
proofs of those premises or a continuum proof from passing software tests.
The actual coefficient-cutoff shape, full same-field cone/scale acceptance,
physical pulse/covariance bounds and complete residual closure remain open.

## 1. Third-angular heat target, including the whole infinite tail

Retain `lambda=exp(-4T)`, `h=lambda^2`, `T>=64`, `Xstar>=100`, `XK>=Xstar`,
`f=(1+eta^2)^-1`. Thus `lambda<2^-256`, `h<2^-512` for the actual family.
Some inequalities below hold on the larger auxiliary box `0<h<=.01`,
`Xstar>=100`; its corner is an arithmetic control, not the selected profile.
The actual heat matrix still uses the retained `lambda<=.01`.

The positive Gamma factor is

\[
\mathcal H(Z)=\Gamma(1+h)^{-1}\int_0^\infty
 e^{-v}v^h(1+Zv)^{-h}dv,\qquad Z=2(1-\eta^2)/X\ge0.
\]

Its derivatives obey `|Hcal^(j)|<=(h)_j(1+h)_j`. For j=1,2,3 and
h<=.01 these are bounded by `(2h,4h,32h)`. The exact angular chain rule
retains the second derivative of Z; for example

`Hcal_etaetaeta=Hcal''' Z_eta^3+3Hcal'' Z_eta Z_etaeta`.

Since `|Z_eta|<=4/X`, `|Z_etaeta|=4/X`, for X>=1 raw derivative caps for
`Hcal-1` through order three are `(4,8,72,2240)h/X`. The factorial C3
norm is therefore below `512h/X`. Endpoint derivatives at eta=+-1 are
one-sided; the Gamma integral is never evaluated at Z<0. Its value there
is one, but its first eta derivative need not vanish.

The actual heat cutoff is `chi(y)=sigma((y-.2)/.3)` with y=log(X/Xtail).
It is eta-independent and between zero and one. Put `alpha=chi(1-Hcal)`.
Then `|alpha|3<=512h/X<1` for X>=XK>=100, and
`|(1-alpha)^2-1|3<=1536h/X`.

On the heat support and its ENTIRE future, the unedited Eref is
eta-independent: interpolation and angular bumps have already ended.
Retain the original positive bound

`Eref(XK*x)<=2eK*x^(-1/2-h)`, x>=1, `eK=Eref(XK)`.

Positive integration, not a uniform bound times an infinite length, gives

\[
|\Delta C_p|_3\le3072he_K^2/X_K,\quad
|\Delta S|_3\le3072he_K^2,\quad
|\Delta I|_3\le1024\sqrt2e_K\sqrt{X_K}.                 \tag{1}
\]

For the last bound the integrand is proportional to
`h*x^(-1-h)` and `h int_1^infinity x^(-1-h)dx=1`. That finite contribution
must not be lost by putting h=0 before integration. All three discrepancies
are complete functions of eta, not endpoint or finite-radius surrogates.

Normalize the three compensation rows by `e_star^2 f^2`,
`Xstar e_star^2 f^2`, `sqrt(2)Xstar^(3/2)e_star f`. The FULL C3 reciprocal
norms are `|f^-1|3=5`, `|f^-2|3=24`, not their C1 values 4 and 12.
The retained attenuation gives

`rH=eK sqrt(XK)/(e_star sqrt(Xstar))<=h^4`,
`rE=rH sqrt(Xstar/XK)`.

Thus the target row norms from (1) are at most

\[
\left(73728h^9/X_*,\ 73728h^9/X_*,\ 5120h^4/X_*\right)
\le\delta_3:=8192h^4/X_*.                              \tag{2}
\]

No later pressure memory is truncated to the compensation support. The
pressure, minus-energy and renormalized-I signs and normalizations are the
same as in [root_foundation](../root_foundation/README.md).

## 2. Bound and identify the same accepted heat root

The exact heat B,Q are eta-independent, with `|B^-1|<15`, `|Q|<29`.
They therefore have the SAME bounds in the factorial C3 product algebra.
The quadratic contraction with (2) gives

\[
|c_H|_3\le30\delta_3=245760h^4/X_*
        <\mathfrak c_H:=2^{18}h^4/X_* .                \tag{3}
\]

Indeed `8*15^2*29*8192 h^4/Xstar<1` even at h=.01,Xstar=100.
The old accepted C1 root has pointwise size at most `6400h^4/Xstar`
(the later sharper radius is 1920 times that scale). Both old and new roots
are inside the common POINTWISE ball of radius mathfrak cH. On this ball
the exact fixed-point map has Lipschitz constant at most
`2*15*29*mathfrak cH<1`. Hence both roots agree at every eta; (3) is a
bound on the unchanged field, not a root reselection.

The new, larger C3 ball is NOT claimed to be contained in the old C1 ball.
The proof uses common pointwise uniqueness instead. Both endpoint values
and one-sided jets are retained. This does not replace the original inverse
or moment-identification premises by numerical small residuals.

## 3. Mixed J3 bounds for every later reference factor

Write `Ahi=(C^4096)^256`, `G=exp((C^4096)^32)`, `H=exp(G^4)`.
The earlier fixed-frequency proof gives `Ahi<H`, `H>2^1024`, and native
`|E|J3<H^70` on the compact modulation/first-repair annulus, with unchanged N.
That annulus ends at `Xrep=rho_p exp(5)=Xstar`. At its endpoint the field
and its fixed jets equal the prescribed reference; all modulation/repair
localizers have ended. The value budget implies `E(Xstar,0)<=2Ahi`.

For X>=Xstar, factor the reference swirl into its positive scalar slope
profile S, angular interpolation gamma, two relative-edit factors, the
original terminal taper and the heat multiplier. Each factor uses the
original supports. They need not be simultaneously active; taking their
product norms is a conservative bound.

The scalar S starts at at most 2Ahi and is nonincreasing after Xstar.
Its logarithmic derivative b has `|b|<=2`, `|b_y|<=64`, `|b_yy|<=20000`:
the later slopes are between -1 and zero, their unit-length transitions
use the same sigma, and the constant pulse, hold and waiting stages have
zero slope derivatives. All stage endpoints, including the scalar waiting
time, are eta-independent by the retained root construction. The exact
identities `S_yy=S(b^2+b_y)`, `S_yyy=S(b^3+3bb_y+b_yy)` yield

`|S|J3<=2Ahi[1+2+68/2+20392/6]<8192Ahi`.

Before interpolation gamma=f. During it,
`gamma=exp[-log 2+theta(y) log(2/(1+eta^2))]`,
`theta=1-sigma(y/1000)`, and afterwards gamma=1/2. Retain
`|log(2/(1+eta^2))|3<8`; the step caps 32,1216,220000 imply
`|theta|J3<2`. The product algebra gives `|log gamma|J3<17`, so
`|gamma|J3<exp(17)<2^34`.

Each width-.3 bump in the heat and angular relative edits is an UNNORMALIZED
sigma-prime translate, with mass .3. The already retained fourth-derivative
step bound 1e100 gives `|beta|J3<1e103<2^412`. The two angular coefficients
have C3 norm at most `2^32 lambda^30`; the three heat coefficients obey (3).
Consequently both relative-edit J3 norms are below one for the actual
`lambda<2^-256`, `h<2^-512`, and each factor `1+edit` has J3 norm below 2.
The heat edit and its compensation remain separate, nonzero moment memories.

Keep the EXACT terminal prescription from `outer_pressure_pilot/schedule.py`:

\[
f_o(y)=1-\rho\sigma((3-y)/2),\quad\rho=h/1000,\qquad
\text{terminal multiplier}=f_o(y)/f_o(0),\quad f_o(0)=1-\rho.
\]

It is one before terminal y=1 and increases on the original two-unit
interval [1,3]. It is not replaced by `1-rho sigma(y)` or by a unit-width
transition. Its J3 norm is bounded by

`[1+rho*(1+32/2+1216/8+220000/48)]/(1-rho)<2`.

It is actually still one on the needed terminal interval through y=1/2;
the bound also covers its complete later continuation.

For the heat multiplier a mixed, not just angular, estimate is needed.
The positive argument `Z=2(1-eta^2)/X` has `|Z|J3<20/X`. Finite Taylor-jet
composition at the current NONNEGATIVE Z, with Gamma derivative caps
`(2h,4h,32h)`, gives, for X>=100,

\[
|\mathcal H-1|_{J_3}
 \le4h/X+2h(20/X)+2h(20/X)^2+(32h/6)(20/X)^3
 <64h/X.                                               \tag{4}
\]

This is finite-order real jet calculus; no complex or negative-Z extension
of the Gamma integral is assumed. The actual heat cutoff has J3 norm below
2e6, including its inverse width .3. Thus the heat multiplier has J3 norm
at most `1+128000000h/X<2` for the actual h. It equals one before its support.

Multiplying S, gamma and the four bounded factors gives the native E bound

\[
|E|_{J_3}<2^{13+34+4}A_{hi}<2^{60}A_{hi}<H^2<H^{70}
\qquad (X\ge X_{rep}).                                \tag{5}
\]

The whole annulus before Xrep and (5) therefore use the same explicit bound
`M_E=H^70` after the axis-side interval. This includes compensation, the
entire pulse swirl, angular correction, waiting, heat and terminal taper.
It is a quantitative mixed bound, not a compactness-only finiteness claim.

## 4. Complete the leading swirl inputs to the first-stress estimate

Let Kaxis be the same activated-collar/common radial-order-0,1,2,
angular-C3 envelope used in the preceding global axial checkpoint.
Only after that axis-side interval use `F=E/sqrt(2X)` and X>=Xa.
The logarithmic J3 norm of `(2X)^-1/2` is less than
`2/sqrt(2Xa)<=2(1+Xa^-1)`. Define

\[
M_F=2(1+X_a^{-1})H^{70},\quad
F_3=\max(K_{axis},M_F),\quad
F_{X,2}=\max(K_{axis},M_F/X_a),
\]
\[
F_{XX,1}=\max(K_{axis},3M_F/X_a^2).                    \tag{6}
\]

They bound `|F|3`, `|F_X|2`, `|F_XX|1` on [0,Xc]. The factor three retains
`F_XX=X^-2(F_yy-F_y)`; no division at zero or radial analyticity at a flat
join is used. Here `Xc=Xtail exp(1/2)` is unchanged.

For the actual swirl exponent `b=-1-h`, its shifted exponent `b-D=-3/2`
also has modulus below two. The complete coordinate-derivative estimate
from [global_axial_source_envelope](../global_axial_source_envelope/README.md)
therefore supplies

\[
\boxed{C_F=64X_c^2 F_{XX,1}+704X_cF_{X,2}+1664F_3
       \ge\sup_{[0,X_c]}|Z_{b-D}Z_bF|_1.}              \tag{7}
\]

In the existing first-stress interface use `f=F3`, `fx=F_X,2`, `CF` from
(7), and the preceding explicit axial inputs `u,ux,v,w,CU`. This now
supplies ALL its primitive leading-field jet/source inputs, relative to
their stated lemmas. `complete_first_interface` composes them with the
same reconstructed first pressure, viscosity terms and compact Cstar formula.
The auxiliary rational composition is not the selected huge profile.

This still does not certify an actual old scale. The first-coefficient jets
retain the actual inner-cutoff norms L1,L2 from the cutoff audit. Supply
valid bounds for those unchanged shape norms and then verify

`2h*(ell_scale-2)>=log2(Cstar(L1,L2)/r_Kc)`

together with the same-tilt, freezing, geometry, inverse and other acceptance
gates. No cutoff shape is silently chosen here. The signed and absolute
pressure-source bounds remain distinct and include the full future Xv;
neither a pressure datum nor memory between reserved supports is reset.

## Evidence, sources and remaining limits

Exact rational controls check the Gamma caps, complete-tail normalization,
root contraction/common pointwise uniqueness, the original terminal taper
and all reference product constants. Symbolic chain rules detect omitted
second-argument and inverse-radius terms. Ten finite positive-integral Gamma
controls include both poles; 45-digit quadrature is rounded to 20 significant
digits in their displayed diagnostic values. These are numerical controls,
not continuum proofs or samples of the selected enormous field.

Sources are the retained manuscript (A.10)-(A.13), (A.34)-(A.46),
`root_foundation`, `annular_exterior_foundation` and the preceding global
checkpoints. [Manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
SHA-256: `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its proposed theorem is not accepted as established. Historical source
bindings and frozen evidence are reproduced under original comparators.

Next: resolve the unchanged coefficient-cutoff shape budget and verify the
full same-field Cstar/tilt/scale gates. Uniform physical pulse bounds,
complete wave-covariance mean corrections and the complete viscous residual
remain separate. No unconditional pulse bound, residual improvement or
blow-up has been proved.

```sh
python3 research/navier_stokes_cascade/results/global_swirl_envelope/test_audit.py -v
python3 research/navier_stokes_cascade/results/global_swirl_envelope/audit.py \
  --verify-record research/navier_stokes_cascade/results/global_swirl_envelope/evidence.json \
  --output /tmp/fresh_global_swirl_envelope.json
```
