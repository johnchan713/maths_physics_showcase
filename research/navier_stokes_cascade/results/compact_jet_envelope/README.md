# Actual compact derivative and denominator envelope

Status: `actual-compact-input-envelope-bounded`.

This checkpoint bounds the actual, unmodulated joined leading profile on the
compact region selected by the [stress proposal](../joined_stress_construction/README.md).
It supplies the second missing input estimate, following the
[actual C2 moment transfer](../angular_c2_transfer/README.md). In particular the
proposed input envelope **A = C^4096** is sufficient for the specific norms and
positive quantities below. Narrow radial cutoffs are differentiated with their
inverse-width factors intact.

This is an analytic estimate with reproducible outward checks of constants and
independent equation paths on manufactured inputs. It relies on the preceding
analytic core and joining arguments; it is not independent peer review or a
formal proof-assistant certificate. The fixed-phase loop, modulation errors,
repair-to-state errors, and proposed frequency are **still unverified for the
actual profile**. No full stress construction, PDE solution, smooth force, or
blow-up is claimed.

Source: the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
especially (4.9), (4.11), (4.15)-(4.16), and Appendix B. Page 28 was visually
checked against the implementation. PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Parent: `028333027590ae2bad71f5be74d699f4cb979c3d`. Earlier audit files and the
goal document remain frozen.

## 1. Norm and exact scope

Write y=log X and use the algebra norms

\[
\|F\|_{\mathcal J_k}=
\sum_{r+s\le k}\frac{\sup|\partial_y^r\partial_\eta^s F|}{r!s!},
\qquad k=1,2.
\]

The product rule gives \(\|FG\|_{\mathcal J_k}\le
\|F\|_{\mathcal J_k}\|G\|_{\mathcal J_k}\). Either derivative maps J2 to
J1 with norm at most 2. Tuples use the maximum component norm. The angular
endpoints use one-sided derivatives. All patch boundaries are independent of
eta. These are logarithmic radial derivatives, **not unweighted X derivatives**.

Use exactly the previously selected parameters and locations:

\[
X_a=4/\Lambda,\quad \delta=C^{-20},\quad
t_c=\delta/(10\sqrt{\log C}),\quad X_-=X_a e^{t_c/2},
\]
\[
X_R=110(CP)^{10},\quad T_w=240T_d,\quad
\rho_p=X_R e^{T_d+T_w-23},\quad X_{\rm rep}=\rho_p e^5.
\]

The input rectangle is \([X_-,X_{\rm rep}]\times[-1,1]\). For integrating
moments, field upper bounds extend back to 0<X<=Xrep; a positive E floor is
claimed only away from the axis. The inner offsets are kept symbolic, never
added to a rounded log Xa. The region ends at the end of the first reserved
repair patch, before the heat-compensation patch, and contains the modulation
support ending at \(X_+=X_R e^{T_d+3}\).

The inherited geometry and new logarithmic checks give

\[
C^{-1}<X_a<X_-<X_{\rm rep}<C^{12},\quad
1<\rho_p<C^{12},\quad P,K,\Lambda,B<C,
\]
\[
\lambda^{-1},h^{-1}<C,\quad C>10^9,\quad
\log(X_{\rm rep}/X_a)<13\log C<C.
\]

In exponent bookkeeping, products add exponents and a sum of at most C terms
bounded by powers of C costs one extra power. This is conservative arithmetic,
not a substitute for deriving the input estimates.

## 2. Fields, including the thin transitions

### Core and reference

The frozen core argument bounds raw mixed derivatives of u and log Phi through
radial order 2 and angular order 6 on 0<=Y=Lambda X<=4.1 by K. Thus converting
to \(D_y=Y\partial_Y\) costs at most 50K at the orders needed here.
The fixed factor has \(\partial_\eta\log\phi_*=\Lambda\zeta_*\), with
the needed derivatives bounded by Lambda K. The exact axis factor is

\[
E=C^{-1}\sqrt{2X}\,\phi_*\Phi,\qquad
U=4\eta+j+\Lambda^{-1}u.
\]

The preceding C2 transfer proves \(\|CE/\sqrt X\|_{C^2_\eta}<e^B\)
and E<1 in the core. The finite derivative bounds and the formulas for
differentiating an exponential therefore give \(\|E\|_{\mathcal J_2}<C^{16}\)
and \(\|U\|_{\mathcal J_2}<C^2\). In particular 50K/Lambda<.01 is checked
outward. No inverse X is introduced by these logarithmic derivatives.

On [Xa,Xi], Xi=110, the directly constructed reference obeys
\(\|p_{1,r}\|_{C^3_\eta},\|n_{s,r}\|_{C^3_\eta},
\|\log(CE_r)\|_{C^3_\eta}<B/100\),
\(\|E_r\|_{C^3_\eta}<e^B/C<1\), and p1,r>=X/10.
Tapering native slopes multiplies them by cutoff values in [0,1]. Hence
\(\|U_r\|_{C^2_\eta}<10\), \(|D_yU_r|<1\), and
\(|l_r|<100K<C\). These assertions do not differentiate the cutoff.
The axis pressure has C1 norm below 15P^2. Splitting the reference Cp integral
at Xa, using its axis factor and log(Xi/Xa)<C, gives
\(\|\Pi_r\|_{C^1_\eta}<C^4\). Also |Wr|,|Hc,r|<C.

Now use the ORIGINAL source equations, not a derivative bound on the already
modified profile:

\[
S_{q,r}=-W_rl_r-h(1-2\eta U_r)-H_{c,r}(\log E_r)_\eta,
\]
\[
S_{n,r}=-W_rD_yU_r-A_{\rm phys}(1-2\eta U_r)U_r-H_{c,r}U_{r\eta}
-d\Pi_{r\eta}+4A_{\rm phys}\eta\Pi_r+\eta E_r^2.
\]

Here \(A_{\rm phys}=.5+h\), \(D=.5-h\), \(d=1-\eta^2\), and
\(L=1-2h\eta^2\ge.98\). The bounds give |Sq,r|<C^4 and |Sn,r|<C^6.
Since \(n_{s,r}=N_{s,r}/L\), equation (4.9) implies

\[
D_yp_{1,r}=XS_{q,r}/L-l_rp_{1,r},\qquad
D_yn_{s,r}=S_{n,r}/L-n_{s,r}.
\]

Both radial derivatives are less than C^8 in absolute value. This avoids
assuming a width-independent second derivative of the reference fields.

### Actual activation and final cutoffs

Put F=log(CE). Before the final cuts,

\[
F_y=\tfrac12-\kappa p_{1,r}/2,\qquad
U_y=-\kappa Xn_{s,r}/2,\qquad
|\kappa_y|\le64\delta^{-1}<C^{21}.
\]

The final axial taper has the same derivative bound; the subsequent a taper
is a convex interpolation from kappa0 p1,r to .8. Pure angular derivatives
through order 2 have the preceding C2 bounds. Differentiating these displayed
radial slopes once gives, conservatively,

\[
|F_y|,|F_{y\eta}|,|U_y|,|U_{y\eta}|<C^2,\qquad
|F_{yy}|<C^{24},\quad |U_{yy}|<C^{25}.
\]

The actual continuation has E<1. Use
\(E_{ij}=E(F_{ij}+F_iF_j)\), for i,j in {y,eta}, rather than
exponentiating the enormous full J2 norm of F. This gives J2 bounds below
C^28 for both fields. Dropping delta^-1 would be invalid even though the
field values and angular errors are small.

### Shape transition, earlier joining annulus, and outer prefix

The shape formula is

\[
F=y_i/10+(1-\sigma(y_i/T_{sh}))\ell_i
+\sigma(y_i/T_{sh})\log f,\quad U=G_i,\quad f=(1+\eta^2)^{-1}.
\]

With Tsh=2000B, sigma'<=64 and |sigma''|<=20000, the radial and mixed
derivatives of F through order 2 are less than 1, while its pure angular
derivatives are bounded by 408B<C. E<1 gives an E J2 bound below C^4.
The following ideal power has E=Pf(X/XR)^.1.

In the earlier joining annulus, the SAME exact coefficients have factorial
C2 norm <=5.72e-14. Its bumps have bounds 64, 250000, 2e8 for value, first,
and second y derivatives. Product differentiation gives, even summing all
five supports, a correction J2 norm no larger than

\[
5(5.72\,10^{-14})(64+250000+10^8)<1.
\]

The relative E factor exceeds 1/2. The axial restoration involving
\(\|g\|_{C^2_\eta}\le3\,10^{-18}\) obeys the same conclusion using
the step bounds. Thus the join has field J2 norms below 20P<C^2.

After the join and through Xrep, E is Pf times a scalar amplitude at most
2, with alpha=Dy log E in [-.51,.1] and |alpha_y|<=39.
The axial field is k(y)eta, where
\(k=4[1-\sigma(\log(1+y)/64)]\) during its taper.
Consequently |k|<=4, |k'|<=4, |k''|<24. It is zero on the subsequent
negative-slope and power stages. Direct differentiation gives
\(\|E\|_{\mathcal J_2}<120P<C^2\), \(\|U\|_{\mathcal J_2}<32<C\).
Together these regional estimates prove the uniform bound

\[
\boxed{\|E\|_{\mathcal J_2},\|U\|_{\mathcal J_2}<C^{64}.}
\]

## 3. Positive denominators

On [Xa,Xi], |log(CE)|<B, hence E>=C^-1 exp(-B)>C^-2. The shape formula
preserves this lower bound; its extra y_i/10 is nonnegative. The ideal
power increases, and the earlier relative join loses at most a factor 2.
On the outer prefix a conservative full-length estimate is

\[
E\ge\tfrac12 P\exp[-.5-T_d/2-.51-.51T_w]
>e^{-123T_d}>C^{-1}.
\]

The last two comparisons are checked in logarithms; no very small E is
materialized. Thus E>=C^-4 throughout the input rectangle. In the actual
activation and continuation a=kappa p1,r>=C^-20 X/10>C^-22.
The final convex interpolation and all later stages preserve a>=C^-22.

Reciprocal differentiation, with the field J2 bound, gives

\[
\|E^{-1}\|_{\mathcal J_2}
\le C^4+C^8\|E\|_{\mathcal J_2}+C^{12}\|E\|_{\mathcal J_2}^2
<C^{144},
\]

and, retaining the sharper first-order estimate,

\[
\|E^{-1}\|_{\mathcal J_1}<C^{73},\quad
\|E^{-1}\|_{\mathcal J_2}<C^{144},\quad
\|H^{-1}\|_{\mathcal J_1}<C^{75},\quad
\|H^{-1}\|_{\mathcal J_2}<C^{146}.
\]

Also \(\|X^{-1}\|_{\mathcal J_2}<C^2\),
\(\|X\|_{\mathcal J_1}<C^{13}\), and
\(\|L^{-1}\|_{\mathcal J_2}<2<C\).
For the future repair, e(eta)=E(rhop,eta); its inverse has the same angular
C2 bound. The five normalization multipliers are
\((\rho_p e)^{-1},(\sqrt2\rho_p^{3/2}e^2)^{-1},
(\sqrt2\rho_p^{3/2}e)^{-1},(\rho_p e^2)^{-1},e^{-2}\).
Since rhop>1, each has C2 norm below C^288. Combining the first two
normalized rows and dividing by the POSITIVE lambda gives an operator bound
below C^290, and normalized actual moment norms below C^546. This controls
the normalization, not the still-missing O(1/N) discrepancy.

## 4. All five original cumulative moments

Integrate the actual fields, not a sampled target or normalized surrogate:

\[
M=\int_0^X U\,dx,\ I=\int_0^X\sqrt{2x}E\,dx,\
J=\int_0^X\sqrt{2x}UE\,dx,\ S=\int_0^X(U^2-E^2/2)\,dx,
\quad C_p=\int_0^X E^2/(2x)\,dx.
\]

Pure angular derivatives commute with these integrals. The logarithmic
radial identities and their first derivatives are

| Moment | Dy moment | Dy^2 moment |
|---|---|---|
| M | XU | X(U+Uy) |
| I | sqrt(2) X^(3/2) E | sqrt(2) X^(3/2) (3E/2+Ey) |
| J | sqrt(2) X^(3/2) UE | sqrt(2) X^(3/2) (3UE/2+Uy E+U Ey) |
| S | X(U^2-E^2/2) | X(U^2-E^2/2+2U Uy-E Ey) |
| Cp | E^2/2 | E Ey |

Differentiate the middle column in eta for the mixed derivative. Xmax<C^12
and the field J2 estimates give safe full J2 bounds C^78,C^86,C^150,C^144
for M,I,J,S respectively. For Cp, its core integral is less than 1 in angular
C2 by the retained sqrt(X) factor. The remaining logarithmic length is less
than C, so its angular C2 bound is at most 1+C^129/2; the displayed radial
identities give a full J2 bound below C^132. Thus every moment has J2 norm
below C^256. A uniform E bound alone would give a divergent axis integral
and is NOT used there.

The frozen positive-mixture axis datum has mass less than 5 and integrands
\(g_\theta=(1+\eta^2)^{-2\theta}\), 0<=theta<=1. Bounds
|g|<=1, |g'|<=2, |g''|<=10 give
\(\|\Pi_{ax}\|_{C^2_\eta}\le40P^2<C^4\).
Hence \(\Pi=\Pi_{ax}+C_p\) has J2 norm below C^257. The common axis
datum is retained, including its angular derivatives.

## 5. Propagating through the original stress coordinates

Use (4.16) exactly, with H=sqrt(2X)E and
\(W=1-(2D\eta M+dM_\eta)/X\):

\[
Q_s=-W+\frac{(1-h)I-D\eta I_\eta-dJ_\eta+2(h-D)\eta J}{XH},
\]
\[
N_s=-WU+\frac{D(M-\eta M_\eta)+4h\eta S-dS_\eta}{X}
+4A_{\rm phys}\eta\Pi-d\Pi_\eta,
\quad p_1=XQ_s/L,\quad p_2=XN_s/(LE).
\]

The explicit positive-exponent ledger gives J1 bounds as follows. Constants
and all coefficient functions in these formulas have J1 norm less than C.

| Quantity | Safe exponent q in J1 norm < C^q |
|---|---:|
| W | 262 |
| Qs, Ns | 337, 327 |
| p1, p2 | 351, 414 (use 512 for both) |
| a=1-2Ey/E, bs=2Uy/E | 140, 139 (use 256 for both) |
| 1/a | 301 |
| ts=-bs/a | 557 |
| vs=a(1+ts^2) | 1371 |
| Pc=p1+p2 ts, Jc=p2-p1 ts | 1070 |

Only first derivatives of stress coordinates are claimed. Angular moment
derivatives in (4.16) use the second angular derivatives proved above; a
first-derivative input must not be silently substituted. The sign ts=-bs/a
is the source convention; the shear vector is (a,-bs).

## 6. Original margins and decision

The preceding region arguments give Pc-2>.1 on this compact interval. On
activation after X-, the exact cutoff lower bound ea>=C^-401 and the
previous comparison ratio <C^-28<1/2 give third and fourth cone gaps at
least C^-401 and C^-802 wherever vs>=2. Use C^-1024 for both. The intervening
inner, shape, and matching regions have vs<1; original third/fourth positive
gaps are not required there. Their absence is why the loop is needed.

The first outer ramp has bs=0, a<=2 and p1>.5XR; any vs=2 point has both
remaining gaps greater than 1. The existing finite-radius axial certificate
has p1>35, |bs w|<.525247, bs^2<=2.048e-25 and Pc/p1>.737376.
For D0=bs^2/2 and A0=1-(2+D0)/p1, the exact fourth gap factors as

\[
2p_1^2(A_0+D_0/2)(A_0-bs\,w-D_0/2)>1.
\]

The factors exceed .94 and .417; the third gap also exceeds 1. In bounding
the first factor from below, drop its positive D0/2 term: substituting an
UPPER bound on D0 into that positive term would not give a lower bound.
On the two following bs=0 stages, p1>12 and the intermediate quadratic
bracket exceeds .3135, again giving both gaps greater than 1.
The retained left collar has vs-2>.29, and the right collar lies in the
power stage with vs-2=2lambda>C^-1, represented without rounding to zero.

Therefore the field/moment J2 norms, the input-coordinate J1 norms, the
listed positive reciprocals and normalization operators, and the VALUE
inverses of the required original cone and boundary margins all fit below
\(\boxed{\mathcal A=C^{4096}}\). This does not assert J1 bounds of every
cone-gap reciprocal or any derivatives of an as-yet-unconstructed loop.

Only `actual_compact_jet_envelope` is newly closed. The next step is to bound
the fixed-phase loop and inverse phase map, then transfer the modulation
through (4.16) to Cstate and D. D must still include the inverse-lambda loss.
The repair-to-state Lipschitz bound and its positive-field neighborhood also
remain to be established. An enormous proposed N cannot replace these proofs.

## Reproduce and review controls

```sh
python research/navier_stokes_cascade/results/compact_jet_envelope/test_audit.py -v
python research/navier_stokes_cascade/results/compact_jet_envelope/audit.py \
  --verify-record research/navier_stokes_cascade/results/compact_jet_envelope/evidence.json \
  --output /tmp/compact_jet_envelope_reproduced.json
```

`bounds.py` checks scalar comparisons at 80 and 110 interval digits and
records the exponent ledger. `jets.py` carries raw mixed derivatives without
inventing unavailable third derivatives. `review.py` compares the five exact
polynomial moment primitives and the integrated stresses against the ORIGINAL
source ODEs and independent differentiation. Deliberate controls expose lost
width factors, a missing exponential product term, frozen angular data,
dropped pressure terms, and a reversed shear sign. Manufactured fixtures are
never represented as numerical values or global suprema of the actual core.
Regression tests also retain the two issues corrected during this review:
integer division must not introduce binary floats into exact/high-precision
jets, and a positive term's upper bound must not be used as its lower bound.
