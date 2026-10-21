# Independent flat activation, five-moment and cone-collar recheck

Parent: `ce175f68ed606c259616c62eaa795ae83ce4500f`.
Status: `core-flat-activation-five-moment-and-cone-collar-independently-rederived`.

This checkpoint completes the next local validation step after the angular
source and two-sector exit recheck.  It independently rederives the flat
Volterra quotient, every forward moment difference, the complete pressure
terms in the stress map, the physical `C` gauge, and the quantitative cone
collar at the already frozen width `t1=C^-20`.

The conclusion is relative to the exact fixed point and exit gaps bound by
the preceding rechecks.  Exact arithmetic tests freeze the algebra and
scale inequalities but do not replace the continuum argument below.  No
continuation, matching root, compact modulation, repair, global pulse,
residual improvement, or blow-up conclusion is promoted here.  Historical
evidence is unchanged.

## 1. Positive flat Volterra operator, including the endpoint

For `0<=s<=1/4` write

\[
 \sigma(s)=e^{-1/s^2}g(s),\qquad
 g(s)=\{e^{-1/s^2}+e^{-1/(1-s)^2}\}^{-1},
\]

with `g(0)=e`.  For every nonnegative integer `m`, define

\[
 {1\over\sigma(s)}\int_0^s\sigma(u)u^m b(u)\,du
       =s^{m+3}{\cal V}_m[b](s).
\]

The substitution `v=u^-2-s^-2` gives the exact positive formula

\[
 {\cal V}_m[b](s)=\frac1{2g(s)}\int_0^\infty e^{-v}
 (1+s^2v)^{-(m+3)/2}
 g\!\left({s\over\sqrt{1+s^2v}}\right)
 b\!\left({s\over\sqrt{1+s^2v}}\right)\,dv.              \tag{1}
\]

Since `2<=g<=6` on this interval, (1) has operator norm at most `3/2`
on values and on every eta jet: `g` has no eta dependence.  Polynomial
`v` factors produced by fixed `s` derivatives remain integrable against
`e^-v`, so differentiation under the integral gives the smooth endpoint

\[
 {\cal V}_m[b](0)=b(0)/2.                                 \tag{2}
\]

For `m=0`, direct differentiation also gives

\[
 {\cal V}_0[b]'(0)=b'(0)/2,\qquad
 {\cal V}_0[b]''(0)=\{b''(0)-3b(0)\}/2.                  \tag{3}
\]

The `-3b(0)/2` radial-Jacobian term is retained.  Assigning zero to the
flat `0/0` quotient, or deleting this term, would invalidate the edge jets.

## 2. Field differences and all five forward moments

Put `y=t s`, `ea=(1-kappa0)sigma(s)`, and
`alpha=t ea s^3`.  Subtracting the activated slopes from the exact seed
slopes and applying (1) with `m=0` gives

\[
 \Delta\log E=\alpha b_E,\qquad \Delta U=\alpha b_U.     \tag{4}
\]

The entire `exprel` function then factors `Delta E`, `Delta H`,
`E_r/E-1`, and `Delta(E^2)` by the same `alpha`, with no division by a
flat number.  One further integration uses `dy=t ds` and (1) with `m=3`.
It gives exactly

\[
 (\Delta M,\Delta I,\Delta J,\Delta S,\Delta C_p)
    =t^2ea s^6(b_M,b_I,b_J,b_S,b_C).                      \tag{5}
\]

The five coefficients before applying `V_3` are:

| row | exact coefficient after removing `alpha` |
|---|---|
| `M` | `X bU` |
| `I` | `X bH` |
| `J` | `X(Ur bH + Hr bU + alpha bU bH)` |
| `S` | `X(2Ur bU + alpha bU^2 - bE2/2)` |
| `Cp` | `bE2/2` |

The pressure row has no `X`: in logarithmic radius,
`(Cp)_y=E^2/2`.  Conversely, deleting the `alpha bU bH` term from `J`
or the `bE2` terms from `S,Cp` would destroy the exact subtraction.
All eta derivatives commute with (1) and (5).

## 3. Complete moment-map and pressure subtraction

Use the original forward map

\[
 W=1-(2D\eta M+dM_\eta)/X,
\]
\[
 N_Q=(1-h)I-D\eta I_\eta-dJ_\eta+2(h-D)\eta J,
 \qquad Q_s=-W+N_Q/(XH),                                  \tag{6}
\]
\[
 N_s=-WU+\{D(M-\eta M_\eta)+4h\eta S-dS_\eta\}/X
       +4A\eta\Pi-d\Pi_\eta,                             \tag{7}
\]
\[
 p_1=XQ_s/L,\qquad p_2=XN_s/(LE),\qquad\Delta\Pi=\Delta C_p.
\]

The reciprocal subtraction in (6) is

\[
 \Delta Q_s=-\Delta W+{\Delta N_Q\over XH}
 -{N_{Q,r}\Delta H\over XH_rH}.                           \tag{8}
\]

The complete axial subtraction is

\[
 \begin{aligned}
 \Delta N_s={}&-W_r\Delta U-U_r\Delta W-\Delta W\Delta U\\
 &+\{D(\Delta M-\eta\Delta M_\eta)
       +4h\eta\Delta S-d\Delta S_\eta\}/X\\
 &+4A\eta\Delta C_p-d\Delta(C_p)_\eta,
 \end{aligned}                                             \tag{9}
\]

and the field denominator in `p2` contributes

\[
 \Delta p_2={X\over L}\left({\Delta N_s\over E}
             -{N_{s,r}\Delta E\over E_rE}\right).        \tag{10}
\]

Thus both exceptional pressure terms, the `H` reciprocal difference, and
the `E` denominator difference survive.  Equations (4)--(10) imply

\[
                         \Delta p_s=t ea s^3 b_p.          \tag{11}
\]

The stronger factor in (5) cannot be assigned to (11), because the local
field-denominator variations are already order `t ea s^3`.

## 4. Physical `C` gauge and explicit envelopes

Gauge

\[
 \bar E=CE_r,\quad \bar H=CH_r,\quad \bar q=q_r/C.
\]

Then `E_r q_r=bar E bar q=-2(U_r)_y`.  Consequently `M,S,W,N_s,p1`
have zero physical `C` power, `I,J,N_Q` have power `C^-1`, `C_p` has
power `C^-2`, and `p2,q_r,t_s` have one positive `C` power.  In
particular, pressure memory becomes smaller under this gauge; it is not
discarded.

Set `J0=1+Lambda+K`, `B=10^20 J0^64`, and `R=exp(B/100)`.  Recounting the
products in (4)--(10), with the positive Volterra norm `3/2`, gives the
following sufficient bounds after removal of the displayed factors:

| normalized difference | bound |
|---|---:|
| `Delta logE/(t ea s^3)` | `R` |
| `Delta U/(t ea s^3)` | `R^2` |
| `C Delta E, C Delta H` divided by `t ea s^3` | `R^4` |
| `Delta M, C Delta I, C Delta J, Delta S, C^2 Delta Cp` divided by `t^2 ea s^6` | `R^8` |
| `Delta W, C Delta NQ` divided by `t^2 ea s^6` | `R^10` |
| `Delta p1/(t ea s^3)` | `R^18` |
| `Delta Ns/(t ea s^3)` | `R^16` |
| `Delta p2/(C t ea s^3)` | `R^24` |
| `(Er/E-1)/(t ea s^3)` | `R^2` |

For example, the largest exact `J` integrand terms cost `R^5`, `R^3`,
and `tR^6`; the already selected `tR^3<.1` makes the last at most a
small multiple of `R^3`.  Hence `R^8` is not being inferred from a
fixture.  In (8), the three terms cost at most the displayed `R^18` after
the `X/L` factor.  Equation (9), including both pressure terms, fits
`R^16`; (10) then fits `CR^24`.

## 5. Cone coefficients and the unchanged width

The cancellation-safe identities are

\[
 {P_c-v_s\over ea}=v_r+t s^3D_c,\qquad
 {J_c\over ea}=t s^3J_c^*,\qquad
 {v_s-kv_r\over ea}=t s^3V_c.                              \tag{12}
\]

After division by `t s^3`, the worst product is

\[
 {t_s\over C}\,{\Delta p_2\over C t ea s^3},
\]

which costs `C^2 R^(3+24)=C^2R^27`.  The other terms cost at most
`R^18`, `C R^24`, or `C^2R^6`.  Allowing sums and eta product rules,

\[
 |D_c|,|J_c^*|,|V_c|\le K_{act}:=C^2R^{80}.                \tag{13}
\]

This count keeps the `C^-2` pressure increments in (9); it does not replace
them by zero.

The selected amplitude satisfies

\[
 \log C=10000B+100\log(1/\epsilon),\qquad \log R=B/100.
\]

Therefore `R^80<C`, so `Kact<C^3`.  With the unchanged
`t1=C^-20`,

\[
 t_1R^3<.1,\qquad K_{act}t_1<C^{-17},\qquad
 t_1C^2R^8<.05.                                           \tag{14}
\]

The last bound propagates the exact exit growth `v_r>2.3` to
`v_r>=9/4` on the native collar.  Put

\[
 s_c={1\over10\sqrt{\log C}},\qquad t_c=t_1s_c.
\]

Because `s_c<=1/4`,

\[
 ea\le\exp(16/9-1/s_c^2)=e^{16/9}C^{-100},
 \qquad ea(C^3+1)<1/8.                                    \tag{15}
\]

Using (12)--(15), `v_r<=C^3`, and `Kact t1<C^-17` gives

\[
 v_s\ge17/8>2,\qquad {P_c-v_s\over ea}\ge9/8.            \tag{16}
\]

Moreover `J_c/ea=t_1s^3J_c^*`; after squaring, `R^160<C^2`
and the factor `v_s-2<=C^3` leave exactly `C^-31`.  Hence

\[
 { (v_s-2)J_c^2\over2(P_c-v_s)^2}
 \le {32\over81}C^{-31}s^6<1/2.                           \tag{17}
\]

This proves the local quantitative cone collar at the already frozen
width.  It neither selects nor validates the subsequent continuation and
matching roots.

## 6. Reproduction and frontier

From the repository root:

```bash
python3 research/navier_stokes_cascade/results/core_activation_volterra_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_activation_volterra_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_activation_volterra_recheck/evidence.json \
  --output core_activation_volterra_recheck_replay.json
```

The next obligation is to independently rederive the prescribed
continuation and both matching-root constructions on the same selected
profile.  Only after those roots are fixed can compact modulation, exact
repair, and the annular/exterior cone join be independently audited.
