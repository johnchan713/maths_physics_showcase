# Independent angular source and two-sector core-exit recheck

Parent: `3476acc9ccdb5921e20964f86180128fac094ad6`.
Status: `core-angular-source-and-two-sector-exit-independently-rederived`.

This checkpoint independently rederives the next part of the frozen local
core argument after the all-index fixed-point ledger: the scalar comparison,
the positive angular source without division by `chi`, and the angular and
pressure-driven axial exit sectors.  It uses the sharper shared-pressure
separation already proved in `core_pressure_gap_recheck`, but otherwise
restarts the estimates from the exact source formula.

The result is an analytical local lemma for the exact fixed point.  The
rational program checks every displayed constant and freezes provenance; it
is not a continuum proof by sampling.  This checkpoint does **not** reprove
the flat activation Volterra estimates, continuation, matching roots,
compact modulation, exact repair, outer joins, complete physical residual,
or blow-up.  No historical evidence is modified.

## 1. Scalar comparison and the weighted parameter derivative

Put

\[
 f_0(z)=\sum_{n\ge0}\frac{(-z/2)^n}{n!(n+1)!},\qquad
 z=Y\chi,\qquad 0\le Y\le 41/10.
\]

With `t=z/2<=41/20`, the alternating terms decrease after the first
pair.  The odd partial sum is therefore a lower bound:

\[
 f_0(z)\ge1-\frac t2+\frac{t^2}{12}-\frac{t^3}{144}
 \ge \left.1-\frac t2+\frac{t^2}{12}-\frac{t^3}{144}
       \right|_{t=41/20}>\frac{53}{200}.
\]

The derivative series has decreasing absolute terms from its first term,
so `-1/4<=f0'<0`, while `f0<=1`.  Consequently

\[
 \left|Y\frac{f_0'(Y\chi)}{f_0(Y\chi)}\right|
 \le \frac{(41/10)(1/4)}{53/200}=\frac{205}{53}<4.       \tag{1}
\]

The exact angular weight identities give

\[
 |H_*\chi_\eta|\le106\chi,\qquad
 |\chi_\eta|\le106/\sigma.
\]

It is essential to retain the first, weighted identity before estimating
the parameter derivative.  From (1),

\[
 |D_Y\log f_0(Y\chi)|\le4\chi,qquad
 |H_*\partial_\eta\log f_0(Y\chi)|\le410\chi.            \tag{2}
\]

For the exact fixed point write `Phi=f0+r`, where the all-index contraction
gives the required `Y,eta` jets of `r` bounded by `e` and `Phi>=1/4`.
The reciprocal-difference identity yields

\[
 |D_Y\log\Phi-D_Y\log f_0|\le32e.
\]

For the weighted parameter derivative the two error terms cost

\[
 4\cdot64e+4\cdot410e=1896e<2000e.                       \tag{3}
\]

Thus

\[
 |H_*\partial_\eta\log\Phi|\le500\chi+2000e,
 \qquad |\Phi_\eta/\Phi|\le500/\sigma.                  \tag{4}
\]

Using only the unweighted estimate in the first inequality of (4) would
lose the vanishing at `H_*=0` and would not justify the source bound below.

## 2. Complete source ledger

The exact normalized source is

\[
 S_q=-W(1+D_Y\log\Phi)-h(1-2\eta U)
     -H_c(\Lambda\zeta_*+\Phi_\eta/\Phi).                \tag{5}
\]

On the real interval, `-W_*>2.9`, `|W|<=5`, `|U|<4.1`,
`|W-W_*|<=2 Ubar/Lambda`, and
`|zeta_*|<=sqrt(chi)/sigma`.  Substitution of (2)--(4) into
(5) gives the ten-group ledger

| group | lower contribution |
|---|---:|
| main angular term | `+Lambda L chi` |
| constant part of `-W` | `+2.9` |
| small `h` term | `-0.1` |
| `W` scalar logarithm | `-20 chi` |
| weighted `H_*` logarithm | `-500 chi` |
| `W` fixed-point error | `-160 e` |
| `H_*` fixed-point error | `-2000 e` |
| constant `W-W_*` error | `-2 Ubar/Lambda` |
| `H_c-H_*` times `Lambda zeta_*` | `-(Ubar/sigma)sqrt(chi)` |
| `H_c-H_*` times `Phi_eta/Phi` | `-500 Ubar/(Lambda sigma)` |

No pressure term is removed from the fixed-point source: its effect is
already in the exact `u` entering `W,H_c`.  The ledger gives

\[
 S_q\ge\Lambda L\chi+2.8-520\chi
 -(\bar U/\sigma)\sqrt\chi-2160e
 -2\bar U/\Lambda-500\bar U/(\Lambda\sigma).             \tag{6}
\]

Young's inequality, with `L>=.98`, gives

\[
 (\bar U/\sigma)\sqrt\chi
 \le .03\Lambda L\chi+9\bar U^2/(\Lambda\sigma^2).      \tag{7}
\]

Let `H=K^(1/16)>1`.  The frozen parameter definitions imply
`Ubar<=H/10^6`, `E Fb<=H/10^6`, `sigma^-1<=H/10^6`, and
`Lambda=10^44 H^128`.  Hence the constant remainder in (6)--(7) is at most

\[
 \left(\frac{2162}{10^6}+\frac{500}{10^{12}}
             +\frac9{10^{24}}\right)\frac{H^4}{\Lambda}
 <3\,10^{-47}<.1.                                       \tag{8}
\]

Also `520 chi <= .02 Lambda L chi`.  Combining (6)--(8), without
ever dividing an error by `chi`, proves

\[
                         S_q\ge2.5+.95\Lambda L\chi.     \tag{9}
\]

In particular (9) remains meaningful at `chi=0`.

## 3. Positive forward shear

The exact radial equation and the regular-axis condition give

\[
 \frac{p_1}{X}=\frac{Q_s}{L}
 =\frac1L\int_0^1v\frac{\Phi(vY)}{\Phi(Y)}S_q(vY)\,dv.  \tag{10}
\]

Because `1/4<Phi<2`, the quotient is at least `1/8`.  Equations
(9)--(10), `L<=1`, and the constant part of (9) give

\[
 \frac{p_1}{X}\ge\frac18\frac52\int_0^1v\,dv
 =\frac5{32}>\frac1{10}.                                \tag{11}
\]

The radial logarithmic estimate in (1), plus the `32e` error, also gives
`0<p1<=10` throughout this core.  At `X=4/Lambda`, (11) implies

\[
                         a_r=p_1\ge2/(5\Lambda)>0.        \tag{12}
\]

## 4. Angular and axial exit sectors

If `chi>.99` at `Y=4`, then `t=2chi` lies in `(1.98,2]`.  The series

\[
 f_0(z)+zf_0'(z)=\sum_{n\ge0}\frac{(-t)^n}{(n!)^2}
\]

has decreasing absolute terms after `n=1`.  Truncation after the positive
fourth term is therefore an upper bound.  Its derivative is negative on
`[1.98,2]`, and at `t=1.98` it equals exactly

\[
 1-t+\frac{t^2}{4}-\frac{t^3}{36}+\frac{t^4}{576}
 =-\frac{75535511}{400000000}.                            \tag{13}
\]

Since `f0<=1`, (13) makes the comparison shear exceed `2.376`.
The actual logarithmic shear differs by at most `64e<.001`, so

\[
                            p_1>2.3\qquad(\chi>.99).      \tag{14}
\]

If `chi<=.99`, the independently sharpened pressure separation gives
`|Z_*|>delta_*`, where `delta_*=jP^2/100`.  The fixed-point error satisfies
`|n_s-Z_*/L|<=2e<delta_*/2`, hence `|n_s|>=delta_*/2`.
With physical axial shear `p2=q=-2U_y/E=X n_s/E`,
`E^2<=8X exp(4 Lambda M)/C^2`, `p1<=10`, and `X=4/Lambda`,

\[
 \frac{p_2^2}{p_1}
 \ge\frac{C^2\delta_*^2}{80\Lambda e^{4\Lambda M}}.     \tag{15}
\]

This is where the physical `X/E` factor matters; omitting it would give the
wrong power of the core scale.  From
`log C=10000B+100log(1/epsilon)`,
`B=10^20(1+Lambda+K)^64`, `M<Lambda`, `P>=16`, and
`Lambda>10^44`, the logarithm of the right side in (15) is greater than

\[
 20000\Lambda^{64}-4\Lambda^2-\log\Lambda-105>\log(2.3).
\]

Thus `p2^2/p1>2.3` in the complementary sector.  Together with
(12)--(14), both sectors prove the exact native-exit statement

\[
 a_r\ge2/(5\Lambda)>0,\qquad
 v_r=a_r+q_r^2/a_r>2.3.                                  \tag{16}
\]

This is still only the leading local core exit.  It does not promote the
inherited activation or compact-repair claims.

## 5. Reproduction

From the repository root:

```bash
python3 research/navier_stokes_cascade/results/core_growth_source_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_growth_source_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_growth_source_recheck/evidence.json \
  --output core_growth_source_recheck_replay.json
```

The first command runs exact rational regression checks.  The second binds
the note, checker, protocol, three frozen dependency records, and scope
flags into a fresh replay without overwriting frozen evidence.
