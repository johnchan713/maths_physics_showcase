# A global third-angular axial bound and explicit signed pressure-source input

Parent: `0ba0d99a92c7e6de0158f91a11a49709f805ec64`.
This checkpoint supplies the global `K_U=sup_[0,Xv]|U0|_3` input left by
[signed_pressure_source_moments](../signed_pressure_source_moments/README.md),
relative to the retained exact core and regional estimates. It joins the
already proved compact modulation/repair C3 bound to an explicit axis-side
bound and NEW third-order bounds for the same later axial pulse strength.
It does not repeat the compact-loop proof or substitute a sampled profile.

The complete post-pulse energy, including angular corrections and the
unchanged waiting/tail stages, is used in the strength equation. Its tiny
affine end-bump terms are retained. All parameters, roots, gauges, support
endpoints, width C^-20 and frequency remain unchanged. Historical evidence
is not edited. Passing tests are not a continuum proof, and the retained
core/regional premises are not independently proved unconditional here.

## 1. Raise the SAME reference angular root to C3

Use factorial norms `|f|_k=sum_(j<=k) sup|partial_eta^j f|/j!`.
The actual angular map from root_foundation has CONSTANT B,Q,
`||B^-1||<3`, `||Q||<32`, and target
`z=(-deltaI exp[-beta(Tu-3)],0)`, with `beta=1-lambda`,
`Tu=30 log(1/lambda)`, `lambda=exp(-4T)<=exp(-256)`, `T>=64`.
Its exact interpolation discrepancy is

\[
 \delta I=e^{-\beta T_f+a}\delta_0+
 {a\over\beta}\int_0^1e^{-\beta T_f(1-z)+a\sigma(1-z)}\sigma'(z)dz,
 \quad a=\log{2\over1+\eta^2},\quad T_f=1000.             \tag{1}
\]

Here delta0 is eta-independent with `|delta0|<=1/beta`.
On [-1,1], `|a|<=1`, `|a'|<=1`, `|a''|<=2`, `|a'''|<=12`,
so `|a|_3<8`. The Banach exponential bound, `beta>1/2`,
`int sigma'=1` and `exp(8)<2^16` give `|deltaI|_3<2^21`.
Since `30 lambda log(1/lambda)<1`,

\[
 |z|_3\le2^{29}\lambda^{30},\qquad
 |c_{ang}|_3\le2^{32}\lambda^{30}=:\mathfrak c.          \tag{2}
\]

Indeed `exp[3beta+30lambda log(1/lambda)]<exp(4)<2^8`, and the
quadratic contraction condition is at most `2^41 lambda^30<1`.
This is a C3 algebra proof, not relabelling the old C1 target.

The C3 root is the SAME root: the old root has pointwise norm
`<=4200 lambda^30<mathfrak c`; the new root has norm <=mathfrak c.
The exact pointwise map has Lipschitz constant at most
`192 mathfrak c<1` on that larger ball. Thus any two roots in it agree.
It is not necessary, and would be incorrect, to assert the new larger
C3 radius is contained in the old C1 radius. No parameter is reselected.

## 2. Bound the COMPLETE post-pulse energy in C3

Keep `R=int_(Xend)^infinity E^2 dX/(Xend Eend(eta)^2)` for the POWER
reference before the heat attachment. The retained full value bound is
`0<R<=256T`. Its relative angular edit is below .01.

On interpolation the normalized angular shape is
`g_theta=((1+eta^2)/2)^theta`, `theta=2sigma(y/Tf) in [0,2]`.
Writing p=1+eta^2 gives the exact derivatives

\[
 g'=2\theta\eta g/p,\quad
 g''=g[2\theta/p+4\theta(\theta-1)\eta^2/p^2],
\]
\[
 g'''=g[12\theta(\theta-1)\eta/p^2+
                 8\theta(\theta-1)(\theta-2)\eta^3/p^3]. \tag{3}
\]

Their absolute upper bounds `(1,4,12,56)` imply `|g_theta|_3<32`.
Thus interpolation contributes at most `32Tf<=512T`.
All UNEDITED subsequent normalized energy has the form `C p^2`, C>=0,
because all those endpoints, the reference wait and scalar shapes are
eta-independent. Since `|p^2|_3=24` and `max p^2=4`, its C3 norm is six
times its value at eta=1. Pointwise relative-edit positivity yields
`R_unedited<2 R<=512T`; hence the unedited subsequent C3 contribution
is at most `3072T`. Together the unedited norm is below `4096T`.

For each of the two angular bumps, its positive weighted linear mass is
<=.3 and its quadratic mass <=20. The scalar attenuation factors are
<=1. Therefore the COMPLETE signed edit has C3 norm at most
`24(1.2 mathfrak c+40 mathfrak c^2)<=29 mathfrak c+960 mathfrak c^2<1`,
using `mathfrak c<2^-40`. No signed bump energy is omitted. Consequently

\[
                         \boxed{|R|_3\le8192T.}          \tag{4}
\]

This argument keeps the infinite tail and waiting-stage mass; it does not
multiply an infinite interval by a uniform value. The later Gamma heat
edit and its compensation have exact zero net S change as ETA FUNCTIONS,
so (4) supplies the same strength equation after that edit as well.
It is not a pointwise power-tail formula for the heat-modified profile.
Nor is (4) asserted for the heat-modified R by itself: it is the COMPLETE
strength equation, with both the earlier compensation and later heat edit,
that stays identical as a function of eta and trial strength.
The angular root and post-swirl energy are independent of axial strength,
so this estimate is obtained BEFORE differentiating the strength root.

## 3. The retained positive strength and both affine end corrections

The axial_root_foundation solve has
`c_i=a_i Amp+b_i(eta)`, where a_i are eta-independent and
`b_i=v_i eta(1+eta^2)`. Its uniform trial bound includes Amp=0, and gives
`|a_i|<=lambda^20`, `|b_i|_1<=lambda^20`.
Since the polynomial has C1 norm 6 and C3 norm 10,

`|b_i|_3<=2lambda^20`.

The incoming s0 has exact form `s0=scoef eta^2(1+eta^2)^2-Gp/2`.
The retained prefix bounds give `0<=scoef<1/4`,
`Gp/2<=241T+17/6`. The polynomial has C3 norm 76, so `|s0|_3<512T`.

Keep the full energy equation

\[
 F=k\,Amp^2+B(\eta)Amp+D(\eta)=0,\quad
 k=K+\lambda\sum w_i a_i^2,
\]
\[
 B=2\lambda\sum w_i a_i b_i,\quad
 D=-C+\lambda s_0-\tfrac\lambda2e^{-26}R+
                         \lambda\sum w_i b_i^2.        \tag{5}
\]

Here K is the ORIGINAL shape integral, C=(1-exp(-26))/4, and `0<w_i<20`.
The parent gives `K<1/4`, hence `k<=1/4+40lambda^41<1`, together with
`Amp<1.2` and inverse Jacobian `<100/49<3`.
The odd cross term B is generally nonzero; neither parity nor evenness
of Amp is assumed. Equations (4),(5) imply

`|B|_3<=160lambda^41`, `|D+C|_3<=8192Tlambda=:epsilon`.

The latter follows from `512+4096+1<8192`; the end-square contribution
`160lambda^41` is below `Tlambda`. Also `160lambda^41<=epsilon`.
Monotonicity of `T exp(-4T)` on T>=64 gives
`epsilon<2^19*2^-256=2^-237<2^-20`.

For normalized jets u_n=Amp^(n)/n!, the exact SAME-root recursion is

\[
 (2kAmp+B)u_n=-k\sum_{j=1}^{n-1}u_ju_{n-j}
                    -\sum_{j=1}^{n}B_ju_{n-j}-D_n.      \tag{6}
\]

Using u0<=2 and the above coefficient jets proves successively
`|u1|<=9epsilon`, `|u2|<=279epsilon`, `|u3|<=15939epsilon`.
These are normalized jets; no factorial is lost. Thus

\[
 |Amp|_3\le1.2+16227\epsilon<2,\qquad
                         |c_i|_3\le4\lambda^{20}.       \tag{7}
\]

Only fixed third order is claimed. This is not an order-uniform estimate
and does not change the previously selected root or ignore its affine terms.

## 4. Global radial cover and pressure-source interface

On the later pulse `U=eb f exp[-(.5+lambda)y](Amp R0+sum c_i beta_i)`.
Here `eb<1`, `|f|_3<16`, `R0<=11`, `beta_i<=64`.
Main and end-bump supports are pairwise disjoint. Equations (7) give

\[
                              |U|_3\le16\cdot22=352.    \tag{8}
\]

The other reference axial taper is `U=k(y)eta`, |k|<=4, hence norm <=8;
all later post-pulse/heat stages have U=0.

For the axis-side cover use the formula K in collar_operator_envelope
equation (4), now through `ell_axis=tc/2`, not just tc/4.
The actual chosen C>1e9 implies log C>1, so `ell_axis<t/20`.
Its identical proof and zero-free strip apply through the modulation START,
with `4exp(ell_axis)<4.1`; the modulation is flat at that endpoint.
No coefficient cutoff is enlarged and no radial analyticity is assumed.
With delta=rho/(800S), define

`K_axis=K(ell_axis)*(1+delta^-1+delta^-2+delta^-3)`.

Angular Cauchy cancels the factorials in |.|_3 and gives the actual U bound
on `[0,Xa exp(tc/2)]`, including the axis.
The existing higher_jet_budget C3 bound `2A_hi`, `A_hi=(C^4096)^256`,
applies on `[Xa exp(tc/2),rho_p exp(5)]` to the SAME modulation/repair root.
Matching-root uniqueness preserves its identity; that historical annular
bound is not silently extended to later radii. The reference taper, zero
power stages, later axial pulse and zero post-pulse field cover the rest.

\[
 \boxed{K_U\le K_U^*:=\max\{K_{axis},2A_{hi},352\}.}    \tag{9}
\]

This is an explicit finite expression in the fixed construction parameters,
not their numerical materialization or an independent unconditional proof
of all inherited core/regional estimates. The original Xv stays unchanged.
Substitute `u2=K_U^*`, `v2=36K_U^*`, nu=1 and that Xv in the preceding
Q0,Q1 formulas: this supplies the two full signed pressure-source bounds
for the unchanged first-coefficient repair. The axis-viscosity term remains.
No incoming pressure datum, Pi1/Pi2 reserved-mean memory or fifth moment
is reset; signed cancellation is still NOT an absolute-integral bound.

## 5. Checks and remaining limits

Exact rational budgets verify the contraction, same-root ball comparison,
complete energy constants, normalized third-root recurrence and region join.
Symbolic differentiation checks all three interpolation-power derivatives
and every term of the third strength derivative. Five manufactured positive
quadratic roots independently reproduce their first three normalized jets,
including eta=+/-1; they are not the actual enormous-parameter root.
Parent tests and evidence are replayed under their original comparators;
no frozen record is rewritten. Passing tests are not a continuum proof.

Next quantify the separate pointwise w0/shifted-axial sources, compact
leading derivatives and actual cutoff/repair-bump envelopes in Cstar.
The SAME fixed-tilt test and all scale/freezing/geometry/inverse gates
remain. No historical first-band, unconditional pulse, complete residual
improvement or blow-up conclusion follows from (9).

```sh
python3 research/navier_stokes_cascade/results/global_axial_angular_bound/test_audit.py -v
python3 research/navier_stokes_cascade/results/global_axial_angular_bound/audit.py \
  --verify-record research/navier_stokes_cascade/results/global_axial_angular_bound/evidence.json \
  --output /tmp/global_axial_angular_bound_fresh.json
```
