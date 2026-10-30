# Q21: affine extended-domain to whole-path source transfer

Parent: `2f844abbf80dda5931a2143fe22ffb987698efdc` (Q20).
The next source bound is reduced to a precise extended-domain inequality.
This checkpoint proves the derivative transfer along the original lift; it
does not prove that inequality for every actual generated source. Historical
geometry, fields, evidence, pressure and widths remain unchanged.

## 1. Use the original lift, with frozen discrete data

In `auxiliary_geometry_audit` (15)--(16), put b=sqrt(2)-1,
v_t=(b,1), lambda_t(H)=(b H_1+H_2)/(1+b^2), Delta=i-i_* in {0,1,2,3}.
On each fixed common-torus lift the path is

    v(H)=(lambda_t(J^Delta H-c_color-k_copy)+r0)/c_i,
    H_w=H+Tg^(-Delta)c_i(w-v(H))v_t, 0<=w<=L,
    L=2r0/c_i, c_i<=S_i^-1<=S_*^-1.                  (1)

The last bound uses the coarsest common band, S_i>=S_*. All labels,
covering indices, c_i, r0, centers and copies are fixed when coefficient
jets are taken. Physical slow coordinates (R,Z,T) are unchanged along
the path. Formula (1) is evaluated separately on each lift, never by
averaging preimages. Its band image stays in the same closed enlarged
rectangle, with the same transverse coordinate and pulse time w.

Since lambda_t J^Delta=Tg^Delta lambda_t, differentiation gives

\[
 D_HH_w=\Pi=I-v_t\lambda_t
 ={1\over1+b^2}\begin{pmatrix}1&-b\\-b&b^2\end{pmatrix},
 \quad \partial_w H_w=Tg^{-\Delta}c_i v_t.             \tag{2}
\]

All second and higher derivatives of this affine map vanish, including
mixed H,w derivatives. Each column of Pi has absolute sum at most 2,
since 0<b<1. The w column has sum at most 2/S_*. Slow-coordinate columns
are identities. The Pi projection fixes the transverse direction and
annihilates v_t; it is not the incompressibility projection in the ODE.

## 2. A quantitative whole-path jet estimate

Let F(x,H) be a source on the entire enlarged extended rectangle, where
x denotes the unchanged slow variables. Suppose for all |I|=m<=R

\[
 |D^I F(x,H)|\le\epsilon^a\sqrt\zeta\,d^{-p_m}
                  P(v(H),x)A_m(S_*).                  \tag{3}
\]

Use a common finite pole if needed. Chart/band normalization factors,
all logarithmic degrees and original profile constants belong in A_m.
Assume (3) at every lift and point reached by (1), not only on a physical
trace. Let k auxiliary derivatives, r path derivatives and l slow
derivatives be applied to F(x,H_w). The affine chain rule gives

\[
 |D_x^lD_H^k\partial_w^r F(x,H_w)|
 \le 2^{k+r}S_*^{-r}\epsilon^a\sqrt\zeta\,d^{-p_m}
                   P(w,x)A_m(S_*),\quad m=l+k+r.       \tag{4}
\]

For multiindices the same bound holds by expanding each derivative into
a linear combination of native derivatives and multiplying the column
absolute sums. No higher derivatives of the map appear. The radial edge
and epsilon factors are unchanged because the slow coordinates are fixed.
Crucially P is evaluated at the path point: v(H_w)=w. There is no ratio
P(w)/P(v(H)) to bound and no reciprocal Gaussian minimum.

Equation (4) requires only native source jets through order m. Its factor
is uniform in band, lift and copy, and introduces no q loss or positive
scale power. If the native estimate already used a covering normalization,
that cost remains in A_m and is not silently discarded. The statement is
in coefficient coordinates; physical graph derivatives still incur their
separate established costs. Moving labels or differentiating c_i as a
continuous variable are outside (4).

Q20 now converts the weighted derivatives in (4) into derivatives of
F(x,H_w)/P(w,x), using the existing log-Gaussian jets. Thus a proved
extended W bound supplies the normalized whole-path inverse input with
no additional source derivative order. Finite constants/log degrees can
increase under the Q20 convolution. The original inverse-normal and
propagator estimates are still needed to solve the response.

## 3. Source support and full pressure

The same path formula preserves the transverse coordinate. A source which
vanishes, with the required jets, on an entire transverse/slow path gives
zero zero-entry response there by uniqueness. Vanishing at just the entry
or at one physical point is insufficient. Endpoint cutoffs must still be
applied to both potential and pressure; a zero-entry solve need not vanish
at the outgoing endpoint.

Keep the unprojected source in

    pi=i(n.Kbg.t-nprime.t+n.f)/(k m |n|^2).

Equation (4) applies to full F, not merely its tangent projection. For
f=n g, the projected velocity forcing may be zero while pressure is
nonzero. The complete mean A0/A1 velocities, moving pressure moments and
Pi1/Pi2 memory remain in their respective source formulas; this coordinate
transfer does not cancel any of them. Additive Gaussian/base/exterior tails
remain outside the forwarded regular source.

## 4. The remaining actual inequality and scope

The unclosed input is now exactly (3), for each actual recomputed source
and each finite Q14 derivative/pole request, on the full common extended
rectangle and all accepted bands. Together with the retained inverse
hypotheses, (4) and Q20 discharge its conversion to the whole-path norm.
The domain coverage is supplied by the frozen geometry; the magnitude of
the actual source is not supplied by that geometry.

This strengthens the interface between Q19's typed source grammar and the
original first-sourced inverse (9)--(11). It does not prove (3) by treating
an assigned epsilon rate as a norm bound. It also does not establish every
Q16 component rate, complete residual decay, endpoint closure, force
regularity or blow-up. The next obligation is to bound the actual extended
source formula (starting with the first sourced block) in the corrected
half-weight space, with all linear pressure/viscous and mean-wave terms.

Exact rational fixtures test affine composition, path time, transverse
invariance, mixed derivatives, uniform column bounds and lift offsets.
They use rational b in (0,1); the analytical proof covers the actual
b=sqrt(2)-1 without approximating it. They are not PDE or continuum source
enclosures. The parent evidence and full inherited hash graph are replayed.

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q21_lifted_path_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q21_lifted_path_join/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q21_lifted_path_join/evidence.json --output FRESH.json
```
