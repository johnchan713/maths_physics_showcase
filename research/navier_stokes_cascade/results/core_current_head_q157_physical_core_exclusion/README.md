# Q157: same-physical core path and conditional competitor exclusion

Parent: Q156, `5384be92235f17e51a1e5206b08bb767831a7fde`.
The live branch, Q156 evidence and CI, the frozen same-field core theorem,
Q149 physical-primary binding, Q153 physical diagonal, Q155 localization
and the certified complete-potential support theorem were read first.
Historical evidence, profiles, roots, pressure gauge, supports, accepted
bands, diagonal cutoffs and `C^-20` width are unchanged.

Q156 correctly refused to import an older quantitative core path without
first identifying it with the later physical field. This checkpoint supplies
the missing identification. It does not identify the canonical and physical
primary pulse families: the core path only sees their shared finalized Borel
background, while every physical diagonal increment vanishes in the same
annular support hole.

The conclusion is conditional on all named continuum providers. Exact tests
reproduce arithmetic and frozen provenance; they are not a continuum proof.

## 1. Fixed path and positive leading coefficient

Keep the unchanged parameters

\[
 h={1\over200},\qquad A={1\over2}+h={101\over200},
 \qquad X_a={4\over101},\qquad X_{\rm in}={1\over101}<X_a. \tag{1}
\]

For `tau=1-t`, choose

\[
 x_\tau=(\sqrt{2X_{\rm in}\tau},0,0).                 \tag{2}
\]

At `z=0`, the implicit chart equation gives exactly

\[
 q=\tau,\qquad\eta=0,\qquad X=X_{\rm in}.             \tag{3}
\]

The protected Cartesian core ansatz then has no radial contribution to the
azimuthal component and gives

\[
 q^A u_{\theta,0}(x_\tau,1-\tau)
 =E_0(X_{\rm in},0)=:e_0.                             \tag{4}
\]

The frozen exact-core provider gives `Phi>33/125`; the exponential factor
is positive and the amplitude denominator is positive. Hence

\[
 e_0=\sqrt{2X_{\rm in}}{\varphi_*\Phi\over C}>0.       \tag{5}
\]

This is the leading velocity coefficient identity. The local shear and
stress-growth diagnostics are not used as a substitute for (4).

## 2. Why the result applies to the physical family

Three exact interfaces are needed.

First, Q149 binds the inherited physical primary family to the **same**
finalized leading profile and rule-defined Borel background used by the
positive-core theorem. It preserves the original profile, Borel scales,
supports and pressure gauge. Thus the background estimate remains

\[
 \left|q^A u_{\theta,{\rm bg}}(x_\tau,1-\tau)-e_0\right|
 \le C_{\rm bg}q^{2h}                                \tag{6}
\]

for a finite `C_bg` and sufficiently small `q`. Positive-order background
coefficients are small after normalization; they are not incorrectly set
to zero.

Second, Q153 applies the certified complete-potential diagonal to that same
physical family. The theorem retains the closed annular support

\[
                         X_a\le X\le X_b              \tag{7}
\]

for every stage potential, direct swirl, pressure difference and synchronized
stage-zero annular piece. Cutting a potential before curl, differentiating
it, adding the direct swirl and applying the scale cutoff do not enlarge its
closed `X` support. Since `X_in<X_a`, every physical diagonal increment is
identically zero in a neighborhood of (2). This statement uses support and
does not compare canonical and physical roots or frames.

Third, Q155 chooses the localization cutoff `c` equal to one on the late
protected core. There `grad c=0`, so

\[
 u=\operatorname{curl}(c\mathcal A)+c\mathcal B
   =\operatorname{curl}\mathcal A+\mathcal B=u_0.     \tag{8}
\]

The complete pressure gauge is unchanged throughout, although pressure is
not needed for the pointwise velocity lower bound.

Combining (4), (6), (7) and (8) gives, for the exact Q153--Q155 physical
field,

\[
 \tau^A u_\theta(x_\tau,1-\tau)=e_0+R_{\rm bg}(\tau),
 \qquad |R_{\rm bg}(\tau)|\le C_{\rm bg}\tau^{2h}.    \tag{9}
\]

## 3. One explicit positive threshold

Let `q0>0` be the range of (6), `tau0>0` the late time-cutoff width and
`r_L>0` the spatial plateau radius. Because `1/(2h)=100`, define

\[
 q_{\rm core}=\min\left\{q_0,
        \left({e_0\over2C_{\rm bg}}\right)^{100}\right\}, \tag{10}
\]

with the harmless replacement of the ratio by `min(1,ratio)` if it exceeds
one. The path lies in Q155's plateau when

\[
 \tau\le\min\left\{{\tau_0\over2},
                    {r_L^2\over8X_{\rm in}}\right\}.  \tag{11}
\]

Thus the combined threshold

\[
 \tau_*=\min\left\{q_{\rm core},{\tau_0\over2},
                    {r_L^2\over8X_{\rm in}}\right\}>0 \tag{12}
\]

gives the quantitative same-physical-field estimate

\[
 \boxed{
 |u_\theta(x_\tau,1-\tau)|
 \ge {e_0\over2}\tau^{-101/200}}
 \qquad(0<\tau\le\tau_*).                             \tag{13}
\]

No numerical lower bound for the enormous selected profile or Borel scales
is asserted. Positivity and finiteness in (5), (6) and (12) are analytical
consequences of the named conditional contracts.

## 4. Join to Q156 in the correct limit order

Let `v,P` be a locally smooth solution with Q155's same compact smooth force
and zero datum, and suppose

\[
                 \sup_{0\le t\le T}\|v(t)\|_2<\infty  \tag{14}
\]

for every fixed `T<1`. Q156 recovers `grad(P-p)` without pressure decay and
proves, by the expanding-ball estimate,

\[
                 v=u\quad\hbox{on }\mathbb R^3\times[0,1). \tag{15}
\]

The order is essential: fix `T`, recover pressure, send the comparison
radius to infinity, and only then let `tau` decrease along (2). Equations
(13)--(15) imply

\[
 |v_\theta(x_\tau,1-\tau)|
 \ge {e_0\over2}\tau^{-101/200}\longrightarrow\infty. \tag{16}
\]

A solution smooth through `(0,1)` is bounded on one compact neighborhood of
that point, contradicting (16). Therefore, relative to every retained
continuum premise, there is no globally smooth same-force, zero-datum
competitor in the fixed-interval finite-energy class.

This is a **conditional competitor-exclusion implication**. It is not an
unconditional admissible-force theorem, an unforced singularity result or a
verified Navier--Stokes blow-up. No future velocity or pressure is
constructed, and passing tests do not promote conditional providers to
continuum facts.

## 5. Next obligation

The next justified checkpoint is a current physical dependency and
quantifier audit for Q147--Q157. It must list every still-conditional
continuum primitive on this exact profile, physical family, diagonal,
pressure gauge and localization, then prove or precisely obstruct the first
unresolved item. Further downstream implication replay cannot make the
conditional exclusion unconditional.

Reproduce with fresh output paths:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q157_physical_core_exclusion/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q157_physical_core_exclusion/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q157_physical_core_exclusion/evidence.json \
  --output q157_replay.json
```
