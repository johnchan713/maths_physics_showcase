# Same-profile global pressure, moment, and exterior-stress recheck

Parent: `fbefd07f4fb1b3738a00429dd282ce6f4e5c715e`.
Status: `same-repaired-profile-global-pressure-moments-and-exterior-stress-identities-independently-rechecked`.

The parent fixes one repaired leading profile and joins its inner, compact,
annular, and exterior value cones.  This checkpoint follows that same
profile through the global moment ledger.  It corrects a notation trap at
the interface between finite matching and the infinite heat tail, retains
the complete pressure memory, and rederives the exact exterior defect map.

The conclusions are conditional on the exact roots and common-tail facts
already frozen in the dependencies.  Exact-rational checks verify the
displayed algebra and negative controls; they are not a proof-assistant
certificate.  This checkpoint does not prove positive coefficient supports,
uniform physical-pulse estimates, complete covariance corrections, the full
Navier--Stokes residual, or blow-up.

## 1. The finite global heat object is renormalized

Write

\[
 A=\tfrac12+h,\quad D=\tfrac12-h,\quad d=1-\eta^2,
 \quad L=1-2h\eta^2,
\]

where the selected parameter has `0<h<1`.  In the terminal heat region,

\[
 E=c_\infty X^{-A}{\cal H}(2d/X),\qquad
 H=\sqrt{2X}E,\qquad H_{\rm pow}=kX^{-h},
 \quad k=\sqrt2c_\infty,
\]

and `k` is independent of `eta`.  The Gamma representation is

\[
 {\cal H}(Z)=\frac1{\Gamma(1+h)}
   \int_0^\infty e^{-v}v^h(1+Zv)^{-h}\,dv,              \tag{1}
\]

with exact derivatives

\[
 {\cal H}^{(j)}(Z)=\frac{(-1)^j(h)_j}{\Gamma(1+h)}
   \int_0^\infty e^{-v}v^{h+j}(1+Zv)^{-h-j}\,dv.        \tag{2}
\]

Taylor's formula at `Z=0`, with `Z=2d/X`, yields

\[
 H-H_{\rm pow}=O_h(X^{-1-h})                              \tag{3}
\]

and the same integrable order for every fixed `eta` derivative.  At the
axis, `H=2XF` with smooth `F`, while `Hpow=kX^-h` is integrable because
`h<1`.  Consequently

\[
 R(\eta)=\int_0^\infty(H-H_{\rm pow})\,dX
 =\lim_{X\to\infty}\left[I(X)-\frac{XH_{\rm pow}(X)}{1-h}\right]             \tag{4}
\]

is the finite global angular moment.  The raw cumulative `I(infinity)`
diverges and is never set to zero.  The terminal densities satisfy
`E^2/2=O(X^-1-2h)` for `S` and `E^2/(2X)=O(X^-2-2h)` for `Cp`, so their
tails converge as well.

This convergence is not uniform at `h=0`.  The pointwise heat edit tends to
zero as `h` decreases, but its infinite renormalized integral has a finite
angular limit.  Thus setting `h=0` before integrating is not an allowed
step.

## 2. Correction: endpoint `I` is not the global row `R`

Both matching maps act on the finite cumulative row

\[
 (M,I,J,S,C_p),                                           \tag{5}
\]

not on `(M,R,J,S,Cp)`.  The annular transform and inverse are

\[
 T_a(M,I,J,S,C_p)=(M,J-4\eta I,I,-(S-8\eta M),C_p),       \tag{6}
\]

\[
 T_a^{-1}(r)=(r_1,r_3,r_2+4\eta r_3,8\eta r_1-r_4,r_5). \tag{7}
\]

For positive `lambda`, the compact transform and inverse are

\[
 T_c(M,I,J,S,C_p)=(M,(M-J)/\lambda,I,-S,C_p),             \tag{8}
\]

\[
 T_c^{-1}(r)=(r_1,r_3,r_1-\lambda r_2,-r_4,r_5).         \tag{9}
\]

For two profiles, splitting (4) at the repair endpoint gives the exact
conversion formula

\[
 \Delta R=\Delta I(X_{\rm rep})+
   \int_{X_{\rm rep}}^\infty\Delta(H-H_{\rm pow})\,dX.  \tag{10}
\]

The same repaired profile and its comparison target have the same future
heat tail and the same power subtraction.  Therefore the exact endpoint
match `Delta I(Xrep)=0` implies `Delta R=0`.  Without both common-tail facts,
the implication would be false.  This is the correction recorded here;
historical evidence is left unchanged.

## 3. All five global defects survive the edit ledger

Use global defect order `(M,R,J,S,Cp)`, with pressure defect
`Pi_axis+Cp`, not `Cp` alone.  On the fixed profile:

| Exact operation | Global consequence |
|---|---|
| Axial and strength roots | Set `M=J=S=0`; leave `I,Cp` unchanged |
| Angular roots and waiting interval | Set the terminal `R=0`; preserve total `Cp`; include their energy in the strength equation |
| Heat edit and its separate compensation | Cancel the signed changes in `R,S,Cp`; leave `M,J` unchanged |
| Shape and annular root | Restore all five finite endpoint defects |
| Modulation and compact repair | Restore all five finite endpoint defects again |

The last two rows give global restoration only through (10).  Every row is
an equality of functions of `eta`, including its fixed derivatives, rather
than a finite set of angular samples.  Hence the completed leading profile
has

\[
 M_\infty=R=J_\infty=S_\infty=0,
 \qquad \Pi_{\rm axis}+C_{p,\infty}=0.                   \tag{11}
\]

Keeping the original gauge gives the backward pressure formula

\[
 \Pi(X)=-\int_X^\infty\frac{E(x)^2}{2x}\,dx.             \tag{12}
\]

The heat edit, compensation, shape repair, and compact loop may each carry
nonzero partial pressure memory between their supports.  Equation (12)
follows from the total exact cancellation; no intermediate pressure is
reset.

## 4. Exact exterior defect join

Compare two profiles with the same terminal local fields but cumulative
defects `(m,R,j,s,p)`, where primes denote `eta` derivatives and `p` is the
pressure defect.  Put

\[
 B_m=2D\eta m+d m'.                                      \tag{13}
\]

Direct subtraction of the profile identities gives

\[
 \Delta V_0=-\frac{B_m}{L},                              \tag{14}
\]

\[
 \Delta(r^2T_\theta)=
 \frac{HB_m+(1-h)R-D\eta R'-dj'+2(h-D)\eta j}{L},        \tag{15}
\]

\[
 \Delta(rT_z)=
 \frac{D(m-\eta m')+4h\eta s-ds'
       +X\{4A\eta p-dp'\}}{L}.                          \tag{16}
\]

These are exact identities, not large-`X` estimates.  All five functional
defects in (11), together with their angular derivatives, therefore remove
the radial-velocity and both weighted-stress memories.  The pressure term
and its derivative are essential: at `eta=0`, a defect with `p(0)=0` and
`p'(0)=1` contributes `-X` to (16).  Likewise `U=0` does not erase a
nonzero `m'`; (14) detects it immediately.

Thus the same-gauge backward stresses join the inner/compact/annular value
construction to the terminal heat exterior without discarding pressure
memory.  This closes the global leading moment/stress obligation, but says
nothing yet about positive-order supports or a complete residual estimate.

## 5. Exact frontier

The next checkpoint must independently audit the finite coefficient solves,
their common supports, and the untouched reserved velocity and exceptional
pressure means.  After that remain uniform actual-pulse existence,
inverse/gap/derivative/support estimates at every physical scale; complete
wave covariance and mean corrections including radial and axial velocity,
curl remainders, and both pressure tails; and the full residual including
viscosity, every cutoff, and every interaction.  No unconditional pulse
bound, residual improvement, force, or blow-up conclusion is made.

Reproduce from the repository root:

```bash
python3 research/navier_stokes_cascade/results/core_global_pressure_moment_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_global_pressure_moment_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_global_pressure_moment_recheck/evidence.json \
  --output core_global_pressure_moment_recheck_replay.json
```
