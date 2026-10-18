# Explicit Gaussian curvature, primary Haar mass and cutoff-tail constants

Parent: `a70298fb0db13af479b2547b52e5bf1454b99e6d`.
Status: `current-fixed-tilt-gaussian-haar-and-tail-constants-from-spectral-floors`.
The parent CI run was still progressing without a failed step when inspected
on 18 October 2026. Historical checkpoints and frozen evidence are unchanged.

The current finite residual iteration already exists conditionally. This
checkpoint supplies a quantitative input it still uses: the Gaussian
envelope, the positive primary radial amplitude, its actual Haar mass and
the all-band cutoff-tail rate. These constants are derived from the SAME
fixed tilt and current leading spectral floors, not fitted to sampled
pulses. No pulse, positive root, width C^-20, radius, frequency, pressure
gauge, support or starting-band choice is changed.

The argument is relative to the reviewed core/profile and actual moving-
frame comparison estimates. It removes an unspecified Gaussian-constant
step; it does not discharge every core/finite-map dependency or establish
an unconditional complete residual or blow-up theorem. Exact rational and
numerical manufactured controls below are separate from the analytical
inequalities. No actual giant selected pulse is numerically evaluated.

## 1. Curvature of the SAME growing envelope

Write L=Ls, epsilon=Q^h, S=ell^2 and k=ceil(epsilon^-1/2). The original
tilt u is fixed and strictly positive. On the whole closed pulse set
`a=1/2+v/L`, `s_sigma=sigma*u*a`, 0<=v<=L. The reference rate is

\[
 r(v)=\lambda_0\left[(1+u^2a^2)^{-1/2}
                    -{1+u^2a^2\over(1+u^2)^{3/2}}\right],
 \qquad P'/P=r,\quad P(L/2)=1.                         \tag{1}
\]

Both signs have this SAME rate. The current leading/frame foundation
provides positive uniform bounds lambda_-<=lambda0<=lambda_+ on the
enlarged charts. For example, its leading lambda_star/lambda^upper bounds
and 1/4<=q/Q<=4 permit the conservative choices
`lambda_-=lambda_star/16`, `lambda_+=16 lambda^upper`, since h<.01.
These are finite parameter-dependent bounds, not numerical evaluations.

Differentiate (1), keeping the whole interval a in [1/2,3/2]:

\[
 -Lr'(v)=\lambda_0u^2a
 \left[(1+u^2a^2)^{-3/2}+{2\over(1+u^2)^{3/2}}\right]. \tag{2}
\]

For x>=1, x^-2<=x^-3/2<=x^-1. Thus the explicit positive constants

\[
 c={\lambda_-u^2\over2}
       [(1+9u^2/4)^{-2}+2(1+u^2)^{-2}],
\quad C={3\lambda_+u^2\over2}
       [(1+u^2/4)^{-1}+2(1+u^2)^{-1}]                  \tag{3}
\]

satisfy c<=-Lr'<=C on the WHOLE pulse. Since r(L/2)=0, integrating twice
on either side of the center proves

\[
 e^{-C(v-L/2)^2/(2L)}\le P(v)
 \le e^{-c(v-L/2)^2/(2L)}\le1.                        \tag{4}
\]

This proves the quantifiers in the Gaussian estimate; checking a grid of
v values would not. An independent closed expression useful for controls is

\[
 \log P=\lambda_0 L\left[
 {\operatorname{asinh}(ua)-\operatorname{asinh}(u)\over u}
 -{(a-1)+u^2(a^3-1)/3\over(1+u^2)^{3/2}}\right].    \tag{5}
\]

There is a genuine fixed-tilt requirement. Equation (2) is at most
`(9/2)lambda_+ u^2`, so any uniform positive lower curvature tends to
zero as u tends to zero. At u=0, r=0 and P=1. No Gaussian curvature or
Haar-mass constant uniform down to zero tilt follows from spectral
positivity alone. The original construction keeps its positive tilt.

## 2. Explicit amplitude factors from the actual frame comparison

Retain the reviewed frame estimates, in Euclidean operator norm,
`||E||<=CE/L`, `|d-dref|<=Cd/L`, with CE>0, Cd>=0. Put
`mu=lambda_-/sqrt(1+9u^2/4)` and require the ORIGINAL Riccati acceptance
condition L>=8CE/mu. Then rho=CE/(mu L)<=1/8.

For the real primary frame solution, zplus(0)=P(0), zminus(0)=0. The
ratio q=zminus/zplus obeys the reviewed outward barrier. Its boundary
derivative is bounded above by `-47CE/(64L)<0`; hence |q|<=rho. The
logarithmic growing coordinate satisfies
`|(log(zplus/P))'| <= [Cd+CE(1+rho)]/L`.
Consequently, with D=Cd+9CE/8 and radial amplitude x=zplus+zminus,

\[
 x_-P\le x\le x_+P,\qquad
 x_-={7\over8}e^{-D},\quad x_+={9\over8}e^D.         \tag{6}
\]

The scalar exponential stays positive, closing continuation of the ratio
argument. This is a lower bound for the uncut radial ODE amplitude. The
complete oscillating/cutoff velocity can vanish. Full moving-normal and
source-normal pressure are retained in the underlying projected system:
`pi_m=i(n.Kbg*t_m-nprime.t_m+n.f_m)/(km|n|^2)`.
No angular-component inverse or zero pressure is assumed.

## 3. Positive Haar mass and direction concentration, with radius retained

Use the SAME transverse cutoff with
`(r0/2)e^(-2/3)<=int chi^2<=r0`, the SAME temporal cutoff equal to one
for |v/L-1/2|<=1/5, and the unchanged positive Haar Jacobian J. The
manuscript cutoff support is inside |v/L-1/2|<1/3; the older geometric
replay's 2/5 bound is a safe enlargement, not a new cutoff.
For L>=25 the interval |v-L/2|<=sqrt(L) lies in its plateau. With
`ci=2r0/L`, the actual primary mass is

\[
 h_\sigma={J\over2}\left(\int\chi^2\right)c_i
                  \int\psi^2x^2dv
 ={r_0^2\over\sqrt L}j_\sigma,
\quad j_-\le j_\sigma\le j_+,                       \tag{7}
\]
\[
 j_-=J_{min}e^{-2/3}x_-^2e^{-C}>0,\qquad
 j_+=J_{max}x_+^2\sqrt{\pi/c}.                       \tag{8}
\]

For the lower bound use a length 2sqrt(L) interval and (4),(6). For the
upper bound integrate the Gaussian square on the whole line. No angular
Haar factor 1/2 or transverse area is discarded. In particular h is of
order r0^2/sqrt(L), hence r0^(3/2)S^-1/2 for L comparable to r0 S.
The roots and inverse retain their existing radius costs.

The same argument gives the useful explicit concentration bound

\[
 {\int |v-L/2|\psi^2x^2dv\over L\int\psi^2x^2dv}
 \le {x_+^2 e^C\over2x_-^2c\sqrt L}.                 \tag{9}
\]

Indeed the numerator is at most x_+^2 L/c. If the reviewed column
direction satisfies `|D(v)-Dcenter|<=D1|v-L/2|/L+D0/L`, its normalized
Haar column error is at most Cdir/sqrt(L), with
`Cdir=D0+D1*x_+^2*exp(C)/(2*x_-^2*c)` for L>=1. Insert this Cdir in
the existing primary-cone acceptance inequality
`L >= [4*Cdir*max(Amin^-1,u^-1)/gamma]^2`.
This identifies a finite value constant in that budget; it does not
replace the target or claim that a sampled/numerically specified first
band passes it. CE,Cd,D0,D1 and the original acceptance inequalities
remain the reviewed frame/column inputs. No new jet-dependent band is used.

The positive quantities can be far below floating-point range. Keep their
logarithms, for example
`log(j_-)=log(Jmin)-2/3+2log(7/8)-2D-C`.
Then the primary inverse estimate is the existing
`||Hpulse^-1||<=3*sqrt(L)/(r0^2*j_-*min(Amin,u))`.
Floating underflow of j_- does not make the analytical inverse singular.

## 4. A concrete all-band rate for every fixed cutoff-tail jet

Every nonconstant temporal cutoff derivative, and (1-psi) on the pulse,
lies where |v-L/2|>=L/5. Equation (4) and L>=2r0 ell^2 give

\[
 P\le e^{-cL/50}\le e^{-a_0\ell^2},\qquad
                         a_0=cr_0/25>0.             \tag{10}
\]

For a fixed jet, let the reviewed coefficient/physical bounds contribute
`K ell^B Q^-H`, Q=2^-ell. To obtain a desired Q^N bound use, for ell>=1,
`log ell<=ell` and complete the square:

\[
 \ell^B Q^{-H}e^{-a_0\ell^2}
 \le \exp\!\left({[B+(H+N)\log2]_+^2\over4a_0}\right)Q^N.
                                                               \tag{11}
\]

Thus every fixed tail jet is flat on the SAME entire accepted band
domain. The multiplicative constant depends on that jet, N and the frozen
parameters. No derivative-dependent starting band or numerical huge band
is necessary. Both the complete potential/pressure cutoff commutator and
the source tail retain these bounds. They stay exact additive residuals;
the full-path Gaussian hypothesis still governs every forward inverse.
Actual cutoff/A0/A1 velocities stay in nonlinear products. The complete
background Pi1/Pi2 pressure memory, axial velocity and moving pressure
moment are not deleted by this envelope calculation.

## 5. Evidence scope and next obligation

`gaussian_bounds.py` checks rational curvature enclosures, barrier budgets,
physical rescaling and tail-square completion. `gaussian_checks.py` uses
independent quadrature of (1), the closed formula (5), and noncommuting
normalized frame ODEs. It integrates the actual manufactured cutoff/Haar
mass rather than replacing it by envelope mass. These finite controls
cannot certify continuum suprema or the enormous selected construction.
The all-interval analytical proofs are (2)--(11).

The audit binds the frozen dependency sources. Separate fresh replays of
the relevant pulse foundations accompany publication. Next discharge the
remaining actual frame/column comparison constants and join the current
core/root/inverse bounds, full covariance, fresh mean pressure and every
viscous/cutoff/interaction term into one dependency-explicit residual
estimate. Infinite summability and every unconditional construction or
blow-up conclusion remain outside this checkpoint.

```sh
python3 research/navier_stokes_cascade/results/gaussian_mass_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/gaussian_mass_foundation/audit.py \
 --verify-record research/navier_stokes_cascade/results/gaussian_mass_foundation/evidence.json \
 --output /tmp/gaussian_mass_fresh.json
```
