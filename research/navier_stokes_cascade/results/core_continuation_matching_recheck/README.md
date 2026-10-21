# Independent continuation and same-profile matching-root recheck

Parent: `96865d42b7977ee08cc27ae886d6f07aff0acb8e`.
Status: `core-continuation-and-both-same-profile-matching-roots-independently-rechecked`.

This checkpoint completes the next obligation after the flat activation
collar.  It rederives the prescribed continuation from the unchanged core,
checks that the annular and compact five-row roots lie in the old uniqueness
balls, and traces one actual pre-modulation profile through the intervening
outer prefix.  Every pressure statement keeps `Pi_axis+Cp`; no interface
resets pressure memory.

The argument is relative to the exact core and the frozen existence/value
floors of the cited defining maps.  The rational inequalities below are
analytical estimates.  The executable exact-arithmetic controls detect
transcription, sign, scale, row, and scope errors, but they are not a formal
continuum proof.  This checkpoint does not yet apply compact modulation or
the exact repair to obtain a repaired cone estimate, and it proves no
unconditional pulse, residual, force, or blow-up conclusion.  Historical
evidence is not edited and the width remains `C^-20`.

## 1. Correct continuation after the reference taper

Let `y=log(X/Xa)` and `t1=C^-20`.  The reference is flattened by prescribing
`(log F_r)_y=0` after its taper.  Since `E_r=sqrt(2X)F_r`, this means

\[
 (\log E_r)_y=\tfrac12,\qquad a_r=1-2(\log E_r)_y=0.       \tag{1}
\]

It does **not** mean `p1,r=0`.  Thus the first-collar shortcut
`(log E_a-log E_r)_y=(1-kappa)p1,r/2` cannot be continued past the region
where the reference is still the native field.  The actual continuation is
instead defined from the original core trace by

\[
 \log E_a(y)=\log E_c(X_a)+\frac12\int_0^y(1-a(v))\,dv,
 \qquad U_a(y)=U_c(X_a)+\int_0^y b(v)\,dv.                 \tag{2}
\]

The four successive prescribed slopes are

\[
\begin{array}{ll}
 a=\kappa p_{1,r},&b=-\kappa Xn_{s,r}/2,\\
 a=\kappa_0p_{1,r},&b=-\beta\kappa_0Xn_{s,r}/2,\\
 a=(1-\theta)\kappa_0p_{1,r}+(4/5)\theta,&b=0,\\
 a=4/5,&b=0.
\end{array}                                                \tag{3}
\]

All cutoff endpoints are flat.  The last line gives
`(log E)_y=(1-a)/2=1/10`, exactly the slope used at the shape entry.  There
is no inverse `kappa` or inverse `p1,r` in (2).

The five actual forward densities in logarithmic radius are

\[
 (M,I,J,S,C_p)_y=
 (XU,\sqrt2X^{3/2}E,\sqrt2X^{3/2}UE,
 X(U^2-E^2/2),E^2/2),                                    \tag{4}
\]

with `Pi=Pi_axis+Cp`.  Equation (4) is integrated from the original five
core traces.  Reference moments are not reused after changing `E,U`, and the
eta-dependent axis datum is not discarded.  The finite product, primitive,
exponential, reciprocal, and Cauchy-strip estimates of the retained
continuation theorem then give every fixed mixed jet through `Xsep`.  Width
derivatives cost powers `C^(20r)`; they are finite at the selected parameters,
not uniform in `C` or in derivative order.

## 2. The first matching root is the same annular root

After the shape, extend the actual fields to `XR exp(-8)`, restore the axial
field on `-8<log(X/XR)<-7`, and use the existing five disjoint bumps on the
annular patch.  In the original normalized row order the exact map is

\[
 B_ac+Q_\eta(c,c)=z,qquad
 T_a(M,I,J,S,C_p)=(M,J-4\eta I,I,-(S-8\eta M),C_p).       \tag{5}
\]

The inverse of `T_a` is

\[
 (r_1,r_2,r_3,r_4,r_5)\mapsto
 (r_1,r_3,r_2+4\eta r_3,8\eta r_1-r_4,r_5),              \tag{6}
\]

so the transformed zero set is exactly the original five-row zero set.  The
U-square energy coefficient is eta dependent.  Its normalized derivative
list is `(4,8,8,4,1)` and vanishes only after order four; freezing this
coefficient would give the wrong root jets.

For `epsilon=10^-16`, the actual target has factorial `C2` norm
`delta<0.286 epsilon`.  The positive-mass estimates give
`beta=||B_a^-1||<=610` and `||Q_eta||<=50`.  Hence

\[
 8\beta^2q\delta=4.256824\,10^{-9}<1,
 \qquad 2\beta\delta=348.92\epsilon.                      \tag{7}
\]

The new ball lies strictly inside both the old `572 epsilon` factorial-`C2`
ball and the old `6000 epsilon` `C1` uniqueness ball.  It therefore identifies
the same accepted root, rather than choosing a new numerical branch.  Its
pointwise Jacobian inverse is at most `1220`; differentiation of the same
functional root retains all eta derivatives of `Q_eta`.

## 3. One profile reaches the second matching root

The radii have the strict order

\[
 X_{sep}<X_Re^{-8}<X_Re^{-7}<X_Re^{-6}<X_Re^{-5}<\rho_p
 <\rho_pe^5=X_{rep}.                                     \tag{8}
\]

After the annular root, the outer prefix consists of finitely many positive
exponentials of integrals of fixed smooth radial cutoffs, times smooth
functions of eta.  Flat endpoint jets agree, and positive value floors give
fixed reciprocal jets.  Consequently the same profile has finite mixed jets
on the whole compact input rectangle `[Xminus,Xrep]x[-1,1]`.

This trace keeps a frequently lost moment.  On the power interval `U=0`, but

\[
 M=\mu\eta,qquad
 \mu=4X_R\left[e+\int_1^{T+1}e^y\chi_A(y)\,dy\right]>0.   \tag{9}
\]

Zero axial velocity does not reset the accumulated axial moment.  Heat
compensation starts at `Xrep`, so no later heat-tail premise is needed to
supply the compact input.

On the first reserved patch the exact compact map is

\[
 B_cc+Q(c,c)=z_N,qquad
 T_c(M,I,J,S,C_p)=(M,(M-J)/\lambda,I,-S,C_p).             \tag{10}
\]

Because the selected `lambda=exp(-4T)` is strictly positive,

\[
 (r_1,r_2,r_3,r_4,r_5)\mapsto
 (r_1,r_3,r_1-\lambda r_2,-r_4,r_5)                      \tag{11}
\]

is the inverse transformation.  A `lambda=0` limit can bound weights but
cannot repair independent `M,J` data.  Both field corrections use the same
original amplitude `e(eta)`; only after this complete normalization are
`B_c,Q` eta independent.

The continuum bump bounds give `beta<=800`, `q<=25000`, and
`delta<=D/N` with `D=H^16`, `N=1+floor(H^32)>H^32`.  Thus

\[
 8\beta^2q\delta<128000000000H^{-16}
 <2^{37}H^{-16}<2^{-219},                                \tag{12}
\]

using `H>2^16`.  The same root has full sum-`C1` norm at most
`1600D/N<1600H^-16`, lies in the old `2000D/N` uniqueness ball, and has
pointwise inverse at most `1600`.  The `1/lambda` target loss remains; matrix
preconditioning does not erase it.

## 4. Complete pressure row and exact frontier

Both (5) and (10) contain `Cp` as the fifth independent row.  At each right
endpoint all five moment differences vanish as functions of eta, so all fixed
eta jets vanish too.  Intermediate pressure memory may be nonzero.  It is
carried forward and canceled only by the exact endpoint identities.

If the fifth row were omitted, a remaining pressure defect `p(eta)` would
contribute

\[
 4(1/2+h)\eta p-(1-\eta^2)p_\eta                         \tag{13}
\]

to the axial numerator (and the corresponding `X/L` factor to the stress).
For `p(0)=0,p_eta(0)=1`, (13) equals `-1`: a zero pressure value at one angle
does not remove the defect.

This checkpoint therefore fixes the prescribed continuation and the two
matching roots on one selected pre-modulation profile.  The next obligation
is to audit the signed variance root, moving phase, fixed-frequency compact
modulation, and application of this same exact repair, retaining the moment
gain and full pressure memory, and then to join its cone estimates to the
inner and annular/exterior bounds.  Global pressure/moment identities,
coefficient and reserved-mean supports, uniform physical pulses, covariance
means (including axial velocity and pressure), and the complete viscosity,
cutoff, and interaction residual remain separate obligations.

Reproduce from the repository root:

```bash
python3 research/navier_stokes_cascade/results/core_continuation_matching_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_continuation_matching_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_continuation_matching_recheck/evidence.json \
  --output core_continuation_matching_recheck_replay.json
```
