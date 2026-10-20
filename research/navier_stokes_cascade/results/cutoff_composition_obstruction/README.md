# Cutoff composition survives the moving angle and zero-mean integration

Parent: `117160fffb9824fb26831b72ab8ce03638eba082`.
Status: `local-loop-composition-obstruction-not-an-actual-profile-counterexample`.

This checkpoint closes an alternative explanation for the open variance-cutoff
budget: reparametrization and periodic integration do NOT automatically remove
that budget. A smooth local input family with uniformly small native angular
C4 norms has arbitrarily large first angular derivatives of both the variance
target and the actual moving-angle, zero-mean loop primitive. Thus native
envelopes plus cutoff range/plateaus cannot imply the parent's condition (S),
or even a parameter-only C1 primitive bound. This is an analytical obstruction
to that inference, NOT a counterexample to the retained global profile.

No actual cutoff is replaced. No global pressure-compatible realization of the
test family is asserted. Full pressure/axis/Cp/axial memory and every earlier
root, support, width C^-20, frequency and frozen record remain unchanged.
The parent's conditional fourth-angular implication remains valid. Actual
shape/composition acceptance, radial/state cone joining, FINAL global jets,
physical pulses, complete covariance/means/residual and blow-up remain open.

Source: retained [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
printed pp.159--160, (C.4)--(C.10) and its fixed-origin circle lift.
PDF SHA256: `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its theorem is not used as a premise.

## 1. Uniformly bounded local inputs, not a new physical construction

Fix ANY 0<delta<=1, including the unchanged deltaL=exp(-A^20), and integer
n>=1. Let sigma be the explicitly defined squared-flat step used only to
construct AUXILIARY examples: on (0,1),

`sigma(x)=exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2))`,

extended by 0 and 1. In particular sigma(1/2)=1/2, sigma'(1/2)=8.
Define a permitted variance cutoff by

`zeta_n(v)=1-sigma(n*(8*(v-2)/delta-3/2)+1/2)`.

It is smooth, nonnegative, monotone, <=1, equals one for
v<=2+delta/8 and zero for v>=2+delta/4. On eta in [-1,1] set

`v_s=2+3delta/16+(delta/8)*eta`, `a=E=1`,
`b_s=sqrt(v_s-1)`, `t_s=-b_s`, `p=p2=0`.

Here a+b_s^2/a=v_s exactly. The normalized C4 norms of v_s,b_s,t_s
are <3, independently of n: v_s>=2, its slope is <=1/8, and
|binom(1/2,j)|<=1 controls the square-root derivatives. The other inputs
are constants. These obey the generic native envelopes of the preceding
loop estimate for A>=2^16. This tests only that scalar/loop input interface;
p2=0 is an auxiliary scalar slice, NOT a reset or computation of the actual
global pressure. It does not supply radial matching or full cone premises.

Let D=2+delta/2-v_s, rho=zeta_n(v_s)^2 D and s=zeta_n(v_s)*sqrt(D),
extended by zero. D>=delta/4 wherever the cutoff is nonzero. At eta=0,
write x=2+3delta/16, r=5delta/64, b=sqrt(x-1), V=x+r. Then

`s=sqrt(5delta)/8`,
`s'=-sqrt(delta)*(2sqrt(5)*n+1/(8sqrt(5)))`.

Consequently `|s'|>=2sqrt(5)*n*sqrt(delta)`. For each fixed positive
delta and each finite proposed bound this diverges as n increases. In the
retained hierarchy, (S) necessarily fails for this AUXILIARY family when
`log n>A^32+A^20/2-log(2sqrt(5))`. No such n is assigned to the actual cutoff.

## 2. Exact moving primitive: integration does not rescue a uniform bound

At p2=0 the regular signed variance branch gives d0*mu=sqrt(2)*s and

`t(theta)=t_s+sqrt(2)*s*sin(theta)`, `v=v_s+s^2`,
`w=(1+t^2)/(2pi*v)`, `Phi(theta)=integral_0^theta w`,
`gA(phi)=-(v/(1+t(Theta(phi))^2)-1)/2`.

The lift is strictly increasing and Phi(2pi)=1, exactly. Let P(phi) be
the primitive of gA anchored at phi=0. Direct integration, before inverting
the angle, gives

`P(Phi(theta))=[2*t_s*sqrt(2)*s*(1-cos(theta))-(s^2/2)*sin(2theta)]/(4pi*v)`.

Equivalently `P(Phi(theta))=(Phi(theta)-theta/(2pi))/2`.
The zero-mean primitive is Acal=P-integral_0^1 P; its angular derivative
has the SAME phase oscillation as P_eta. Choose the fixed external phase
phi0=Phi(eta=0,theta=pi). Differentiation at fixed phi0, including
Theta_eta=-Phi_eta/w, yields

`P_eta(phi0)=sqrt(2)/(pi*x)*[t_s'*s+t_s*s'-t_s*s*v'/V]`.

At the stated base point the bracket equals

`b*sqrt(delta)*(2sqrt(5)*n+1/(8sqrt(5)))*(x-r)/V`
`+7*delta^2*s/(1024*b*V)`.

Both terms are positive. Using b>=1, x<3 and (x-r)/V>1/2 gives

`P_eta(phi0)>sqrt(10)/(3pi)*n*sqrt(delta)`.

Since P_eta(0)=0, subtracting the phase mean cannot cancel this difference:

`sup_phi |Acal_eta(0,phi)|>sqrt(10)/(6pi)*n*sqrt(delta)`.

Thus even the integrated, mean-zero moving primitive has no uniform C1
bound over the permitted cutoff class. For instance its H bound fails in
this family if `log n>log H+A^20/2+log(6pi/sqrt(10))` at deltaL.
This is NOT evidence that the actual primitive exceeds H. A frozen-angle
calculation would give P_eta smaller by the factor x/V: the exact moving
derivative is (V/x) times the frozen one, not the frozen derivative.

## 3. What an actual composition proof must still control

Put h=sqrt((v_star-v_s)/a) on the cutoff support. Raw fourth angular
derivatives satisfy the following identity, with v_j=D_eta^j v_s,
h_j=D_eta^j h and zeta_j the jth derivative in its SCALAR argument:

`s''''=h*zeta4*v1^4`
`+zeta3*(6*h*v1^2*v2+4*h1*v1^3)`
`+zeta2*(h*(3*v2^2+4*v1*v3)+12*h1*v1*v2+6*h2*v1^2)`
`+zeta1*(h*v4+4*h1*v3+6*h2*v2+4*h3*v1)+zeta*h4`.

The highest shape derivative has coefficient h*v1^4, strictly positive
at every nondegenerate angular crossing. This identity does not rule out
special cancellations or degeneracy in the ACTUAL profile; those would
need proof. Smooth flat endpoints alone do not provide them. One must
either enclose the unchanged shape budget Z in the parent's sufficient
inequality, or bound the complete actual composition, including these
terms, to prove |s|Ceta4<=G. The present example rules out obtaining that
bound from the generic native C4 envelopes alone.

## 4. Reproduction and limits

`composition.py` verifies exact symbolic center jets, the mean-one lift,
anchored primitive, moving derivative, positive remainder, and all fourth
Leibniz/Faa di Bruno terms. Manufactured 80-digit diagnostics independently
invert a quadrature-defined lift at nearby eta, integrate gA*w, and compare
its fixed-external-phase derivative with the exact formula. They also
reject the frozen-phase value. These are floating, non-interval diagnostics,
NOT evaluations of the actual enormous hierarchy or continuum proofs.

Run from this directory with the pinned requirements installed:

```sh
python -m unittest -v test_audit.py
python audit.py --verify-record evidence.json --output /tmp/cutoff-composition-fresh.json
```

Use a fresh output path. Evidence verifies historical source bindings and
hashes every existing frozen record without rewriting it. Scope-promotion,
missing-chain-term, frozen-phase and evidence-mutation checks are included.
Next: identify/enclose the actual cutoff or prove actual crossing/composition
control; then join radial/state and FINAL global derivative gates before
advancing unconditional all-scale pulses, means and full residual estimates.
