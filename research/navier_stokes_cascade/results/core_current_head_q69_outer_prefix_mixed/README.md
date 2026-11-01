# Q69: mixed derivatives on the native pre-modulation outer prefix

On [R,Xrep], R=X_R, the broader representation c f^theta used in Q68
specializes to E=c(y) f(eta), U=k(y) eta, y=log(X/R), f=1/(1+eta²).
The angular interpolation is later than this prefix. Relative to the same
exact annular match and frozen Md64 schedule, the following bounds have no
unknown derivatives of a ramp, field, or incoming moment on this interval.
They do not describe the modulated field.

## Four original stages in log radius

Set t=y minus the left endpoint in the table, b=(log c)' and T=exp(64)+10.
The native slopes are those of E, i.e. subtract 1/2 from the H slope.
See the [first and axial stages](../axial_stress_audit/diagnostics.py)
(rhs in independent_ode), the [native taper](../actual_cutoff_pullback/README.md),
the [intermediate stages](../intermediate_decay_audit/README.md) and the
[regional endpoints](../regional_cone_join/README.md).

| y interval | b | k |
|---|---|---|
| [0,1] | 1/10-(3/5)sigma(t) | 4 |
| [1,T+1] | -1/2 | 4(1-sigma(log(1+t)/64)) |
| [T+1,T+2] | -1/2-lambda sigma(t) | 0 |
| [T+2,log(Xrep/R)] | -1/2-lambda | 0 |

Here lambda=exp(-4T)<1/2. On the axial stage k is already zero
for t>=exp(64)-1. Thus all seams match to every order by the original
flat plateaus. Neither T nor any physical width is an adjustable
high-order parameter. In particular a fictitious linear cutoff at width
T would give the wrong derivative even at an interior point.

## Uniform finite-order bounds, including the logarithmic taper

Use the [Q63 bound for the actual squared-flat step](../core_current_head_q63_repair_anisotropy/README.md):
B_j=162·3260^j(j!)³ bounds its raw j-th derivative. On a linear ramp
of log width w the corresponding cost is B_j/w^j; BOTH slope ramps
above have w=1. Write [n,j] for unsigned Stirling numbers of the first
kind, with [0,0]=1. For n>=1 on the axial taper,

    |k^(n)(t)| <= 4 (1+t)^(-n) sum_(j=1)^n [n,j] B_j/64^j.   (1)

Indeed derivative n of F(log v) is v^(-n) times the falling-factorial
operator D(D-1)...(D-n+1)F at log v, D=d/d(log v).
Apply this with v=1+t and F(s)=sigma(s/64); the constant drops out.
Taking (1+t)^(-n)<=1 gives the uniform raw bound K_n, with K_0=4.
It applies on all four stages (k is constant elsewhere). This estimate
requires no analyticity of sigma at a flat endpoint.

Define H_1=1 and H_j=B_(j-1)/j! for j>=2. These bound the normalized
positive-order y derivatives of log c: |b|<=1 and the slope jumps of
size 3/5 and lambda are bounded by one. Set

    D_0=1,  D_n=(1/n) sum_(j=1)^n j H_j D_(n-j).             (2)

Applying c'=bc and (1/c)'=-b/c to normalized derivatives at each
point shows |c^(n)|/n! <= c D_n, and likewise for 1/c with its own
base value. This finite recurrence (the exponential series recurrence)
uses the actual value at the center; it does not exponentiate a norm
of high derivatives. For B in [R,Xrep] and Y=log(B/R), c(0)=P gives
c<=P exp(Y) and 1/c<=exp(Y)/P on [R,B], since |b|<=1.

For all i,j>=0 set S_i(V)=R^(-i)/i! sum_(r=0)^i [i,r] V_r.
The same falling-factorial identity now changes y to physical X.
On eta disks of radius 1/2 around [-1,1], |f|<=2, |eta|<=3/2
as in Q68. Cauchy's estimate in eta, after real radial differentiation,
gives the normalized mixed tables on [R,B]:

    |partial_X^i partial_eta^j E|/(i!j!)
       <= 2^(j+1) P exp(Y) S_i((r! D_r)_r),
    |partial_X^i partial_eta^j U|/(i!j!)
       <= v_j S_i(K),   v_0=1, v_1=1, v_j=0 for j>=2.      (3)

The U bound uses the exact real affine factor, not its disk supremum.
The real floor E>=P exp(-Y)/2 supplies the reciprocal in Q67.
All estimates include the flat endpoints; no radial Cauchy radius
is asserted. The sums specify symbolic constants without constructing
astronomically large stage lengths or evaluating an actual pulse.

## Feed the mixed fields into the original forward map

Use (3) in [Q67](../core_current_head_q67_premod_memory/README.md) for
each requested finite index set. For an incoming pure eta row of each
of M,I,J,S,Pi use [Q68](../core_current_head_q68_outer_prefix_traces/README.md)
equation (2), multiplied by 2^j at the actual cut. Every positive radial
row of a primitive uses its normalized local density row (i-1,j) divided
by i. For state output through total m and angular l the velocity inputs
through m+1,l+1 and incoming eta rows through l+1 are therefore explicit.
In addition to a,b_s,p2 this bounds p1=(-w+N_Q/H)/L using Q67's w,
H=sqrt(2X)E and N_Q=(1-h)I-Dphys eta I_eta-(1-eta²)J_eta+
2(h-Dphys)eta J, Dphys=1/2-h. This costs only the same first angular
derivative of a primitive. In p2 retain X(4Aphys eta Pi-(1-eta²)Pi_eta)
with the full Q68 trace. A zero local U does not remove that term.

This discharges the native mixed-input obligation on this outer prefix
only. It makes no assertion that the resulting large constants meet a
physical smallness threshold. The earlier core continuation, shape and
annular joining regions still need compatible quantitative mixed bounds;
then the normalized repair defect and physical inverse/pulse bounds must
be checked on the same bands. The axial and radial tail estimates,
complete residual improvement and blow-up are not established here.

Run test_audit.py and audit.py --verify-record evidence.json --output a fresh
path. The exact controls verify calculus and indexing, not continuum existence.
