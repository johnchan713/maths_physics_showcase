# Q90: first-real entrance compatibility and full normal forcing

Parent: `f47b6620304deb4b09e6c8c11b316968bac05a56` (Q89).
The first-real source is still at the conditional order B=1-kappa.
This checkpoint resolves its entrance and normal-pressure *interfaces*;
it does not infer an actual full-path source norm from Q89's unweighted
residual estimate. All original parameters, cutoffs, pressure gauges,
backgrounds, roots, support separations and accepted domain Dv remain fixed.
Historical evidence is unchanged.

## 1. A finite entrance-jet equivalence for the complete inverse

Freeze one label, lift and nonzero harmonic m; omega=km. On its complete
interval [0,L], use the original constrained system (full-path ODE
reconstruction, equations (1)--(3)):

    t'+Kt+m^2 d t+i omega n pi=-f,  n.t=0, t(0)=0,
    pi=i[(K^T n-n').t+n.f]/(omega |n|^2).                 (1)

Here n is smooth and |n| has the retained positive lower bound. Primes
are pulse-time derivatives; do not freeze the moving normal. For every
finite r>=0, assuming sufficient smooth coefficient/source jets, exactly

    f^(j)(0)=0 for 0<=j<=r
      iff t^(j)(0)=0 for 0<=j<=r+1
          AND pi^(j)(0)=0 for 0<=j<=r.                  (2)

Proof: eliminate pi to obtain t'=H t-Pi_perp f with
H=-K+n(n^T K-(n')^T)/|n|^2-m^2 d I. Repeated Leibniz gives

    t^(j+1)(0)=sum_i=0^j binom(j,i) H^(i)(0)t^(j-i)(0)
                         -(Pi_perp f)^(j)(0).           (3)

Vanishing f jets and t(0)=0 imply all the stated t jets; (1) implies the
pressure jets. Conversely differentiate the complete first equation (1)
through r. Every term vanishes under the right side of (2), hence so do
all f jets. This converse would be false with pressure omitted. The same
induction works for mixed slow/transverse parameter jets on an entrance
collar. Compatibility at one fixed parameter value does not supply the
needed parameter derivatives.

For zero extension of BOTH response and pressure through every finite
order, full f must therefore be entrance-flat. Tangent forcing alone
controls only the velocity: take constant n=e1, K=0 and f=n g with
nonzero entrance g. Then t=0, pi=i g/omega. All velocity entrance jets
vanish, but pressure cannot be continuously extended by zero. This is an
exact omission control, not a counterexample to the actual source.

## 2. Apply the support grammar to the first-real source

Use precisely f_R from Q89, excluding its additive temporal/base ledger.
Its regular terms are full linear wave remainders, both mean-wave
orientations and nonzero wave-wave products. The actual complete cut
velocities and wave pressures, not raw homogeneous fields, enter this
recipe. On each original temporal endpoint collar they are zero with all
jets under the retained cutoff/support contracts. Differential operators,
phase factoring and multiplication by smooth same-state means preserve
this property. The regular split writes its linear part with the temporal
cutoff; its separately retained commutator is not forwarded to the inverse.
The same argument applies if the cutoff is only flat at the endpoint
rather than identically zero on a collar, by Leibniz at the endpoint.

Original and signed waves must use the same retained endpoint extension
on the common lifted rectangle. Other labels contribute no overlapping
wave factor only by the original CLOSED enlarged support separation; a
single physical trace is insufficient. These are the hypotheses under
which the argument applies, not newly certified support inequalities.
Nonlocal means can fill radial gaps and retain pressure memory, but every
term of f_R has a wave factor. They cannot create an endpoint source where
that factor and all its jets vanish. Smoothness of those mean multipliers
on the whole collar remains required.

Since the original P is positive and smooth at each fixed-label endpoint,
f_R/P is also flat there. This statement uses no band-uniform lower bound
on endpoint P. Uniform magnitude estimates for its jets are a separate
issue in section 3. Thus the first-real *entrance compatibility* follows
conditionally from the actual cutoff grammar, without a new PDE estimate.

Flat source at the EXIT does not force a zero raw exit response. Duhamel
requires the additional vector moment

    t(L)=-int_0^L U(L,s) Pi_perp(s) f(s) ds=0            (4)

and its parameter jets for an uncut zero extension. For constant H=0 and
a positive interior smooth bump in a fixed tangent direction, the source
is flat at both ends but its integral is nonzero. A polynomial v^2(1-v)^2
is an exact finite-jet diagnostic of the same failure through order one.
The existing construction instead applies the original cutoff to complete
potential AND pressure, retains the full cutoff residual, and estimates
its Gaussian tail separately. No new outgoing moment condition is assumed.

## 3. The precise full-source inequality; no separate normal-source oracle

For each finite native request R, the remaining actual inequality is

    |D^I f_R(x,H)| <= epsilon^B sqrt(zeta) delta^(-p_R)
                           P(v(H),x) F_R(S),  |I|<=R,   (5)

on EVERY point of the original enlarged lifted rectangle and EVERY accepted
band, with finite band-uniform polynomial F_R, finite poles and all lifts
specified BEFORE taking the supremum. This includes the full unprojected
vector. Q26 supplies fixed-q domain closure, not its magnitude. Q89's
unweighted bound and a physical trace cannot imply (5).

Given (5), Q21's affine pullback costs only 2^(k+r) S_*^(-r) on a mixed
jet with k auxiliary and r path derivatives. Q20's triangular Gaussian
quotient conversion then bounds g=f_R/P without extra source derivatives,
reciprocal min(P), or a new smaller-q threshold. Both conversions require
the same original Gaussian and the full domain in (5). Their finite
coefficient costs are retained in the new bound G_j(S).

Put u=t/P, a=(K^T n-n')/|n|^2 and b=n/|n|^2. Exactly

    pi/P=i(a.u+b.g)/omega.                              (6)

If a_j,b_j,U_j,G_j bound all raw mixed coefficient jets of order j in
Euclidean/operator norms on that domain, Leibniz and Cauchy--Schwarz give

    |D^I(pi/P)| <= |omega|^-1 sum_i=0^j binom(j,i)
                         [a_i U_(j-i)+b_i G_(j-i)].     (7)

No extra source derivative is required for pressure through order j.
The a coefficient DOES require normal jets through j+1, K through j and
the existing inverse-normal jets; these are not free. For the source-normal
summand alone use b_i G_(j-i). Thus a proved FULL (5), plus the retained
normal coefficients, already controls normal forcing. A bound only on
Pi_perp f does not. The pressure's source-normal contribution has sufficient
order B+1/2 because 1/k<=epsilon^(1/2); this says nothing about the response
term until its own inverse bound is supplied.

## 4. Frontier and checks

Established analytically: (2), its conditional application to the first-real
cutoff grammar, and the explicit pressure coefficient envelope (7). The
unclosed numerical/continuum input is still (5) for the ACTUAL first-real
fields, with the requested whole-path background/propagator norms. Neither
finiteness of the source expression nor named epsilon orders certify it.
A full proof must derive its finite polynomial from the current complete
wave and mean coefficient providers, including spatial cutoffs, viscosity,
full wave pressure, axial velocity and Pi1/Pi2-dependent means. Only then
apply the first sourced inverse and recompute the changed covariance,
means and stage-specific compact tails. No unconditional pulse bound,
complete residual improvement, infinite realization or blow-up follows.

Exact manufactured tests check moving-normal pressure, the finite entrance
recurrence, a normal-only omission, a nonzero exit moment and mixed pressure
jet majorants. They are not selected-field continuum certificates. Run
`python test_audit.py -v` and `python audit.py --verify-record evidence.json
--output FRESH.json`. The Q89 evidence graph is replayed unchanged.
