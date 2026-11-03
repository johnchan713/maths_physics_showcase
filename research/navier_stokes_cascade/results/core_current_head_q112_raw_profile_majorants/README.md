# Q112: closed raw radial-power and moving-cutoff majorants

Parent: `ec0141cbcf42074bb55c227eada4129987aa0f30` (Q111).
On 3 November 2026 its research CI 41549497447 and field CI 41549497511
were still running; Q110 field CI 41542142363 had passed. The inspected
remote head was Q111. Q44 and Q46 already give finite power/cutoff jet
recurrences. This checkpoint proves closed all-order bounds dominating
those formulas and binds the RAW e=0 profile slots in Q110/Q111. No new
cutoff is selected and no actual A0 gain is certified. Historical evidence,
fields, root, accepted cover and width C_setup^-20 are unchanged.

## 1. The high-order input profile is an explicit power, not a new norm callback

The original raw compact kernel on the positive U hull is

    A_e f=k_e chi_U int_all f_e(U+a) f(U+a,y+Ma vr) da,
    f_e=R^(e-d+1)/d, k_e=d R^(d-1-e), U=R^d.             (1)

The original radial operator is e=0. In particular

    f0=R^(1-d)/d=U^(1/d-1)/d,
    k0=d R^(d-1)=d U^(1-1/d), g0=1.                    (2)

The labels d,M,c are fixed under coefficient differentiation. This does
not freeze the slow Z/T motion of chi. Use the SAME current positive
hull a<=R<=b and 1/3<d<5/2. Q44/Q111's rational enclosure is
u=min(1,a^3)<=U. Put H2=max(1,a^-2,b^2). Then for every N>=0,

    N_N(f0), N_N(k0) <= 3 H2 N! (3/u)^N.                (3)

Here N_N is the raw through-order seminorm, not a normalized Taylor
coefficient. To prove (3), for a general fixed power R^beta=U^alpha,

    partial_U^j(R^beta)=(alpha)_j R^beta U^-j,
    |(alpha)_j| <= j! (1+|alpha|)^j.                    (4)

The last inequality follows term by term from
|alpha-l|<=|alpha|+l<=(1+|alpha|)(l+1). For the two powers in (2),
|beta|<2, |alpha|<2, R^beta<=H2, and each leading scalar is <=3.
Because u<=1, j! (3/u)^j<=N! (3/u)^N for j<=N. This proves a raw
through-N bound on the whole fixed hull. The other coordinate derivatives
of these powers vanish. The e=0 profiles are not substituted for e=1/2.
No phase, Gaussian, mean, A1 velocity, local pressure or compensator is
removed by this multiplier estimate.

## 2. Retain the physical cutoff's slow motion

Q46 composes the ORIGINAL one-dimensional X profile phi with the moving
physical coordinate X(U,Z,T), X<=Xb. Let
sigma_j bound |phi^(j)|/j! on that same selected physical profile. Those
values remain original native inputs; they are not set to one or inferred
from a finite test. The normalized total-degree coordinate majorant has

    L_j=2/(d j u^j)+ell_j, j>=1,
    ell_j=(1/3) sum_(i=0)^j 32^i 64^(j-i),
    X(t)=Xb exp(L(t)), L(t)=sum_(j>=1) L_j t^j.          (5)

These are precisely the old Q46/Q27 coordinate-jet bounds. Their native
coordinate hypotheses remain. In particular ell_j retains the original
slow scale-motion derivatives. Dropping it would not bound the moving
cutoff. The compact profile has no auxiliary-y dependence.

Since d>1/3 and ell_j<= (2/3)64^j, for t below both radii,

    L(t) <=6(t/u)/(1-t/u)+(2/3)64t/(1-64t).

Let H=1+Xb and tau=u/(512H). This tau is only a generating-function
evaluation point, not a new physical scale, support or epsilon threshold. Then

    first term <=6/(512H-1)<1/(64H),
    second term <=2/(21H),
    L(tau)<1/(8H).

The exponential series implies exp(L)-1<=2L for 0<=L<=1/2. Thus
Xb[exp(L(tau))-1]<Xb/(4H)<1/4. For any FINITE order N the normalized
cutoff coefficients through N are dominated by those of

    C_N(t)=sum_(k=0)^N sigma_k [Xb(exp(L(t))-1)]^k.

All coefficients are nonnegative. Therefore C_N(tau)<=sum_(k=0)^N sigma_k.
For j<=N, the raw derivative bound is j! times its normalized coefficient;
tau<=1 and j!<=N! give

    N_N(chi) <= N! [512(1+Xb)/u]^N sum_(k=0)^N sigma_k.   (6)

Convergence is used only for the known coordinate majorants, not for
phi. A compact flat C-infinity cutoff need NOT be analytic and no analytic
or Gevrey assumption on its shape jets is introduced. The finite sum in
(6) may grow rapidly with N and must be bounded for the ACTUAL original
profile. The formula does not invent those bounds or replace the profile.
Any original width/normalization dependence in sigma_k and in the other
operator factors remains; there is no uniform-in-width claim.

## 3. Populate the raw kernel slots with explicit remaining inputs

The output multiplier in Q110's raw e=0 tail is B_n=N_n(k0 chi_U).
The raw Leibniz rule and (3)--(6) imply

    B_n <= 2^n [3 H2 n! (3/u)^n]
                 [(n+1)! (512(1+Xb)/u)^(n+1)
                                     sum_(k=0)^(n+1) sigma_k],
    A_J=N_J(f0) <=3 H2 J! (3/u)^J.                      (7)

The cutoff input for this PARTICULAR output factor is through n+1,
not J. That does not reduce Q110's much larger source/shape/profile graph:
r0 itself contains waves, pressure-generated means, cut curls and cutoffs
whose derivatives still run through the high-order requests. In particular
profile order 100046 is not erased by (7).

For the existing J=100005 source request, the previously unknown power
factor A_J now has the closed expression (7). Replace that profile slot
by this expression; do not multiply it by an additional independent A_J.
Q111 already retains a J! in its forward chart and Q110 retains a J! for
through-to-raw source conversion. Including (7) yields THREE conservative
source-side factorial instructions with argument J. This is not a sharp
factorial-growth claim and no such giant factorial is evaluated. Output
n! and (n+1)! factors in B_n and Q111's output chart are also retained.
The inherited 2^n raw-product factor is recorded explicitly.

The raw U integration, irrational divisor, fixed-label derivatives,
whole-line smooth extension, source-pole unweighting, original geometry,
length, output edge floor and temporal inverse floor remain unchanged.
An arbitrary KA0 constant is not substituted for B_n A_J. No A1 pressure
or velocity term is dropped: the whole current source is still Q110's
six-mean/total-covariance r0; K0 is removed only under A0, as previously
proved, and remains in pressure and all other components.

## 4. Bind to the current hull; keep actual certification gates

`KernelProfiles` takes Q111.ChartBridge, checks identical source state,
uses its actual d/R hull and the current profile's original Xb. Its
`plan` supplies closed f0,k0 and cutoff recipes, a shape-sum requirement
through n+1, and the three source factorials. Only finite scalar bases
and derivative-order arguments are computed. `materialize` is a capped
moderate diagnostic and requires a supplied nonnegative shape-jet sum
for the cutoff; it does not certify that supplied sum against phi.

`recipe` delegates to Q111/Q110's SAME high-jet archive, bundle-length,
exact pole, frozen-width and original geometry guards before attaching
these profile majorants. It does not enlarge the norm reserve by changing
a maximum attribute, mix state objects, or certify raw source jets from
smoothness alone. The original coordinate inequalities (5), actual
shape-sum bound, actual source-through-J bound, full common-pole/archive
proof and all whole-hull/mixed-trace contracts still need their legitimate
native inputs. Passing these metadata checks is not that proof.

For fixed positive a,b,Xb, the power constants are uniform in d on its
original range and independent of S/epsilon. The cutoff composition
constant has the same fixed-hull character conditional on (5) and the
original shape jets. Uniformity over physical scales additionally needs
uniform NORMALIZED hull/shape controls and all original physical transfer
losses; this checkpoint supplies no actual-pulse theorem.

## 5. Verification and next obligation

Exact rational checks compare (3) against Q44's raw falling-factorial
bounds across fractional exponents and hulls on both sides of R=1.
Independent differentiated fractional-power fields check the raw jets.
Moving-cutoff checks compare (6) against the complete Q46 recurrence and
detect removal of the slow-motion term. Exact inequalities check the
chosen majorant radius; mixed derivatives check the raw product rule.
A typed current graph verifies the original hull, e=0 channel, deferred
large-order request, shape order and all three source factorials while
retaining old archive guards. These are finite diagnostics, not a
continuum proof of the actual shape/source inputs. Frozen replay checks
unchanged Q111 history and all source bindings.

Next establish the actual fixed-shape jet controls and the current
source-through-J/high-pole bounds, now with (3), (6) and (7) made explicit,
or prove a sharper actual A0 cancellation. No actual A0 gain or full cone
join follows yet. Native coefficient/reserved-mean supports, uniform
actual-pulse size/derivative/support across all physical scales, complete
viscous/cutoff/interaction residual improvement, infinite realization and
blow-up remain unproved.
