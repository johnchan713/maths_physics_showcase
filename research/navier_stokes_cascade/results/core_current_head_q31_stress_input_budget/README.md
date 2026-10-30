# Q31: explicit finite-head bounds for the weighted background stress inputs

Parent: `e026133f29dc758d04a269166e0ff8e46a30dd4e` (Q30).
This supplies finite constants for the leading normalized target and higher
background stress used below the first signed-source budget. The historical
weighted-stress foundation proves existence of these bounds; here the
exceptional first coefficient, compact higher coefficients, weight derivatives
and chart conversion are separated into explicit recipes. This is a norm
estimate for the SAME rule-defined stress, not a new Borel selection or a
proof that its associated velocity has zero residual.

All conclusions are relative to the retained exact leading construction,
globally finalized coefficient equations and original combined Borel rule.
The root parameter remains h=exp(-8*(exp(64)+10))>0. Rational h values in the
tests are manufactured controls. No original parameter, fixed support,
pressure gauge, accepted band, historical file or width C^-20 is changed.

## 1. Identify the two distinct quotient families before estimating

Use E=q partial_q at fixed X,eta and D=(E,partial_X,partial_eta).
Write a=t1^2, ya=log(X/Xa), yb=log(Xb/X),

    zeta=exp(-a/ya^2-4/yb^2), delta=min(1,ya,yb),
    Tsum=T0+sum_(n>=1) chi(a_n q) q^(2nh) Tn.

The exact inner/outer factors in `weighted_stress_foundation` imply raw
through-order bounds, for every fixed m,

    |D^I T0| <= A_m zeta delta^(-3m-3),
    |D^I T1| <= B_m zeta delta^(-3m-3),
    |partial_(X,eta)^I W1| <= V_m, Wn=Tn/|T0|.       (1)

Here A_m,B_m,V_m are finite continuum seminorm bounds, not sampled maxima.
W1 extends smoothly to both flat edges. T1/zeta need NOT be bounded:
at the outer edge T1/zeta contains a yb^-3 factor where eta is nonzero.
Thus it must not be inserted as an unweighted compact-quotient input.
For n>=2, both Wn and Zn=Tn/zeta are smooth with the same compact interior
support K=[Xinner,Xv]. This support uses the full forward pressure and
conservative moment identities, including shifted axial viscosity.

Let C_(n,m)>=1, n>=2, bound EVERY D jet through m of both cutoff
increments chi(a_n q)q^(2nh)Wn and chi(a_n q)q^(2nh)Zn after removing
q^(2nh). One sufficient value is the existing
2 K_m(1+2nh)^m times the larger of their two profile seminorms, padded by
one; K_m>=1 bounds the fixed Euler cutoff jets through m. These are stress
increments: no extra cutoff derivative is introduced by this direct sum.
Velocity still needs K_(m+1) and its radial potential commutator.

The existing combined rule, with ALL velocity/potential/physical rows,
requires C_(n,m) a_n^(-nh)<=2^-n for m<=n. The two stress quotients use
the SAME a_n as velocity. This follows from the finite maximum prescribed
in `weighted_background_edges` (8), `background_jet_foundation` (11), and
the rule identification in `core_leading_background_borel_join`; it is
not inferred from the velocity diagonal alone. A previously fixed sequence
not identified with that combined rule would still have precisely these
inequalities to verify. No replacement sequence is selected here.

## 2. An explicit head for every order on one entire q interval

For 0<q<=q0<=1, put N_m=max(2,m) and

    H_m=sum_(2<=n<N_m) C_(n,m)+2^(1-N_m).             (2)

Both sums over n>=2 (the Wn sum and the Zn sum) have every D jet through m
bounded by H_m q^(2h). For the finite head use q^(2nh)<=q^(2h). For the
tail, all its indices obey n>=m, so the rule gives

    sum_(n>=N_m) 2^-n q^(nh)
      = q^(2h) 2^-N_m (q^h)^(N_m-2)/(1-q^h/2)
      <= q^(2h) 2^(1-N_m).                            (3)

This is valid on cutoff transition collars as well as where chi is one.
It proves termwise finite-order differentiation by uniform domination;
local finiteness for q>0 also follows from scale separation. There is no
order-dependent shrinkage of q0. The exceptional n=1 contribution remains
finite-head data at EVERY derivative order. Neither a uniform-in-n
coefficient bound nor convergence of the uncut formal series is needed.

In particular the raw relative perturbation R=(Tsum-T0)/|T0| obeys

    |D^I R|<=q^(2h) R_m,
    R_m=max_(0<=j<=m){K_j(1+2h)^j V_j+H_j}.           (4)

This gives a concrete sufficient value gate R_0 qstar^(2h)<=1/2 for
|Tsum| between |T0|/2 and 3|T0|/2. The existing qstar and its stricter
Q2 velocity/shear/cone gates must already satisfy their combined budgets;
(4) does not authorize silently shrinking an identified historical field.
Taking more derivatives changes constants, not this value gate.

## 3. Weight and reciprocal-weight derivatives have the same majorant

Fix 0<Xlo<=X on the shell. With unsigned first-kind Stirling numbers c(j,k),
define for j>=1

    L_j=(a+4) Xlo^(-j) sum_(k=1)^j c(j,k)(k+1)!,
    E_0=1,
    E_j=sum_(i=1)^j binom(j-1,i-1) L_i E_(j-i).        (5)

The identity partial_X^j=X^-j sum_k s(j,k)partial_(log X)^k
and |partial_y^k y^-2|=(k+1)! y^-k-2 prove
|partial_X^j log zeta|<=L_j delta^-3j: k+2<=3j and delta<=1.
The same bound holds for -log zeta. The differentiated exponential
recurrence therefore proves BOTH

    |D^I zeta|<=zeta E_j delta^-3j,
    |D^I zeta^-1|<=zeta^-1 E_j delta^-3j, |I|=j.      (6)

Non-X derivatives vanish. Formula (6) does not bound zeta^-1 itself by a
polynomial or continue it smoothly to its own edge. Nor is delta ever
differentiated: it is only a pointwise upper-bound notation. Both exact
flat constants a and 4 remain in (5).

Leibniz, (1)--(6) now give explicit constants

    T_m=max_(j<=m) sum_(i=0)^j binom(j,i) E_i A_(j-i),
    S_m=max_(j<=m){K_j(1+2h)^j B_j
                         +sum_(i=0)^j binom(j,i)E_i H_(j-i)},

    |D^I(T0/zeta)| <= T_m delta^(-3m-3),
    |D^I(Tsum-T0)| <= q^(2h) zeta S_m delta^(-3m-3),
                                             |I|<=m. (7)

For the first term in S_m, only E differentiates the cutoff and q power;
profile derivatives act on T1 and are already in B_j. For the compact
tail, write it as zeta times the Zn sum before differentiating. Its true
pole is at most 3j; padding to 3m+3 permits the common first-order pole.
For T_m the poles add to exactly 3j+3, not 6j+3. Applying edge unweighting
first would lose the full zeta weight needed by the mean maps.

## 4. Transfer the bounds without hiding coordinate costs

In the enlarged band chart s=q/Q in [1/4,4], ell=log(s/s_ref),
X=R^2/(2s), eta=Zs^(-(1/2-h)), let Comp_m be equation (7)'s finite Taylor
composition constant in `core_background_finite_head_budget`. Its
log/X/eta increments have zero constant term. Composition of a raw
through-m bound in D with those finite Taylor jets costs

    Jchart_m=m! Comp_m.                               (8)

No external velocity factor s^-A or 1/R belongs in this scalar stress
composition. Each profile derivative of order k<=m in (7) is bounded at
the base point with delta^-3k-3<=delta^-3m-3. Thus every raw slow derivative
through m of T0/zeta is bounded by Jchart_m T_m delta^-3m-3, and every
such derivative of Tsum-T0 by

    4 epsilon^2 Jchart_m S_m zeta delta^-3m-3,          (9)

since q^(2h)=epsilon^2 s^(2h)<=4epsilon^2. This conversion keeps the
entire weighted stress expression inside the differentiation; it does not
freeze zeta or q. The same base-point finite-jet argument proves (8)--(9)
at one-sided eta endpoints. It uses finite real jets, not an analytic
extension of the heat profile or an infinite Taylor series.

These are R,Z,T-chart scalar bounds. In the root recipe's five-variable
chart, AFTER its retained R/U and atlas conversions, a raw exact-order-j
bound C_j supplies factorial-normalized target input 5^j C_j/j!.
Those coordinate conversions must be applied, not silently assigned cost
one. Likewise (9) supplies the higher-stress input below the initial Haar
row budget after its retained physical normalization and atlas conversion.
It does not itself identify the dimensionless Tsum with physical Sigma,
or omit the physical stress prefactor. The original leading target remains
T0, not the fully summed corrected stress.

## 5. Scope, checks and next obligation

This removes an unspecified weighted-stress O_m constant in favor of
(1)--(9) for the combined rule-defined background. Finite profile seminorms
remain analytical inputs supplied by the exact regional/factored providers;
none is numerically evaluated at the enormous root parameter. Existing
velocity norms and their additional cutoff derivative remain separate.
Full Pi1/Pi2 memory, axial velocity, moving pressure moments and every
later coefficient remain in their forward reconstruction. A direct stress
bound proves no equality between summed stress and the complete residual
of a potential-cutoff velocity.

Exact tests independently differentiate the two-sided exponential weight,
verify the Stirling/Bell and sharp geometric-tail formulas, compare a finite
mixed-jet stress fixture against (7), check the scalar chart transfer,
and reject using the n=1 quotient as compact data or changing derivative
conventions. They are diagnostics, not continuum norm certificates.
The evidence replays Q30 and its frozen transitive sources without editing
historical records. At parent inspection on 30 October 2026, focused CI
41076091890 passed; broad CI 41076091871 was still running.

Next join these explicit scalar bounds through the actual physical/atlas
normalizations and canonical entrance/column/curl recipes, checking the
remaining original frame and source input inequalities on Dv. Full actual
source closure, all-scale pulse bounds, complete residual improvement and
blow-up remain unproved.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json` in this directory.
