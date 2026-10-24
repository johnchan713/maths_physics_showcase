# Symbolic weighted compact maps on the original shell

Parent: `6d77c8286d38d94b380bd264fdc06158263960c6`.
On 24 October 2026 its reverse-chart CI step passed; the full run was
still in progress without failures. This checkpoint implements the
weighted compact-map constants without evaluating enormous exponential
weight floors. It gives a sufficient estimate relative to the SAME fixed
profile jets, plateau gap and whole radial hull. It does not certify those
profile/domain premises or the complete current source. No historical
evidence or physical parameter is changed.

## 1. Keep the actual shift while estimating the weight

Use exactly the shifted operators and raw five-variable norms in
`core_fast_inverse_budget`, section 2. At fixed integration shift a,
the source is F(U+a,y+Ma vr). The translation is independent of all
differentiated coefficient variables; M and the covering index are held
fixed. Thus no M is introduced by differentiating this formula. Taking
the full torus supremum removes the translation in an upper bound, not
in the operator identity. The unshifted primitive would give the WRONG
divergence identity and is never substituted for the actual operator.

Write x=log(X/Xa), L=log(Xb/Xa), zeta=exp(-w^2/x^2-4/(L-x)^2),
delta=min(1,x,L-x), with the original w=C_setup^-20. Suppose every
required RAW source jet is bounded pointwise by
epsilon^b zeta delta^-p times its coefficient bound. Derivatives of the
weight are already in those source jets. Let J bound dU/dx on the whole
positive shell at fixed slow variables. All moving-profile and chart
jets must be included in the fixed multiplier norms below.

The actual chi is zero on x<=g and one on x>=L-g, uniformly over the
slow parameters, for a fixed g in (0,L/2]. Its differentiated support
lies in [g,L-g]. Set r=min(1,L/4,g), a=w^2, c=4,
q=max(0,ceil((p-3)/2)), and

    D(v)=1/(2v) sum_(j=0)^q binom(q,j) j! (r^2/v)^j,
    E=max(1,ceil(p/2)! (2/w^2)^ceil(p/2)),
    C=D(a) exp(c/(L-r)^2)+D(c) exp(a/(L-r)^2)
                             +L E exp((a+c)/r^2),
    Z=exp((a+c)/g^2).                                  (1)

The prior two-edge proof bounds each one-sided integral on its directed
collar by C zeta delta^-p. On the interior use its weight floor and
the whole integral bound L E. The actual cutoff selects the correct
direction near each endpoint; both terms on the interior obey the
same safe constant. No zero moment is required for this size estimate.
The sum in C safely dominates the old maximum of three branches.
Changing variables dU=(dU/dx)dx contributes J. This is only an estimate
of the absolute integral, not differentiation of a changed-shift formula.

## 2. Actual complete fast-velocity constants

At order m let

    Lt_m=2^(3m) N_m(g_e)[N_m(1-chi)+N_m(chi)]N_m(f_e),
    La_m=2^(3m) L_U N_m(k_e)N_m(chi_U)N_m(f_e),

with the original f_e,g_e,k_e and U-hull length L_U. These are fixed
profile product constants, not independent current-field operator norms.
The same raw Leibniz proof as the historical unweighted maps gives

    KT^w_e,m = J C Lt_m,
    KA^w_e,m = E Z La_m.                               (2)

For A_e use the WHOLE integral of every differentiated source jet, bound
its weight by E, and recover the output weight only on the actual interior
support of chi_U and all its derivatives. There delta^-p>=1. The latter
constant includes chi through m+1, not merely m. Its total shifted moment
is retained even when the Haar source vanishes. For e=0 pressure, the
same construction applies to gr-rho P with the WHOLE current moving P;
it does not supply a bound for P from a projected target.

Let I bound c_i^-1/S in the retained common atlas, J_m=binom(m+5,5),
Gamma_m=24*2^(m+4). For raw source coefficients Et_(m+4), Ez_(m+4),
and SEPARATELY Ez_(m+5), the historical inverse therefore gives

    Ft_m = J_m Gamma_m I S [Et_(m+4)+(1+KA^w_1,m)Ez_(m+4)],
    Fr_m = J_m Gamma_(m+1) I S KT^w_1,m+1 Ez_(m+5).    (3)

The actual axial velocity is gd-A1 gd, retained by 1+KA, and the radial
field has its original extra epsilon. If the raw source has epsilon
order b, the outputs have tangential order b and radial order b+1;
the edge pole stays p and there is one explicit S factor. Normalized
input recipes need their OWN factorial conversion before (3). The code
does not identify the two different Ez orders or the two norm conventions.

## 3. Exact representation and remaining inequality

`ExpSum` stores finite positive sums of rational c*exp(e). Addition and
multiplication are exact positive expressions, not floating-point or
rational evaluations of exp. All expressions in (1)--(3) are finite
for the original positive fixed data. A diagnostic rational upper bound
using exp(e)<=3^ceil(e) deliberately refuses giant exponents; keeping
the expression symbolic succeeds even then. The old materializer's
refusal is an execution limit, not an analytical obstruction to these
fixed constants. No uniform constant as w,g tend to zero is claimed.

The remaining actual-source inequality is, for every required raw jet,

    |partial^alpha E_current^circ|
       <= epsilon^b zeta delta^-p E_alpha(S)

on the SAME whole shifted integration hull, together with original
fixed-profile jet bounds and a plateau gap g>0 on that hull. Neither
unweighted source smallness nor the tensor chart transfer proves it.
The reserved slow potential hull and its positive weight floor must
likewise be supplied before the slow response can be joined. Expressions
here do not yet extend the historical rational-polynomial row code's
coefficient ring; they are rigorous symbolic sufficient constants for
that pending substitution, not a claimed completed mean bound.

Next join those current source/slow scalar estimates and the full pressure
moment using their actual pole ledger (including product poles 2p and
4p where required). All covariance correlations, axial corrections,
reserved supports, viscosity, cutoffs and interactions remain in the
complete residual obligation. No unconditional pulse bounds, complete
residual improvement, infinite realization or blow-up is proved.

The finite exact controls compare symbolic coefficients with the old
three branches, check full shifted rather than unshifted differentiation,
and detect deleted A1 or stale radial source order. They are manufactured
diagnostics, not continuum suprema. Historical records are preserved.
