# Q125: uniform angular reciprocal cost and the same first cut seed

Parent: `2e8e6230f0401fe7964ca13a6cf520c7e46415ce` (Q124).
On 3 November 2026 the freshly fetched branch agreed with that parent.
Q124 research CI 41641232266 and field CI 41641232330 were in progress;
Q122 field CI 41622757388 had passed. Current notes, frozen evidence,
axis_first_coefficient, first_coefficient_c3_boundary and the coefficient
support foundation were inspected before this continuation.
Historical frozen evidence, roots, widths, supports, pressure datum and
scales remain unchanged. The core_growth_activation conclusion remains
local leading core/first collar with width C_setup^-20.

Q124's all-order reciprocal estimate is valid but much too large. The
already established factorial angular product algebra gives a bound
independent of the requested derivative order. We insert it into the
current full-pressure discrepancy estimate and derive the all-angular-order
FIRST coefficient cut-seed bound from its existing Cauchy cover. This
reduces that first coefficient's incoming repair inputs to actual leading
whole-support norms and the SAME core/collar solution bounds. It does not
supply those selected bounds numerically or extend their strip to all n.

## 1. Use the full angular product algebra

For |f|_m=sum_(j<=m) sup|partial_eta^j f|/j! on [-1,1], the product rule
gives |fg|_m<=|f|_m|g|_m. For every fixed m,

    |2h eta^2|_m <= 2h(1+2+1)=8h<1.

The Neumann series for 1/(1-2h eta^2) converges in this finite-order Banach
algebra. Its sum is the same pointwise reciprocal, so

    |L^-1|_m <= R:=(1-8h)^-1 <=25/23,
    L=1-2h eta^2, 0<h<1/100.                         (1)

This fact was already stated in first_coefficient_c3_boundary (7); it is
not a new regularity hypothesis or a new theorem about arbitrary profiles.
In particular Q124's R_m=4^(m+1) may be replaced by R at EACH occurrence,
including R_(m+1) and R_(m+2). No order-dependent enlargement of a field
or scale is involved. Q124's frozen record remains a correct coarse bound.

Retain all of Q124's full-support hypotheses, input conventions and
identities. Its regular v_i estimate becomes

    v_i^*=4R(m+3+D+lambda_i)u_i,
    u_i>=sup_[0,Xv]|U_i|_(m+2), D=1/2-h.             (2)

For the shifted predecessor at pressure n>=2,

    p=4R(m+4+D+lambda_(n-2))u_prev,
    u_prev>=sup_[0,Xv]|U_(n-2)|_(m+3).

The retained two-angular-operator estimate is

    |T_(b-D)T_b f|_m
      <=16R^2(m+1+|b-D|)(m+2+|b|)|f|_(m+2).         (3)

Substituting (1)--(3) in Q124's Q0,Q1 and all five mu_i gives the new
bounds implemented here. The axis contribution and both shifted axial
terms remain; the latter first occur at pressure order n=2. Ordered
products, full pressure memory through Xv, support endpoints and lower
squares are unchanged. These are signed full moments, not absolute source
integrals or pointwise pressure bounds.

With native input BUDGET VALUES held fixed, the extra operator cost is
polynomial of degree at most three in m: (2) is linear, the mixed/velocity
square terms are at most quadratic, and the shifted term (3) times p is
at most cubic. Actual input norms u_i(m+2), u_prev(m+3), f_i(m) can grow
with m. This does NOT assert polynomial growth, analyticity, or a uniform
bound for the actual coefficient profiles.

The exact-number API can now evaluate signed-moment bounds at m=100005
using rational arithmetic of polynomial size in m, without forming R_m
or m!. Raw discrepancy conversion still costs m! as in Q124; it remains
capped at moderate m, or uses the unchanged binary factorial budget.
The binary reciprocal cost is ceil(log2 R), independent of m. No raw
factorial is discarded and no all-order native input is supplied by this
computational observation.

## 2. The same first-coefficient angular cut seed at every fixed order

Keep the selected delta=rho_axis/(800S), and the existing DISTINCT bounds

    Mcore=(z_core+6 z_core^2/delta) exp(6 z_core^2/delta),
    Mcollar=(z_col+12 z_col^2/delta) exp(12 z_col^2/delta).

These are the original complex-core and actually activated collar bounds;
no numerical value or alternative radial domain is substituted. Their
construction and the same-axis Volterra uniqueness are in
axis_first_coefficient and global_first_coefficient_c3. The core and collar
solutions cover the SAME whole inner interval through Xcut; each has the
retained parameter-holomorphic neighborhood with conservative Cauchy disks
of radius delta/2 about the real eta interval. Holomorphy is in eta. No
radial analyticity across flat joins is asserted.

For G=F1 or U1, Cauchy and the overlapping cover give

    sup_[0,Xcut] |partial_eta^j G|/j!
       <=(2/delta)^j(Mcore+Mcollar).

The coefficient cutoff is eta-independent and 0<=kappa<=1. Therefore the
SAME tentatively cut first coefficient satisfies, at every fixed m,

    J_(1,m) <=(Mcore+Mcollar) sum_(j=0)^m (2/delta)^j. (4)

At m=3 this is EXACTLY the inherited c3(Mcore+Mcollar), not another first
coefficient construction. For m>3 the same disks apply; no new Picard solve,
strip shrinkage, frequency, root or cutoff is chosen with m. Only angular
value norms are concluded here. High radial jets still need the separate
axis/collar differentiated equations and cutoff product estimates.
Nor does (4) imply one common complex strip for all coefficient indices n;
first_seed is specifically bound to n=1 and its existing cover.

If Mcore<=2^b_core and Mcollar<=2^b_col, a deferred binary version is

    b_seed=1+max(b_core,b_col)+ceil(log2(m+1))
               +m ceil(log2 max(1,2/delta)).         (5)

It follows by bounding each of m+1 terms in (4) by the largest padded
power. No giant exponential is evaluated by this interface. Supplied
b_core,b_col must actually bound those original solution constants; the
identity label itself does not prove the cover or these inequalities.

## 3. Close the first incoming-discrepancy interface, with its remaining input

first_discrepancies inserts Jf=Ju=J_(1,m) into the improved general moment
formula at n=1. It requires

    sup_[0,Xv]|U0|_(m+2),
    sup_[0,Xvel]|F0|_m,

on the original full domains, and the same Mcore,Mcollar,delta. There is no
shifted predecessor at n=1. The full leading radial source and axis term
are still present in its two signed pressure totals. All five discrepancy
functions, including later pressure memory, remain in the result.

For a Q121 tuple request r, Q122 needs discrepancy order m=r+1. Thus the
first repair still needs leading U0 through r+3. A local core/collar bound,
a low-order compact bound or a finite supplied grid does not replace that
whole-support high-angular input. At n>=2 the separate shifted predecessor
requires the stated m+3 input as before. The exact coefficients can then
enter the existing repair/tuple/retained-scale chain ONLY after actual
native norms, identity bindings and original scale inequalities are proved.
No accepted selector value is recertified by a smaller upper bound alone.

## 4. Validation, scope and next obligation

Checks compare the new bounds with the unchanged Q124 bounds and run
Q124's independent direct full-pressure/five-density control using the
new tighter caps. Exact reciprocal jets check (1), and a separately bounded
analytic polynomial checks Cauchy's implication (4). The m=3 first seed
recovers its historical formula. Moderate exact results are enclosed by the
binary budgets; high-order counts retain the missing leading input and
avoid factorial materialization. These manufactured controls do not evaluate
the selected giant fields and passing tests is not a continuum proof.
The frozen evidence recursively replays Q124 and its source closure.

Next prove the actual whole-support leading F0/U0 high-angular bounds
required above, and the later finalized lower-coefficient norms needed at
n>=2. Bind them to the existing inner solution constants, supports and
cutoffs, then verify the retained scales. The order-100005 pointwise source,
common-pole/trace and compact-modulation/exact-repair cone obligations remain
open. Complete-wave covariance means including axial velocity and pressure,
uniform actual physical pulses and the full viscous/cutoff/interaction
residual are unproved. No unconditional pulse estimate, residual improvement
or blow-up follows.
