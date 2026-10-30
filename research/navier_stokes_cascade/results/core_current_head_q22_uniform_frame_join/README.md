# Q22: uniform frame budget, and the missing compactness quantifier

Parent: `eb6f2073673911fe69d26c803ff3e7dd8f354cf5` (Q21).
Q21 transfers a source bound along a lift; it cannot supply the magnitude
of that source. Tracing the first source back through the canonical pulse
finds an earlier quantifier that must be made explicit. Q3's compactness
argument proves finite constants on each closed band, not automatically
one constant over all bands. This checkpoint corrects that inference and
provides a quantitative repair using the existing finite-head background
bounds. It preserves the original field, exact root parameter, selected
bands, fixed width and all historical evidence.

## 1. What the whole-path inverse actually needs

For the inherited actual-normal frame, the normalized tangent equation is

    y'=[diag(0,-2 lambda)+E-(d-dref)I-(m^2-1)d I] y,
    y=z/P, P'/P=lambda-dref.

The complete normal/frame connection terms remain inside E. With
Kframe bounding both endpoint frame factors, the established energy proof
bounds the physical propagator by

    Kframe exp(CE+Cd) P(v)/P(w),
    CE=L sup_[0,L] ||E||, Cd=L sup_[0,L] |d-dref|.       (1)

For a band-polynomial estimate it is sufficient to prove, uniformly over
all labels and accepted bands,

    log Kframe + CE + Cd <= A + B log S.                (2)

Here A,B are fixed finite constants. This is a sufficient condition for
this majorant, not a necessary condition for the true propagator: sign
and matrix cancellations can make (1) very pessimistic. The original
frame provider uses the stronger condition that Kframe, CE, Cd are
independent of S. Q3's statements that the paths are compact and their
smooth suprema are finite do not, alone, prove either condition.
A finite polynomial bound for E itself also does not supply E=O(1/L).
Picard existence on each finite interval remains valid.

For an optional sharper energy budget, set
 a(v)=lambda_max(Herm E(v))-(d(v)-dref(v)).
Replacing CE+Cd by integral_0^L max(a(v),0) is sufficient by the same
Hermitian pairing. This checkpoint does not claim that this sharper
quantity has been bounded for the actual field.

## 2. Exact obstruction to the compactness-only implication

Take L=S=ell^2, lambda=2, d=dref=3/2+v/L, and

    P(v)=exp(-(v-L/2)^2/(2L)), E=diag(1,-1).

Then the normalized equation is y'=diag(1,-5)y. In particular

    z(v)=P(v) exp(v) e1,
    z_forced(v)=P(v)(exp(v)-1)e1, z_forced(0)=0           (3)

solve the homogeneous equation and the equation with added source P e1,
respectively. The added source has normalized size one. In the repository's
convention (source on the negative right side), it is f=-P e1. Every
coefficient has uniform fixed-order ordinary v-jet bounds; d>=3/2, the
frame is constant, and P has its maximum one at the midpoint. Nevertheless
CE=L and the Gaussian-normalized propagator is exp(v-w), not polynomial
in S. This is a smooth two-coordinate ODE control, NOT a proposed
Navier--Stokes background or a counterexample to the retained frame lemma.
That lemma's uniform CE hypothesis fails here.

The example does not depend on this particular Gaussian: for any prescribed
positive P with P'/P=lambda-dref, adding this E to the reference generator
still makes its first normalized coordinate grow as exp(v). Thus using
the actual reference envelope would not repair the logical implication.

No fixed epsilon reserve can absorb the model loss. Since epsilon=2^(-h ell),
for any fixed h>0, reserve delta>=0 and polynomial exponent M>=0,

 log[exp(ell^2) epsilon^delta / ell^(2M)]
   =ell^2-h delta log(2) ell-2M log ell
   >=ell^2-(h delta+2M)ell, ell>=1.                     (4)

The last expression tends to infinity; for ell>=2(h delta+2M) it is
at least ell^2/2. The forced factor exp(ell^2)-1 has the same failure
(up to a fixed factor for ell>=1). Tests check exact identities and
rational lower bounds, rather than using floating-point overflow as proof.

## 3. Repair from the original finite-head background theorem

Use the SAME fixed profile, complete coefficient tuple, Euler cutoff
schedule and enlarged slow chart as `core_background_finite_head_budget`.
Its equations (4)--(8) give constants independent of band:

    ||D^<=2(Flead,Glead)|| <= Jlead_2,
    ||D^<=1(F-Flead,G-Glead)|| <= epsilon^2 Jerr_1,
    ||D^<=1(F,G)|| <= Jfull_1.                          (5)

These component-sum majorants include all required components separately.
They follow from finite leading seminorms, the coefficient-one head and
the geometric Borel tail; they are not suprema over increasing path lengths.
The same rule-defined schedule is identified by
`core_leading_background_borel_join` and the Q1 same-object chain. No
replacement h=1/200 is used in the analytical statement. Rational fixture
parameters in tests are manufactured only.

Let the maximum slow-coordinate displacement from the frozen representative
be Dbox S^-3. Assume its connecting segments lie in the enlarged chart,
as required by the inherited phase proof. The mean value theorem in the
THREE slow variables, applied also to F_R and G_R, gives

    |F-Flead(rep)|, |F_R-Flead_R(rep)|,
    |G_R-Glead_R(rep)|
      <=3 Dbox Jlead_2 S^-3 + Jerr_1 epsilon^2
      <= Db_req (S^-3+epsilon^2),
    Db_req=3 Dbox Jlead_2+Jerr_1,
    Mz_req=Jfull_1.                                    (6)

The factor three cannot be dropped for a maximum-coordinate box.
The Mz bound controls both F_Z and G_Z. Frozen reference values are taken
from the leading profile, not from a different full background. All axial
terms remain. Background radial curl commutators and Pi1/Pi2 memory stay
in the complete coefficient construction supplying (5).

Equations (5)--(6) explicitly fill the Db,Mz slots of
`frame_column_foundation` (1)--(9). Its rational formulas then give
CE,Cd,Kframe independent of S, proving (2) with B=0 on their accepted
value domain. Therefore Kframe*3^ceil(CE+Cd) is a finite constant recipe;
there is no need to allocate this enormous integer or evaluate a selected
pulse numerically. Existing tangent and canonical entrance-jet recurrences
can use this bound, with their remaining source/root/support premises.

## 4. Do not silently replace constants inside an old selector

The original selected constants and gates must remain attached to their
original field. A sufficient way to reuse them without changing any
accepted band is to verify

    Db_req <= Db_original, Mz_req <= Mz_original,        (7)

on the same chart and coefficient schedule, with the already retained
floors and original length/rounding gates. `join` reports these two tests;
`require_reuse` refuses to issue a uniform frame recipe when either fails.
It NEVER modifies Db_original, Mz_original, a carrier, length or first band.
The constants in (6) are sufficient upper bounds, so a failed comparison
in (7) does not establish that the actual field violates the original
bound. A sharper direct estimate, or proof of (2) on the SAME accepted
domain, is an alternative. Recomputing larger CE and silently reusing the
old normal/Riccati selector would be invalid.

The same-profile selected suprema and inequalities (7), or an alternative
uniform proof, have NOT been evaluated/certified here. This is the precise
remaining interface for this repair route. The old frame theorem remains
correct under its original uniform primitive hypotheses; the correction
is to inferring those hypotheses from per-band compactness alone.

## 5. Scope and reproduction

Established: the compactness-only implication fails; a band-independent
Db/Mz recipe follows analytically from the retained finite-head theorem;
its guarded substitution gives a uniform frame constant if (7) and the
original gates hold. The actual first sourced block's complete extended
bound remains conditional. Q19's harmonic half weight and full-weight
mean/covariance distinction, Q20's Gaussian quotient terms and Q21's
liftwise derivative transfer are unchanged. Forced pressure still uses
n.f and the moving-normal term; full mean pressure, c_rho P, all six
covariances, A0/A1 and endpoint cutoff tails remain separate obligations.

Next verify the same-profile constant/domain interface (7), or establish
(2) directly, before promoting the whole-path inverse and canonical pulse
bounds used by the first source. Complete representative/endpoint bounds,
residual improvement, force regularity and blow-up remain unproved.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json` here. SymPy is used only for exact
manufactured ODE/mean-value controls. Historical tests and provenance are
replayed; none of these finite checks certify the actual continuum inputs.
