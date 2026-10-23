# Generic signed covariance, complete pressure and current-target gain

Parent: `7b1acd4e8f8d966c6eeec2031fbec1cd988ff394`.
The current prefix has full regular tangential mean rows at
C=3/2-kappa, nonzero regular source at D=C-kappa, and conditional
postslow scalar remainders at C+1. Here the NEXT generic signed wave
is joined to the same original canonical roots, the complete six-entry
covariance, reconstructed pressure and current compact targets.

The conditional projected-target order 2-2kappa was already derived
qualitatively in `current_target_cycle_bridge`. The new result is its
finite positive-polynomial coefficient/pressure recipe, the localized
original-root cross identity including both cutoffs and partition
derivatives, its connection to the recent explicit current-row/moment
recipes, and the newly expanded pressure proof domain. We do not repeat
the historical result as a new unconditional theorem.

Original weighted target/root/common-chart/high-jet/support inputs remain
conditional. The executable examples use manufactured lower norm data
and compact C5 polynomial fields, not the actual enormous flat profiles.
Historical evidence and the fields, original roots, width C_setup^-20,
carriers and accepted common band are unchanged. This checkpoint stops
AFTER the signed wave and its recomputed pressure; it does not implement
or certify the subsequent fresh fast/slow means.

## 1. Linearize with the original roots on the entire localized cover

Let H contain the SAME two canonical column vectors, with their complete
temporal/transverse cutoff integrals, angular half and physical Haar
Jacobian. For one box write a for the original positive amplitudes and
let Xi be its current normalized two-component target. Use exactly

    d=H^-1(Xi/scale), delta a_j=d_j/(2a_j).             (1)

The scale here is the physical covariance factor for that box; equivalently
it can be absorbed into Xi using the original physical target normalization.
Do not use a+delta a in the denominator. The two selected principal
cross-covariance entries are scale*H*(2a*delta a)=Xi. Multiplying by the
original slow squared partition and summing boxes gives the full physical
current target because sum eta_box^2=1 and each normalized box target
has the same physical value. Both signs belong inside a box's H, not in
an extra doubled partition. All derivatives of this identity include
the actual cutoff columns, inverse, reciprocals, normalization and
partition. The signed square is not cancelled by (1).

This is a linear signed map about the fixed positive original roots;
the updated amplitudes need not stay positive. The other four entries
of its principal cross tensor generally do NOT vanish. They cannot be
dropped from pressure or scalar updates. Exact controls independently
integrate two polynomial cutoffs, differentiate a nonconstant squared
partition and detect an updated-root denominator or missing factor two.

## 2. The complete covariance error, including complex quadratures

Let W0 be the ORIGINAL principal wave, U the current COMPLETE wave,
E=U-W0, and V=L+R the new signed complete curl. L is its principal part;
R contains every curl-of-amplitude/cutoff/connection contribution. With
Cov the angular average and Bil(a,b)=average(a tensor b+b tensor a),

    Delta W=Cov(U+V)-Cov(U)=Lambda+Ecov,
    Lambda=Bil(W0,L),
    Ecov=Bil(E,L)+Bil(U,R)+Cov(V).                     (2)

This is an exact identity for all six entries, with any current complex
quadratures and harmonics. No real-primary cosine/sine shortcut is used.
For through-order harmonic-l1 coefficient norms on one common weighted
chart, sufficient positive polynomials are

    Lambda_n <= 2 T_n L_n,
    Ecov_n <= 2 E_n L_n+2 U_n R_n+(L_n+R_n)^2,
    Delta_n <= 2 T_n L_n+Ecov_n.                      (3)

T bounds the original principal wave, E its actual error summands and U
the current complete wave. A bound for U alone does not bound W0 or E:
possible cancellation prevents that inference. One valid route for E
adds the ORIGINAL primary curl remainder, complete first signed update
and actual generated sourced wave; their retained orders are at least B.
Their original weighted/common-chart coefficient bounds must be supplied
by their same-family lower recipes, not by unrelated norm subtraction.

With B=1-kappa, C=B+1/2, D=C-kappa, the relevant wave orders are

    W0,U: 1/2; E: B; L: C-1/2=B; R: C-kappa=D.

Hence Lambda and Delta have order C while the three errors in (2) have
orders 2B, 1/2+D and 2B, all equal to C0=2-2kappa. The bound is for every
fixed finite jet, with the corresponding finite polynomial and pole costs.
It is not a uniform bound in stage, jet, physical radius or width.

`RootSignedEnvelope` calls the ORIGINAL `RootBudget.signed` recurrence
on `TargetBudget.V`, including reciprocal-root and target-chart derivatives.
`PairedTargetBudget` keeps the e=2 angular and e=1 axial constructions
separate, sums their coefficient bounds and takes the larger edge
allowance; the axial component uses its whole pressure moment.
It then applies the SAME first-order complete curl matrices as the primary
budget. Its signed curl at n reads the principal at n+1. A retained
weighted common-chart/label conversion is applied to the through-order
sum. This conversion is not inferred from unweighted coefficient bounds.

## 3. Reconstruct pressure before projecting the new target

Write D_ab=bar(Delta W_ab), with bar the auxiliary Haar average. Mean
velocities are held fixed at this substep. On the same whole-line compact
extension and same moving unit-mass rho, the exact changes are

    dg=-D1 Drr-epsilon partial_Z Drz+Dtt/R,
    dP=int dg, dp=I0(dg-rho*dP),
    dJtheta=int R^2 Dtheta,z,
    dJz=int R Dzz-(1/2)int R^2 dg.                    (4)

All six entries are retained. The actual full pressure still contains
the unchanged mean/base components, Pi1/Pi2 and previous wave pressure.
Formula (4) is a CHANGE, not a replacement of that pressure by a leading
term. Its moving pressure moment satisfies

    epsilon partial_Z int R(Dzz+dp)
       =epsilon partial_Z(dJz+c_rho*dP),
    c_rho=(int R^2 rho)/2.                            (5)

The derivative of c_rho is kept. Haar pressure change is order C; full
pointwise pressure change can be only D because its graph source loses
kappa. Neither is promoted to C0. No compact pressure velocity is erased.

For Qe f=f-sigma_e int R^e f, int R^e sigma_e=1, and
Ie f=R^-e int_left^R r^e f, the current target is Xe=-Ie Qe bar(Ee).
The exact signed principal matching in section 1 and (2) give

    Xtheta,new=-Ecov_rtheta-I2 Q2 epsilon partial_Z Dtheta,z,
    Xz,new=-Ecov_rz-I1 Q1 epsilon partial_Z(Dzz+dp).    (6)

Projection is applied AFTER differentiation. In general

    Qe epsilon partial_Z H
      =epsilon partial_Z(Qe H)
            +epsilon (partial_Z sigma_e) int R^e H.   (7)

Dropping this term or the derivative of c_rho changes the actual map.
Both axial responses in (6) have order C+1>C0. Thus the new CURRENT
projected targets have order C0, conditionally on their original weighted
operator/input inequalities. The current moment identities reconstruct
both new Haar tangential rows at C0 as well, since their moment channel
is order C+1. For row jets request the target recipe one order higher.

In contrast, the scalar CHANGES in (4) are generally order C, not C0;
they can dominate the old C+1 defects. The post-signed radial Haar
remainder -rho*Pnew is therefore only order C. Exact controls show a
nonzero dP and dJz despite the selected principal cancellation. Fresh
post-signed fast and then CURRENT slow corrections are indispensable.
No complete mean/residual improvement follows from (6) alone.

The function `pressure_target_budget` bounds (3)--(6) explicitly. At
target output n it uses D through n+2 to bound pressure through n+1.
For example sufficient source coefficients are
G_n=2(n+1)D_(n+1)+2 r D_n and
G_(n+1)=2(n+2)D_(n+2)+2 r D_(n+1). With M for all moment calls,
Pm_(n+1)=Kp[G_(n+1)+rho M G_(n+1)]. Its axial target coefficient is
Iop*(1+sigma M)*(n+1)*(D_(n+1)+Pm_(n+1)), multiplied by the EXTRA
epsilon. Iop must be the original oriented weighted Ie estimate, with
its radial and differentiated-kernel factors, not an arbitrary unweighted
primitive norm. The theta formula omits only the axial pressure addend.

## 4. The next pressure proof domain must be enlarged explicitly

The recent prefix's pressure/moment proof domains were 6p and 12p.
The next signed target has its own edge allowance d_N from the original
`TargetBudget.edge(N)`. At target output n, full pressure through n+1
uses covariance/curl coefficients through n+2 and signed target input
through n+3. A sufficient common wave pole for these finite calls is

    q=max(3p, 6(n+2)+9, 6(n+2)+d_(n+3)+6).             (8)

The three entries cover the old complete wave, original curl and new
signed curl. Products in (2) then have source pole 2q. The previously
declared 6p/12p operator proof is NOT automatically a proof at 2q.

For illustration only, p=1 with a common incoming row pole p_j=12
gives d_(n+3)=12; at n=0, q=30 and the new pressure source pole is 60,
where the old largest declared domain was 12. These are sufficient
finite bookkeeping bounds, not minimal necessary poles or evaluated
actual high-jet inputs. They do not replace the parent's enormous raw
A0/A1 source derivative requirements.

The gap can be closed at the level of CONDITIONAL operator estimates:
invoke the existing exact kernel `shift_recipe(q/6,L,width,gap)`. Its
6(q/6),12(q/6) domains are q,2q. Its original two-edge/raw-profile
constants dominate the required pressure/moment constants at those
domains. Reevaluate every nonnegative lower/input/response polynomial
at S+K_new. For the new target primitive also retain its separate
`primitive_recipe(2q,L,width)` and the original TargetBudget.F/K
derivative/radial factors. Its three branches have the finite upper
K_Ie=sum(branch coefficient * 3^ceil(branch exponent)). Use the common
estimate argument S+K_new+K_Ie. Positivity gives domination on S>=1;
there is no new S>=K_new gate, actual pulse change or band choice.
The original width, cutoff plateau gap and whole potential hull remain.

`signed_pressure_shift` retains that exact finite power recipe. Only
moderate illustrative constants are materialized in diagnostics. For
actual giant pole/width requests the original recipe and its nested
inputs remain deferred arithmetic instructions, not evaluated bounds.
The pressure TARGET primitive thus retains its own oriented Ie estimate
on the resulting source/pressure pole, with all finite requested jets;
the ordinary algebra of (6) supplies no replacement for that estimate.

## 5. What is and is not joined to the current prefix

`CurrentSignedRecipe.target_inputs` references the SAME
`FinalMeanRowClosure.tangential_recipe` and
`PreservedHaarScalarJoin` whole moment recipes, at the stage after the
latest fast/slow updates. It rejects different mean-prefix objects.
The current complete rows are C; the parent scalar/pressure moments
with their axial epsilon are C+2 and are stronger than needed for
constructing the C-order target. Original bump, oriented primitive,
reciprocal weight and chart estimates still enter `TargetBudget`.
No new fine current-row callback or different positive root is supplied.

The actual weighted target, original W0/error summands, common-chart
and label norms, compact endpoint/support identities and background
premises remain unverified lower inequalities. The recipe refuses to
certify them. Tests compare finite exact algebra and manufactured
polynomial majorants, including conjugate +/-1,+/-2 harmonics; they
are not continuum enclosures.

Next prove the actual post-signed fast/slow difference bounds, refresh
all three scalar rows and the pointwise radial repair, and join the new
raw compact-tail proof domains before declaring the canonical seed at
C0 complete. The complete regular source remains at D until a full
cycle proves more; actual input validity, sufficient physical regular
orders, controlled infinite diagonal/force regularity and blow-up remain
open. In particular no absolute residual decay is claimed here.
