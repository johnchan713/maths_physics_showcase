# Q100: complete sourced-wave row increments and the incoming-final-row gate

Parent: `d5264ba07259a2a22484b03a49e69471fe6a4069` (Q99).
Its research and field-comparison runs 41433668238 and 41433668253 were
still running when inspected on 2 November 2026. Q98's field-comparison
run 41426008465 had passed. This checkpoint derives conditional full
pointwise angular/axial row-INCREMENT coefficients from Q97/Q98. It keeps
the remaining incoming FINAL row bound explicit before any fresh inverse.
It does not mistake Q99's integrated scalar rows for pointwise sources.

## 1. Exact complete row changes at unchanged mean velocities

Put C=3/2-kappa, D=C-kappa and Qcirc=I-Haar_aux. For the complete changed
wave covariance Delta W and Q98 shifted pressure delta pm, exactly

    delta Etheta=D2 Delta Wr,theta+Dz Delta Wz,theta,
    delta Ez=D1 Delta Wrz+Dz(Delta Wzz+delta pm),
    Qcirc E_changed=Qcirc E_FINAL_OLD+Qcirc delta E.       (1)

Here De=Dr+e/R, Dr=partial_R+M d R^(d-1)vr.partial_y,
Dz=epsilon partial_Z. The 2/R and 1/R connections differ and are retained.
No Haar integration is taken before estimating (1). The full shifted
pressure, including the local -Delta Wrr term, the moving compensator
and every old pressure archive, remains as in Q98. Only the pressure
INCREMENT enters delta E; old pressure and old Pi1/Pi2 remain in E_FINAL_OLD.
The nonzero-angular sourced harmonic pressure is a separate harmonic
and does not replace this angular-zero pressure.

Both sides use the same fixed mean velocities. Therefore time, viscosity,
mean transport and assigned stresses have identical values before and
after this wave addition. They cancel in the DIFFERENCE only. They are
not removed from the incoming final row, future mean differences or total
residual. Spatial/temporal cut curls and all complex covariance products
already occur in Q97's COMPLETE tensor; no principal-only tensor is used.

## 2. Explicit weighted coefficient for the full increment

At output n, let F_n=sum_(j<=n) 5^j D_j/j! from Q97's raw covariance
majorants at its exact common pole p. Let

    r_n=sum_(j<=n) a^(-j-1),
    L_n=(n+1)(1+2g_n)F_(n+1), Z_n=(n+1)F_(n+1),
    B_n=(n+1)Kpm_(n+1),
    T_n=L_n+2r_n F_n+Z_n,
    A_n=L_n+r_n F_n+Z_n+B_n.                            (2)

g_n is the original unweighted graph-coefficient norm recipe with its
epsilon^-kappa factor separated; |vr|_1<2. The recipe is read from the
same initial budget at n+1, not assigned to a new arbitrary fixed callback.
The normalized product norm is submultiplicative and a derivative costs
at most n+1 times the through-(n+1) norm. Kpm_(n+1) is Q98's FULL pressure
coefficient requested afresh at output n+1 and pole p. No pole monotonicity
or reuse of a smaller archive is assumed. Thus, conditionally,

    ||delta Etheta||_n <= epsilon^D zeta delta^-p T_n,
    ||delta Ez||_n <= epsilon^D zeta delta^-p A_n,
    ||Qcirc delta Ek||_n <=2 ||delta Ek||_n.              (3)

The graph derivative supplies the only retained loss kappa in this coarse
bound. Axial derivatives really carry epsilon; their favorable grades
C+1, and the radial connection grade C, are conservatively embedded into
D for 0<epsilon<=1. They are not physically replaced by unscaled
derivatives. Haar projection has norm at most two because zeta and the
edge pole weight are auxiliary-independent after Q97's Gaussian handling.
No weighted inequality is differentiated to infer a derivative inequality.
The native derivative contracts are separate inputs at every requested jet.

Whole-hull smooth compact extensions, mixed traces, current weighted
covariance norms, fixed graph norms and Q98 compact operator hypotheses
remain conditional. If traces do not match, Q98's cumulative jump and
Haar(J) pressure terms must be propagated through Dz too; (2) is not an
enclosure for discontinuous piecewise classical data without those terms.
The coefficients are symbolic positive exponential-polynomials, not actual
numerical continuum enclosures or a proof uniform over all physical scales.

## 3. Why the incoming final row cannot be skipped

Q47 FullRows.rows supplies the nonzero-Haar POSTSIGNED, PRE-NEW-FAST
row of the FIRST real mean cycle. The source state here has already
completed those fast and slow changes (Q86 FinalMeans). It is not that
earlier row. The exact retained first-cycle cancellation instead gives

    Qcirc E_FINAL_OLD,theta=Qcirc(Rfast_theta+Rslow_theta),
    Qcirc E_FINAL_OLD,z
       =Qcirc(A1 Qcirc Ez_POSTSIGNED+Rfast_z+Rslow_z).     (4)

Rfast/Rslow are the complete mean differences apart from c Nt of the
tangential/axial increment, which cancels the old row or vanishes for the
slow field. They retain slow time, viscosity, all correlations, full
recomputed T0 pressure, and actual A1-corrected axial velocities. The A1
source in (4) is POSTSIGNED; it is not Q47.initial_A1, whose source is
the original PRIMARY PREFast row. Neither tail may be transplanted into
the other stage. The frozen core_fast_row_gain_join and
core_post_signed_mean_closure already establish this cancellation calculus;
it is not newly claimed here.

The precise remaining sufficient inequality, at each needed jet j, is

    ||Qcirc E_FINAL_OLD,k||_j
          <=epsilon^D zeta delta^-pstar O_(k,j)(S),       (5)

derived from the actual terms of (4) with common-domain operator recipes,
not supplied as a desired-row callback. Then (1)--(3) give total changed
row coefficients O+2(T,A), after embedding the NEW increment's pole p
into pstar with its original finite geometry factor. For delta<=L and
pstar>=p use max(1,L)^(pstar-p); do not presume operator constants are
monotone in pole. The old first-real row grade 2-4kappa, when established
for the appropriate component, is stronger than D for 0<kappa<1/4.
That exponent comparison does not construct its coefficient or close its
compact axial source/support obligation.

Q99's scalar gain cannot prove (5). A compact off-diagonal comparison
h(R)cos(y) has zero auxiliary Haar, hence zero scalar moment changes, but
nonzero D2[h(R)cos(y)]. This is a tensor comparison, not a claim that it
is a realizable canonical covariance. Similarly a nonzero axial pressure
gradient is not erased by zero Haar. The tests retain both comparisons.
SourcedRowChanges deliberately refuses fresh_fast until the incoming
final row composition and its domain are supplied. Its output is explicitly
an increment majorant, not a bound on the entire changed row.

## 4. Requests and next obligation

Output n needs covariance through n+2, Q96 source through n+3, Q94 mean
factory through n+5, conservative target/profile n+19,n+32, and pressure
cutoff shape through n+2. This is one extra pressure output derivative
relative to Q99's integrated-moment calculation; no Haar primitive is used.
For a prospective fast output m, the inherited inverse reads rows through
m+5: the NEW increment alone then needs covariance m+7, mean factory m+10,
target/profile m+24,m+37 and pressure cutoff shape m+7. All older incoming
row/compact-tail DAGs remain additive and unexpanded. The original time
inverse S and e=1 KT1/KA1 domains will also be required, independently
of the e=0 pressure bound. No new accepted-band or width condition is set.

Next compose (4)--(5) at the required weighted poles, then add (3) and
apply the original complete fast inverse. Only after that may its full
old/new mean correlations and fresh slow response be recomputed. Actual
trace/native contracts, unconditional pulse bounds, full residual gains,
infinite realization and blow-up remain unproved. Original fields,
root, width C_setup^-20 and frozen historical evidence are unchanged.

Exact tests subtract independently written full momentum rows including
time and both viscous expressions, verify connections and the pressure
gradient, check two derivative orders of (2), and enforce the next pressure
jet, exact pole, stage and reserve. Manufactured norm leaves are interface
controls, not continuum data. Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json`.
