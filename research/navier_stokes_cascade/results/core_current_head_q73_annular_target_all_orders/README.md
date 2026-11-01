# Q73: an upstream angular bound for the earlier annular target

Q70--Q72 still left higher angular orders open beyond the short taper.
We give explicit (very coarse) constants for every finite mixed order
of the native PREMOD inputs from Xa through the earlier annular join.
Together with Q69 this removes the six/five angular restriction on
[Xa,Xrep]. This is conditional on the frozen core, accepted low-order
annular root and its inverse; no high-order contraction is assumed.

## Extend Q72 without extending the analytic seed in radius

Use exactly [Q72](../core_current_head_q72_short_taper_inputs/README.md)'s
order-m radial norm, constants b,t,B and strips D_s. Set tau=log(Xi/Xa),
Xi=110. On D_b the reference F, reciprocal F and U bounds of Q72
are unchanged on [Xa,Xi]: after B the fields F,U are constant in y,
with matching jets by flatness. This is not true for their moments.
In Q72 (2),(3) replace the primitive factor 2t+1 by tau+1 and the
geometric bound x by exp(1)Xi, retaining q=exp(1)/Xa and the original
five axis/core traces. This gives new source bounds Q1,N on D_(b/2).
This argument does not evaluate the native seed outside Y<=4.1.

The later switches beta,theta and their complements, as well as kappa,
have the same norm cap T_m as Q72 (unit values, width t). Globally on
this continuation the actual shear prescription is
a=(1-theta)kappa p1r+(4/5)theta, U_y=-kappa beta X nsr/2, with the
original plateau extensions. Hence put

    H=(tau+1)(T_m² Q1+(4/5)T_m)/2,
    fa=F0 exp(H), va=V0 exp(H), ua=U0+(tau+1)T_m² x N/2.

These bound F_a,1/F_a,U_a on D_(b/2) by integration of (log F_a)_y=-a/2
and U_y. Compute ACTUAL moments anew using Q72 (2) with fa,ua.
The second use of (3) on s=b/2 then encloses actual states on D_(b/4).
For b_s use the prescribed U_y bound and reciprocal E as in Q72.
No parameter derivative is repeatedly fed back into an unknown field.

## Explicit shape entry, with all angular orders

Write s=b/2, h=s/2, R=X_R. Use the preceding bounds at order zero here.
Let e=sqrt(2Xi)fa and v=va/sqrt(2Xi) at Xi, so |E|<=e and |1/E|<=v
on D_s. The nonvanishing E has a holomorphic log on this convex strip
with the original real branch. For ell=log(C E(Xi)) the real bound is
|ell|<=J6 by the native sixth estimate. Cauchy on D_h gives
|ell'|<=ev/(s-h). Integrating to the nearest real center thus gives

    |ell|<=L0=J6+ev,   |Gi|<=ua,   |g=Gi-4eta|<=ua+8   on D_h.

Both endpoint disks are included. This step avoids inferring bounds on
the imaginary part of a logarithm just from its modulus on a complex strip.
On D_h, |log f|<=4 for f=(1+eta²)^(-1), by Q68. On the shape use the
same T=Tsh, z=log(X/Xi) and (1-sigma(z/T))ell+sigma(z/T)log f
as Q70. The value |E|<=Es=exp(T/10+L0+4)/C follows there. In Q70 (1)
take H_j=1_{j=1}/10+(L0+4) B_j/(j! T^j), with B_j the Q63 raw
step bound. Its exponential recurrence with initial factor Es bounds the
radial derivatives uniformly on D_h. U=Gi is radially constant.
Angular Cauchy now bounds each normalized eta derivative by h^(-j)
times that radial bound, at every real center. The real E floor in Q70
is unchanged. This also bounds the data throughout the ideal extension
and width-one axial restoration: E=P f exp(y/10), |f|<=2 and
U=4eta+g(1-sigma(y+8)), with the same plateau extensions.

## Bound the target before solving for its correction

Let xc=exp(-6), Xc=R xc. Use only the ACTUAL prefix through Xc, before
the relative E and additive U bumps are applied. On D_h uniform value
caps there are

    Ebar=max(sqrt(2Xi)fa,Es,2P),  Ubar=ua+16.

In particular the fixed radial restoration creates no angular singularity.
Apply the five original X densities on [Xa,Xc] with X<=R and X>=Xa.
Starting with the five Q72 trace caps t_i at Xa, take as component caps

    t + R (Ubar, sqrt(2R) Ebar, sqrt(2R) Ebar Ubar,
           Ubar²+Ebar²/2, Ebar²/(2Xa)).

Call their maximum plus 200(1+P²)(1+R)² the number D. This additive
term bounds every IDEAL trace magnitude at Xc: substitute xc<1,
|eta|<=2 and |f|<=2 in the five ideal primitive formulas of the
[annular ledger](../axis_matching_audit/README.md), retaining R for M,S
and R^(3/2) for I,J. For ideal pressure include |Pi_axis|<=20P².
Consequently every actual-minus-ideal difference has modulus <=D,
also after removal of the positive R powers since R>=1. Delta Pi equals
Delta Cp only because both use the SAME axis datum.

The target is the negative normalized difference at Xc, after restoration.
Its numerators in hatted units are delta M, delta J-4eta delta I, delta I,
-(delta S-8eta delta M), delta Cp. Their magnitudes are <=17D.
The inverse denominators are <=25 xc^(-2): |f^-1|<=5, P>=1 and
the positive xc exponents are 1,8/5,8/5,6/5,1/5. Thus a single
uniform cap for the five-component target is Z=500 D xc^(-2) on D_h,
and ||partial_eta^k target||_infty/k! <= Z h^(-k) on the real interval.
This is a size bound; it is NOT the smallness estimate used to select
the root. No value or derivative of that unknown root entered Z.

## Pass to the same real root at every finite order

Use the [earlier eta-dependent quadratic map](../annular_exterior_foundation/README.md),
not the constant quadratic map of the later compact repair. In the
component maximum norm set d0=6e-13 (accepted root) and, for k>=1,

    d_k=2000 [Z h^(-k)+100 sum* d_j d_l],
    sum*: i+j+l=k, 0<=i<=4; omit (0,k,0),(0,0,k).

The frozen inverse bound is 2000. Each normalized coefficient derivative
of Q_eta is bounded by 100 and vanishes above degree four, because its
only varying scalar is proportional to (1+eta²)². The omitted terms
are exactly the Jacobian applied to the kth derivative. Induction gives
||c^(k)||_infty/k!<=d_k. In particular i>0 terms involving the
undifferentiated root must remain. Only real derivatives of c are used;
this argument does not assert its continuation to the entire target strip.

For a requested angular order l insert sum_(k=0)^l d_k for the coefficient
norm in Q70 (2), (ua+8) sum_(k=0)^l h^(-k) for the g norm, and
2 sum_(k=0)^l h^(-k) for the f norm in place of 28. The width-.08
relative swirl and additive axial bumps, and their radial factorial
bounds, do not change. Use the ORIGINAL small real coefficient bound
for the E floor, not the large high-derivative envelope. Q70 (3)
converts the log-radius bounds to X. For state order j take field and
trace order l=j+1; Q68's matched trace at R and Q70's partial integrals
then bound both pressure terms inside the bumps. This last backward use
of the match comes AFTER the target and root bounds, not before them.

This establishes explicit finite-order regional inputs, not smallness
uniform in order or the physical-band inequalities. Next propagate them
through the actual modulation and the normalized defect of the later
five-row repair at the fixed frequency, including its scale multipliers.
The physical pulse, pressure tails, complete covariance/residual improvement
and blow-up are not established by these estimates.

Run test_audit.py and audit.py --verify-record evidence.json --output a fresh
path. The finite algebra checks and provenance replay are not continuum proofs.
