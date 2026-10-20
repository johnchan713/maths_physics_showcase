# Second-coefficient axis jets by nested analytic cores

Parent: `f2acf72d1abd5d68653ad25cddbe58507aa90040`.
Status: `nested-complex-core-second-axis-third-jets-and-conditional-inner-collar-join`.

The parent bounds the complete second inner source but obtains its third
radial jet only at positive radius. This checkpoint supplies the missing
AXIS-side estimate for the regular F2,U2, using a strictly nested complex
core. It joins that estimate to the parent's positive-collar result.
The result is conditional on the SAME retained infinite-core norm and
sparse Picard hypotheses. No independent unconditional foundation, finalized
global H2,3, retained-scale acceptance, physical pulse, covariance/mean,
complete residual improvement or blow-up result is claimed.

All roots, widths (including C^-20), frequencies, cutoffs, gauges, schedules,
supports and historical frozen evidence remain unchanged. In particular we
do NOT assume a flat activation or the smooth coefficient cutoff is radial
analytic. Only the analytic SEED is continued beyond the unchanged core,
and only its restriction to the core is identified with the actual fields.

## 1. A smaller disk supplies the full second source

Keep S>=2, 0<rho<=1, Lambda>0, 0<h<=.01, 0<tc<=.01,
delta=rho/(800S)<=1/1600,
Xa=4/Lambda. The previous first-axis checkpoint constructs the first
analytic seed solution on |X|<5/Lambda, angular strip 3delta/2, with

`|W1core| <= M1core=(z1+B1) exp(B1)`.

Here z1,B1 are EXACTLY the explicit rational factors of
`axis_first_coefficient.core_envelope`, not new small chosen constants.
Its first four components are even holomorphic in xi and hence holomorphic
in X=xi^2. The last two are odd. M1core is kept symbolic.

For the second source use |X|<9/(2Lambda). Disks of radius 1/(4Lambda)
around that disk lie in |X|<19/(4Lambda)<5/Lambda. For r<=3,

`|partial_X^r(F1core,U1core)| <= r! (4Lambda)^r M1core`.

On the angular strip delta/2, angular disks of radius delta/2 stay inside
delta<3delta/2. Thus a common RAW mixed bound through X order three and
eta order two is

\[
 H_{1c}=R_1 M_{1c},\qquad
 R_1=6(1+4\Lambda)^3(1+2/\delta+8/\delta^2).       \tag{1}
\]

The factor 8/delta^2 is 2! times (2/delta)^2; reusing a factorial angular
seminorm as a raw second-derivative bound would lose that factor.

The joint analytic leading reference norm bounds its X derivatives through
three on |X|<=5/Lambda, angular strip 4delta, by

`B3=2(256S+9+S/Lambda)(1+Lambda)^3`.

Indeed each r!20^-r(8/5)^(r+1), r=0,1,2,3, is less than two, and
Lambda^r<=(1+Lambda)^3. Angular Cauchy then supplies the EXTRA eta-third
jet needed by the radial previous-order axial viscosity:

\[
 H_{0c}=B_3(1+\delta^{-1}+2\delta^{-2}+6\delta^{-3}).\tag{2}
\]

The segment tX, 0<=t<=1, remains in the same complex disk, so every
regular average and its mixed derivatives has the same bounds. There is
no singular division by X and no need for a zero-free complex Phi.

Apply the parent's explicit positive pair map `operator_and_source` at

`x=9/(2Lambda), b=1+x, H0=H1=1`, obtaining q,qxi.

Its formulas are absolute triangle/product bounds and are equally valid
for complex X,xi. The enclosing b satisfies b^2>=x. They retain BOTH
ordered (0,1),(1,0) radial interactions, radial viscosity, the shifted
`-Z_(-1-D) Z_-1 v0`, both lower axial viscosities, and every current
pressure column. The current matrix uses lambda2=4h, lower sources use
lambda1=2h, and the current Pi2 exponent is -1+2h.

For the SCALAR, fixed envelope

`Tcore=max(1,H0c,R1*M1core)`,

every operator/source row is <=q Tcore^2 on the outer eta strip delta/2.
The xi derivatives are <=qxi Tcore^2 as well, though the value estimate
alone suffices for the complex-ray existence bound. This does not set the
first coefficient to a small diagnostic profile, or replace M1core by one.

## 2. Second analytic solution, parity and axis jets

Set delta2=delta/4, so the second outer source strip is 2delta2. The
six-variable system has the SAME diagonal D0=diag(0,0,2,0,3,1) and the
same sparse derivative-block pattern at current order two. The regular
complex-ray inverse is still

`(Gg)_i(xi)=xi integral_0^1 t^(d_i) g_i(t xi) dt`.

The first-axis checkpoint's ordered-simplex/sparse-word argument therefore
applies on complex rays. This invokes its retained analytical lemma, not
a new inference from finite Taylor tests or pointwise A1^2=0. The variable
diagonal kernels and angular derivatives of the right A1 preserve the
annihilation. With loss delta2/2, put

\[
 z_c=2bqT_c^2,\quad B_c=6z_c^2/\delta_2,\quad
 M_{2c}=(z_c+B_c)\exp(B_c).                         \tag{3}
\]

This bounds W2core on |X|<9/(2Lambda), angular strip 3delta2/2.
Bc grows as Tcore^4; Tcore already includes M1core. Thus (3) is a
NESTED exponential estimate, not the first-order bound relabelled, and
not a value at Tcore=1 scaled only after exponentiation. All such enormous
exponentials remain symbolic in the implementation.

For P0=diag(1,1,1,1,-1,-1), the literal current matrices satisfy
`A_j(-xi)=-P0 A_j(xi) P0`. The full source is

`f2=(0,0,0,2xi P2,2Stheta2,2[Sz2-2eta X P2/L])`.

P2,Stheta2,Sz2 are holomorphic functions of X on this core. Hence
`f2(-xi)=-P0 f2(xi)`. Analytic-class uniqueness with unchanged ZERO axis
data gives W2core(-xi)=P0 W2core(xi). F2,U2,K2,Pi2 therefore descend to
holomorphic X; the last two components remain odd. Holomorphy in xi
alone, without this parity, would not justify X derivatives at the axis.

At real 0<=X<=Xa, disks of radius 1/(4Lambda) stay in
|X'|<=17/(4Lambda)<9/(2Lambda). Consequently, for r<=3,

\[
 |\partial_X^r(F_{2c},U_{2c})|\le r!(4\Lambda)^r M_{2c}.
                                                               \tag{4}
\]

There is NO inverse-radius loss at X=0. Angular disks of radius delta2/2
give the C_eta^1 factor c2=1+2/delta2. The same factor is valid on the
smaller complex eta strip delta2/2, since those disks stay inside delta2,
strictly inside the solution strip 3delta2/2.

## 3. Identify only the unchanged core, then join the real collar

This identification must use finalized LOWER coefficients. The unchanged
coefficient cutoff is one through Xkeep=Xa exp(tc/8)>Xa, and all I3
repair bumps lie strictly farther right. Thus on 0<=X<=Xa the actual
finalized first coefficient equals its uncut first analytic core solution.
Its forward regular average and pressure reconstruction from the ORIGINAL
axis datum agree there too. Similarly the leading actual core equals its
analytic seed there. Therefore the actual second inner solution satisfies
the SAME order-two Volterra problem on that real core, and equals the
analytic seed restriction by uniqueness. No backward pressure influence
from the future repair is inserted. Conversely, this local identification
does NOT erase the pressure memory in later global moments or set Pi1/Pi2
to zero on the reserved interval I4.

For the positive region use the parent's SECOND-source estimate on the
actual activated interval, with its actual cutoff-dependent Tcollar. It
has a DIFFERENT Picard bound M2collar. Take a=min(1,Xa), so a^2<=Xa;
let (A_r,B_r) be its C_eta^1 affine constants computed at qcol,qcolxi,
delta2 and a, and p=(0,0,2,4). The two domains overlap, giving

\[
 \boxed{J_{2,r}\le c_2r!(4\Lambda)^r M_{2c}
       +T_{col}^{p_r}(A_rM_{2col}+B_r),\quad r\le3.}              \tag{5}
\]

This bounds the PRE-CUT second inner F2,U2 through X order three and
angular order one on the whole real interval 0<=X<=Xa exp(tc/4), and on
the common smaller eta strip delta2/2. Any supplied x,b for that collar
must enclose this actual interval, not merely the shorter analytic core.
The implementation requires x>=4.1/Lambda and b^2>=x. This encloses the
entire retained collar since exp(tc/4)<=exp(1/400)<41/40; the last strict
inequality follows already from exp(u)<=1/(1-u) at u=1/400.
`joined_jets` returns both independent source/exponential descriptions.
It never identifies M2core with M2collar or Tcore with Tcollar.

For the SAME eta-independent second cutoff the exact binomial rule gives
`sum_(j=0)^r binom(r,j) k_j J2_(r-j)`. `cut_rows` records the resulting
coefficients of M2core,M2collar,1 when a valid Tcollar is supplied. In
particular k3 J2_0 remains. Actual k1,k2,k3 have not been evaluated, and
the later exact repair's angular norms and global signed moment data are
not supplied by (5). Therefore (5) is not finalized/global H2,3.

## 4. Evidence and next obligation

Exact rational tests cover nested radial/angular margins, the third leading
raw jet, factorial conversion, both independent exponential envelopes,
cutoff product propagation, every current matrix reflection entry, sparse
diagonal kernels and full source parity. Nine complex manufactured source
samples independently evaluate the parent's full symbolic source against
its complex absolute bounds. Their nonzero axis pressure-source control
prevents a spurious zero quotient convention. Cancellation also verifies
that all four manufactured regular source expressions have no X denominator.
These floating samples are
not interval certificates, actual enormous selected fields or a continuum
proof. The analytical derivation above, not sample passage, supplies (5).

Next quantify the unchanged cutoff jets, global first U1 third-angular and
leading U0 fourth-angular inputs, and the complete second-repair moments
through Xv. Combine them with (5) for finalized H2,3 and retained-scale
acceptance. Complete physical pulse derivatives/supports, full-wave
covariance with axial velocity and pressure means, all viscous/cutoff/
interaction residual terms, and independent inherited foundations remain
open. Neither local core identification nor disjoint repair truncates the
full future pressure ledger.

```
python test_audit.py -v
python audit.py --verify-record evidence.json --output fresh.json
```

The output must be fresh. Historical evidence is unchanged.
