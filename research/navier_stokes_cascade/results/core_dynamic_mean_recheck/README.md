# Dynamic current-mean budgets with old means and changing covariance

Parent: `6eebf5ce34436d4f29af466bd76928fc1937b9c4`.
The initial fixed-W projected lower bounds are preserved in
`core_covariance_mean_recheck`. The general current identities and conditional
cycle estimates already exist in `current_target_moment_audit`,
`current_target_cycle_bridge` and `complex_mean_relaxation`; this checkpoint
does not claim to discover or replace them. Its new result is an explicit
finite-jet upper budget for the dynamic nonlinear defects, a derivative budget
for the entire moving pressure moment, and sharp correlation controls showing
why a quadratic-only estimate needs additional cancellation hypotheses.

## 1. Exact current update after the linear rows

Use the same complete background (b,V,G,p_B), chart, fixed positive tilt and
common primary tail. On the actual reserved patch G=0 and V is the retained
power law; the identities below also allow G nonzero. Let the OLD actual mean
velocities be (beta,v,gamma), and write their Haar means as B0,v0,g0. Actual
A0/A1 velocities remain in these fields. Let the auxiliary-independent slow
increment be (B,u,g), compactly supported in I=[a,d], 0<a<d, with
D1 B+epsilon g_Z=0 and both zero velocity moments. Let
D_ab=bar(Wnew_ab-Wold_ab) for the COMPLETE angular covariance, averaged also
over the auxiliary variables. The six entries contain all complete-curl
cross terms and squares. This notation describes either a general composed
update or the covariance error left by using stale targets. In the prescribed
fresh slow substep W is held fixed and every D_ab is exactly zero.

For (1)--(3), I must contain the supports of ALL increments and covariance
differences appearing there. A fresh slow-only update can use its reserved
potential hull. A wave change on a wider radial shell requires a larger I
and its larger geometry constants; its support cannot be assigned to I4
merely because the slow map is supported there. The base G need only vanish
on the reserved increment support, and the formulas allow it elsewhere.

Assume the five linear rows use the incoming current P,Jtheta,Jz BEFORE the
displayed combined update, with the same base V,G. Put

    S=2(b+B0)B+B^2+Drr,
    H=(b+B0)g+(G+g0)B+Bg+Drz,
    U=2v0 u+u^2+Dtt,
    K=g0 u+v0 g+ug+Dtz,
    Zf=2g0 g+g^2+Dzz,
    tau=epsilon partial_T+epsilon^3 partial_ZZ.

Direct subtraction of the full current conservative equations gives, after
subtracting the linear centrifugal source 2Vu/R,

    Rg=tau B-D1 S-epsilon partial_Z H+U/R
       +epsilon(B_RR+B_R/R-B/R^2).

Integration by parts (both radial vector-viscosity moments vanish) then gives

    Pnew=tau integral B+integral (U-S)/R-epsilon partial_Z integral H,
    Qrem=tau integral R^2 B+integral R(S+U)
         -epsilon partial_Z integral R^2 H,
    Jtheta_new=integral R^2 K,
    Jz_new=integral R Zf-Qrem/2.                              (1)

All integrals are over the complete positive radial axis using smooth compact
extension. Haar averaging removes only the old zero-Haar parts multiplied by
these slow increments; B0,v0,g0 generally survive. For a fast increment the
products must instead be averaged AFTER multiplication, as in the retained
cycle's fast formula; (1) must not be applied to it by substituting zero
averages. These formulas do not assert that the actual incoming defects or
five-map outputs satisfy stronger bounds than the inherited hypotheses.

Drt is absent from these scalar moments, but remains in the target flux:

    delta Atheta=(b+B0)u+VB+Bv0+Bu+Drt
                  -epsilon(u_R-u/R),
    delta Az=H-epsilon g_R.

Together with the complete pressure change, these fluxes enter the exact
target difference -Ie Qe bar(delta Ee). A scalar defect estimate alone is
therefore not a bound for the complete signed target or residual.

## 2. Explicit finite-jet upper budget

Let Omega be a fixed slow domain. Choose a common radial-independent weight
0<omega(Z,T)<=1. For a field f, including auxiliary dependence, define

    |f|_(n,omega)=sum_{i+j<=n} sup_{I x Omega x torus}
                    |partial_Z^i partial_T^j f|/(i! j! omega).

The unweighted norm uses omega=1. For scalars omit I and the torus. These
norms obey |fg|_(n,omega)<=|f|_(n,omega)|g|_(n,omega), since omega^2<=omega;
an unweighted base factor can multiply a weighted field with the same bound.
The normalized Leibniz rule proves this, without a hidden 2^n factor.
Also |partial_Z f|_n,|partial_T f|_n <=(n+1)|f|_(n+1) and
|partial_ZZ f|_n<=(n+1)(n+2)|f|_(n+2). Haar averaging is a contraction.

Fix m>=0 and N=m+2. Let x,y,z bound the N-norms of B,u,g; let b*,G* bound
the UNWEIGHTED base norms and B*,v*,g* the weighted norms of B0,v0,g0.
Let d_ab bound the weighted N-norm of D_ab. Set

    A=b*+B*, C=G*+g*,
    s=2Ax+x^2+d_rr, h=Az+Cx+xz+d_rz,
    t=2v*y+y^2+d_tt, k=g*y+v*z+yz+d_tz,
    zeta=2g*z+z^2+d_zz,
    t_m=epsilon(m+1)+epsilon^3(m+1)(m+2).

Here g* is an old-mean norm, not the increment g. To make the geometry fully
explicit and rational when endpoints are rational, put
L=d-a and M_j=L max(a^j,d^j) for j=-1,0,1,2. Then integral_I R^j<=M_j,
including M_-1=L/a. Applying the preceding inequalities to (1) yields

    p_m=M_0 t_m x+M_-1(s+t)+epsilon(m+1)M_0 h,
    q_m=M_2 t_m x+M_1(s+t)+epsilon(m+1)M_2 h,
    jtheta_m=M_2 k,
    jz_m=M_1 zeta+q_m/2.                                    (2)

These respectively bound |Pnew|_m, |Qrem|_m, |Jtheta_new|_m and |Jz_new|_m.
They use only incoming field/increment/covariance norms, not an assumed bound
on the desired outgoing residual. Radial derivatives do not occur in these
moment budgets because their boundary terms have been proved to vanish.
They remain necessary in the full pointwise pressure and residual estimates.

In the initial zero-Haar, fixed-W case (2) bounds the dynamic C_P,C_Q and
g*u from the previous checkpoint by setting B*=v*=g*=d_ab=0 and separating
the u^2 and g^2 energies. In the current slow state these old-mean terms
cannot be set to zero. Its retained orders B0 in M_2, v0,g0 in M_1,
B in M_(alpha+1), u,g in M_alpha recover the conditional alpha+1 scalar
order of `complex_mean_relaxation`; no improved order is inferred here.

## 3. Whole pressure and a precise all-scale acceptance inequality

Use the same unit-mass moving pressure compensator rho and
c_rho=(integral R^2 rho)/2. The complete new axial residual moment is

    Mz=epsilon partial_Z(Jz_new+c_rho Pnew),
    Mtheta=epsilon partial_Z Jtheta_new.

If c_(m+1) bounds the unweighted norm of c_rho, (2) at order m+1 gives

    |Mz|_m <= epsilon(m+1)(jz_(m+1)+c_(m+1)p_(m+1)),
    |Mtheta|_m <= epsilon(m+1)jtheta_(m+1).                  (3)

Thus the displayed residual-moment bound uses input norms through m+3;
an m+2-only bound is not established here. The entire product is differentiated: its nonzero
c_rho,Z Pnew contribution is included. The full background Pi1 and Pi2,
their physical powers and cutoff derivatives are retained in the incoming
balance. They cancel only on subtraction of the SAME full base equation.
Neither (1) nor (3) licenses setting their actual axial gradients to zero.

For a proposed scalar tolerance L_m epsilon^rho on every accepted chart,
a sufficient CHECKABLE norm gate is

    max(p_m,jtheta_m,jz_m) <= L_m epsilon^rho,               (4)

with the additional right sides of (3) bounded by the proposed residual
moment tolerance. This is a sufficient upper-budget gate, not a necessary
condition: a failure can reflect cancellation lost in triangle inequalities.
It is not a complete pointwise residual or signed-target criterion.

On any fixed interior potential hull a comparison with the inherited zeta
norms is possible if sup(zeta/omega) is finite for that chart. Such a
comparison constant, the shell endpoints, all input norms and c_(m+1) must
be retained in (4). No uniform comparison over moving hulls or physical
scales is established here. The primary common tail is unchanged; obtaining
uniform actual input bounds so that (4) holds there is the remaining
quantitative obligation, not a consequence of passing a finite fixture.

## 4. Sharp correlation rate controls

The old-mean rate is substantive. Choose a nonzero smooth compact psi in I,
f=D1 psi, so integral R f=0. On an interior slow-time plateau take G=0,
beta=v=0, old gamma=epsilon^mu f and increment g=epsilon^alpha f, B=u=0,
independent of Z. Both states have exact divergence and both velocity moments.
Hold W fixed. The radial-source difference vanishes, and EXACTLY

    delta Jz=(2 epsilon^(mu+alpha)+epsilon^(2alpha)) integral R f^2.

For alpha>mu, a uniform O(epsilon^(2alpha)) bound would require
2 epsilon^(mu-alpha)+1<=constant, which is false. This is a continuum
countercontrol to deriving that stronger rate from support, divergence and
moments alone; it is not the prescribed actual five-map or a refutation of
its conditional alpha+1 estimate. Nonzero old Haar means require the cross
budget or an additional cancellation theorem.

For a changing COMPLETE wave w+q,

    delta W=<w tensor q+q tensor w+q tensor q>_theta.

Using the maximum component norm with sqrt(omega) for waves gives each covariance entry
the bound 2|w|_N |q|_N+|q|_N^2 in weight omega. This follows for every fixed
jet from Leibniz and averaging. Hence orders nu and chi yield only
min(nu+chi,2chi) without correlation cancellation. The sharp same-phase
control q=t w gives delta W=(2t+t^2)W, simultaneously for all six entries
of a generic complete curl. Different-phase orthogonality can remove a cross
in a special case; it cannot be silently imported for a generic signed update.
The retained cycle's subtraction of its prescribed principal response is
still needed before claiming a stronger signed error. Axial covariance and
pressure must then be recomputed, not supplied by the two matched entries.

## 5. Evidence and next obligation

Exact rational polynomial controls recompute two current states independently,
including nonzero old Haar means, a positive complete-wave covariance update,
all five scalar dynamic expressions, compact pressure and its moving moment.
Coefficient majorants on a bounded box verify the explicit norm budgets;
these are conservative analytic polynomial bounds, not sampled suprema.
Rational small-scale sequences detect the old-mean and same-phase wave rate
failures. The finite endpoint fixtures are not the actual C-infinity profile.
The continuum proofs are sections 1--4; tests are implementation diagnostics.

Next obtain uniform actual input norm/comparison constants for (2)--(4), and
use the retained current signed cancellation to control Drt and the full
pointwise pressure/residual increments. Generated harmonic supports, full
viscosity and cutoff interactions, full residual closure and blow-up remain
unproved. Historical frozen sources are unchanged. Run standard-library checks:

```sh
python3 research/navier_stokes_cascade/results/core_dynamic_mean_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_dynamic_mean_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_dynamic_mean_recheck/evidence.json \
  --output /tmp/core_dynamic_mean_recheck_fresh.json
```
