# Certified complete-potential field: localization, compact force and energy join

Parent `6a1e9459dd091595f471e904f0e997e241313ea3` gives compatible
one-sided terminal jets on punctured compact sets for the same certified
complete-potential diagonal field constructed at
`747ee32f11aa9b8336b0686b6f66a5affa31b1d7`.
This checkpoint joins that field to axisymmetric localization, zero initial
data, one smooth compact force and the preterminal energy estimate.

The implication remains conditional on the parent's all-fixed-order
physical bounds and uniform complete-residual flatness.  The exact controls
check algebra, derivative orders, endpoint gluing and the energy ledger;
they do not evaluate the infinite field or prove those continuum inputs.

## 1. Localize the same representatives

Retain the exact fixed-gauge decomposition

    u0=curl(A)+B,        B=b(r,z,t)e_theta,        p0=p_complete.       (1)

`A` contains the certified complete-potential diagonal's background, waves, means and every
cutoff-curl representative.  `B` is its separate theta-independent direct
swirl.  `p0` contains both background tails, every changing wave and
covariance pressure, and all axial and moving-moment memory.  Nothing is
regauged or reconstructed from velocity.

Keep

    q_e=min(qstar/16,2^(-ell0-2)),
    tau0=q_e/16,       z0=(q_e/16)^D.                       (2)

The retained cylinder has `q<=3q_e/8`.  The certified enlarged chart range
is `1/4<=q/Q<=4`, so

    4Q<=16q<=6q_e<=3qstar/8<qstar.                          (3)

Choose a Cartesian-smooth axisymmetric cutoff
`c=chi_r(r^2)chi_z(z)chi_t(t)` whose support closure and a neighborhood lie
strictly inside this cylinder.  It is one on a smaller late-time core and
zero for `tau>=tau0`.  Define the representatives on that neighborhood
before zero extension and put

    u=curl(cA)+cB=cu0+d,       d=grad(c) cross A,       p=cp0.          (4)

Both summands in (4) are exactly divergence free: curl is divergence free,
and `div(cB)=r^-1 partial_theta(cb)=0`.  Also
`d dot grad(c)=0` exactly.  A direct cutoff `cu0`, or a theta-dependent
cutoff on `B`, would not have these properties.  The time cutoff gives an
initial neighborhood on which `u=p=0`; hence `u(0)=0`.  On the smaller
late core, every derivative of `c-1` vanishes, so the certified complete-potential diagonal
field and its protected inner behavior are unchanged.

## 2. Complete residual after localization

For

    R_nu(v,p)=v_t+(v dot grad)v-nu Delta v+grad p,

direct expansion of (4) gives

    R_nu(u,p)= c R_nu(u0,p0)+c_t u0+d_t
      -nu(Delta c)u0-2nu sum_i c_i partial_i u0-nu Delta d
      +p0 grad c+c(c-1)(u0 dot grad)u0
      +c(u0 dot grad c)u0+c(u0 dot grad)d
      +c(d dot grad)u0+(d dot grad)d.                       (5)

These are twelve distinct retained terms.  The otherwise present
`(d dot grad c)u0` is zero by the cross product identity.  In particular,
`p0 grad c` contains the complete pressure memory.  The actual viscosity is
`nu=1`; changing it does not inherit the flat-residual theorem.

An order-`m` force jet needs `A` through order `m+3`, `B` through `m+2`
and `p0` through `m+1`.  Uniform terminal Cauchy convergence reserves one
additional time derivative of each.  Every requirement is fixed and finite,
so the parent's punctured endpoint theorem supplies it conditionally.

Set `f_minus=R_1(u,p)` for `t<1`.  Its exact divergence relation is

    -Delta p=partial_i partial_j(u_i u_j)-div(f_minus).     (6)

The force is not required to be divergence free.  Dropping its divergence,
using an unforced pressure formula, or invoking a compact-support-preserving
Leray projection would change the problem.

## 3. Exact terminal condition and blocker

Near `(0,1)`, `c=1` and `d=0`, so `f_minus=R_1(u0,p0)`.  The needed corner
estimate is, for every fixed multiindex `alpha`, time order `j` and `N`,

    sup_corner |D_x^alpha partial_t^j R_1(u0,p0)|
      <= C_(alpha,j,N) (tau+|z|^(1/D))^N.                  (7)

It must hold uniformly over the complete compact-`X` region, with exact
zero residual in the heat exterior.  Pathwise flatness is insufficient.
Away from the origin, the separate required condition is the parent's
uniform fixed-order bounds for `A,B,p0` and their next time derivatives on
each punctured terminal compact.

Equation (7), together with a small corner ball, the punctured compact
cover and zero exterior, gives uniform limits

    F_j=lim_(t up to 1) partial_t^j f_minus in C_c^infinity(K).         (8)

Spatial fundamental-theorem-of-calculus identities pass through the
uniform limits, including on segments through the origin.  The next time
derivative identifies consecutive time jets.  Every spatial derivative of
every `F_j` vanishes at the origin.  This regularizes the force trace only;
it gives no regular velocity or pressure at `(0,1)`.

If (7) fails, or if the punctured bounds fail through the accepted
`|eta|=1` or `T=0` chart boundaries, smooth global force jets do not follow.
These are the precise retained analytical blockers.

## 4. One future force without analytic continuation

Let `sigma=t-1>=0` and take a smooth `chi0` equal to one on `[0,1/2]` and
zero on `[1,infinity)`.  With `b0=2` and increasing dyadic `b_j`, define

    f_plus(x,1+sigma)=sum_(j>=0)
      chi0(b_j sigma) sigma^j F_j(x)/j!.                   (9)

The exact `m`-th time derivative of the scalar factor is

    sum_(0<=k<=m, m-k<=j) binomial(m,k) b_j^k
      chi0^(k)(b_j sigma) sigma^(j-m+k)/(j-m+k)!.          (10)

On the whole support `0<=sigma<=1/b_j`, each term is bounded by a constant
times `b_j^(m-j)`.  At stage `j`, impose the finite collection
`|alpha|+m<=floor(j/2)` and choose `b_j` large enough that each bound is at
most `2^-j`.  For every fixed mixed derivative, all sufficiently large
stages enter this one diagonal and the differentiated tail is uniformly
bounded by `2^(1-J)`.

Thus (9) is one compactly supported smooth force extension with terminal
jets `F_j`; it vanishes for `sigma>=1/2`.  Gluing it to `f_minus` yields,
conditionally,

    f in C_c^infinity(R^3 x (0,infinity)),
    supp(f) subset K x [1-tau0,3/2].                       (11)

No raw terminal Taylor series convergence is used.  The zero-radius
exterior heat and pressure Taylor series therefore cause no conflict.
Equation (9) constructs no future velocity or pressure.

## 5. Preterminal energy and a corrected endpoint claim

For every `t<1`, compact support and `div u=0` justify

    (1/2) E'(t)+D(t)=<f(t),u(t)>,
    E=||u||_2^2,       D=||grad u||_2^2.                   (12)

With `F(t)=integral_0^t ||f(s)||_2 ds`, regularizing the square root of
`E` at the zero initial datum gives `sqrt(E(t))<=F(t)`.  Substitution into
(12) gives the fixed-time estimate

    E(t)+2 integral_0^t D(s) ds <= F(t)^2.                 (13)

Consequently

    sup_(t<1) E(t)<=F(1)^2,
    integral_0^1 D(t)dt<=F(1)^2/2.                         (14)

The two left sides in (14) cannot be added while keeping the single budget
`F(1)^2`; their extrema need not occur together.  An exact scalar rise and
decay control with `epsilon=1/100` gives

    sup y^2+2 integral_0^1 y^2
      >=37349/22350 > 40401/40000=(integral_0^1 f)^2.      (15)

This corrects that tempting strengthening while preserving (13) and both
separate conclusions in (14).  The diagnostic is not a Navier--Stokes
solution.

This checkpoint establishes a **conditional certified-domain same-field localization,
compact-force and preterminal-energy implication**.  It does not discharge
(7), the primitive physical contracts or the protected-core premise.  It
does not provide a future solution, regularity at the terminal origin,
comparison with an arbitrary global competitor, or blow-up.

The next obligation is fixed-gauge pressure recovery and an expanding-ball
comparison with every admissible smooth finite-energy solution having the
same zero datum and the same force, without assuming competitor pressure,
derivative or boundary-flux bounds.  The retained primitive contracts must
also be discharged before any unconditional conclusion.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
