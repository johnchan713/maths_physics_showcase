# Q123: fixed repair bumps to the current retained diagonal, without giant integers

Parent: `91af8ae76d2b5bd82044792ac4d6d34af1f31e21` (Q122).
On 3 November 2026 the freshly fetched branch agreed with that parent.
Q122 research CI 41622757377 and field CI 41622757388 were in progress;
Q120 field CI 41608567714 had passed. The latest notes, frozen evidence,
Q63's actual-bump bound and the higher-cutoff-jet interface were inspected.
The core_growth_activation conclusion remains local leading core/first
activation collar, with the unchanged width C_setup^-20. Historical frozen
evidence, physical scales, roots, supports and pressure data are preserved.

Q122 exposed fixed bump radial norms as inputs. Those bump shapes were
ALREADY bounded at every order in Q63 and cutoff_higher_jet_interface.
This checkpoint does not claim a new proof of that result. It derives a
closed bound including the shifted E-to-F radial factor, and composes the
existing actual-bump estimate with Q122, Q121 and Q120 using only binary
budgets. This removes the independent bump-norm callbacks from that current
positive-order diagonal interface. Actual cut-seed and full-pressure moment
bounds, identity/domain bindings and existing scale intervals remain inputs.

## 1. Same five bumps, same complete logarithmic Jacobian

Use the retained dy-unit-mass squared-flat bumps

    beta(y)=10 sigma'(10(y-center)+1/2),
    sigma(s)=g(s)/(g(s)+g(1-s)), g(s)=exp(-1/s^2),
    y=log(X/rho3).

The five supports are the original disjoint width-1/10 intervals, ending
by y=2.55. No Borel, variance, or compact-X cutoff is substituted. Flat
endpoints retain the derivative support of each original bump.
Q63 proves the uniform translated logarithmic derivative bounds

    ||beta^(j)|| <= B_j=162*32600^(j+1)*((j+1)!)^3.   (1)

The older sharper bounds B0=320, B1=121600, B2=220000000,
B3=18870643200000 remain available and are used at those orders here.
The combined sequence is increasing, including the switch to (1) at j=4.
The inherited exact Bell/flat-step arguments establish these as bounds on
the fixed function, not numerical samples or arbitrary smooth-bump bounds.

Let p>0 be a lower X bound on all actual bump supports, and rp>0 a lower
sqrt(2X) bound there, with rp^2<=2p. For alpha=0 or 1/2 the full identity is

    partial_X^r[X^(-alpha)b(y)]
      =X^(-r-alpha) product_(i=0)^(r-1)(D_y-alpha-i)b(y). (2)

Its coefficients have alternating signs. Their absolute sum is therefore
exactly product_(i=0)^(r-1)(1+alpha+i). Since
1+alpha+i<=(1+alpha)(i+1), the sum is at most (1+alpha)^r r!.
For alpha=1/2 this retains all inverse-radius terms; replacing (2) by
X^-r D_y^r b would lose them. The additional constant 1/sqrt(2) for
E/sqrt(2X) combines with X^-1/2 into the stated rp^-1 bound.

Put a=max(1,p^-1). Positivity, B_j<=B_m for j<=m and (2) give the
THROUGH-order radial bounds

    B_U,m = a^m m! B_m,
    B_F,m = rp^-1 (3a/2)^m m! B_m.                  (3)

Here B_F bounds derivatives of beta/sqrt(2X), not of beta alone. The same
bound applies to each translated bump. Because supports are disjoint, no
factor five enters Q122's mixed derivative repair estimate. The old exact
low-order radial formulas remain sharper and valid; (3) is a conservative
closed majorant, not an assertion of unchanged numerical selector constants.

## 2. Binary budgets for every requested fixed order

No factorial in (1) or (3) needs to be materialized. For m>=4, define

    b_m=8+15(m+1)+3(m+1)ceil(log2(m+1)).             (4)

Then B_m<=2^b_m, since 162<=2^8, 32600<2^15 and j!<=j^j.
For m<=3 use ceil(log2 of the inherited sharper B_m). Write
f_m=m ceil(log2 m), with f_0=0. Sufficient integer bit budgets for (3) are

    u_m=b_m+f_m+m ceil(log2 a),
    v_m=b_m+f_m+m ceil(log2(3a/2))
                  +ceil(log2 max(1,rp^-1)).        (5)

They bound all radial orders through m and include the original physical
support lower bounds. They cannot be applied on a guessed positive hull
or to an unidentified bump. The API checks numerical geometry consistency;
actual support and function identity are still analytical premises.

## 3. Complete current join to the existing a_n condition

For a positive coefficient n and tuple request m, Q122 needs all five
full incoming pressure/moment discrepancy raw bounds through m+1. With
these and the original patch identity, its explicit coefficient bound
C_(m+1) is retained unchanged. Let c be ceil(log2 max(1,C_(m+1))). If

    Fcut through m <=2^s_F, Ucut through m+1 <=2^s_U,

then Q122's disjoint-support product estimate implies

    Fraw_m <=2^f, f=1+max(s_F,v_m+c),
    Uraw_(m+1) <=2^u, u=1+max(s_U,u_(m+1)+c).       (6)

These include the ENTIRE U primitive interval [0,Xhi]. Moment bounds
must retain the full source integrals through Xv, all lower-order products,
radial memory and viscosity. No local velocity plateau resets a discrepancy.

Q121's tuple prefactor is T=T_F Fraw+T_U Uraw, where

    T_F=4Rhi,
    T_U=2+(112+32(D+2nh))*max(1,Xhi)/((1-2h)Rlo),
    tau=min(1,Rlo^2)/8, D=1/2-h.

Thus define

    t=1+max(ceil(log2 max(1,T_F))+f,
             ceil(log2 max(1,T_U))+u),
    b_tuple=t+f_m+m ceil(log2(1/tau)).              (7)

This bounds Q121's full raw tuple without evaluating m!. Finally Q120 gives
the explicit sufficient condition on the SAME supplied existing scale:

    nh ell_n >= n+1+k_(m+1)
                  +m ceil(log2(1+2nh))+b_tuple,    (8)

where ell_n<=log2(a_n), and k is Q120's BOREL Euler cutoff budget.
The Borel cutoff and repair bump are distinct functions despite both using
squared-flat formulas. Each needs its own exact identity binding.

current_row implements (4)--(8) for 0<=m<=n with no factorial or high-power
evaluation. It chooses no scale. The manufactured passing-row test supplies
a diagnostic interval exactly meeting the sufficient inequality; it is not
an evaluation or reselection of an actual a_n. Failure of (8) leaves the
sufficient test unresolved, and is not proof of physical failure. Finite
heads n<m still need Q121's ordinary finite-head estimate; no diagonal row
outside its original range is inferred. Original weighted-stress, support,
normal-floor and selector gates also remain necessary.

## 4. Validation and next obligation

Exact polynomial checks retain every coefficient of the shifted radial
operator and compare its coefficient norm to (3). Independent symbolic
X differentiation checks (2). The inherited actual low-order bump bounds
and radial formulas are compared with the conservative through-order bounds.
Moderate-order exact Q122/Q121 values are bounded by the composed binary
budget. At order 100005 the diagnostic keeps factorials unevaluated and
requests discrepancies through 100006. Those source values are manufactured,
not evidence for actual giant-profile inequalities. Tests and source-hash
replay do not replace the analytical implications above or a continuum proof.

Next bound the actual cut-seed jets and all five full-pressure incoming
moment discrepancies at the required orders and establish their domain and
identity bindings. With the fixed bump constants now wired in, verify (8)
against the retained schedule. Then carry actual constants through the
current native high-pole source and compact-modulation/exact-repair cone join.
Global pressure/moment and reserved-mean support identities, complete-wave
mean corrections including axial velocity and pressure, uniform actual
pulses at every physical scale and the complete viscous/cutoff/interaction
residual remain open. No unconditional pulse bound, residual improvement
or blow-up follows.
