# Coupled current pressure and axial-defect budgets

Parent: `6d5262d06ad6e1c8fa88e3eb4a3f6abc144f8a73`.
The parent gave a quantitative angular scalar recurrence while retaining
current P,Jz bounds as inputs to the slow inverse. This checkpoint supplies
the corresponding pressure and axial scalar recurrences and their full
moving-pressure residual moment. The known postslow C+1 order is not new;
the new result is its explicit coupled polynomial budget, current covariance
and fast-correlation routing, and the input derivative ledger.

Everything is conditional on the retained canonical seed and current wave,
fast/slow inverse and edge estimates. No actual field constants, historical
first band, complete pointwise residual or infinite correction are certified.
Historical frozen evidence, profiles, width, tilt and the common domain are
unchanged.

## 1. Complete wave changes supply both preslow scalar channels

Use the normalized slow norms of the parent. On a positive integration
interval I=[a,d] containing the WHOLE change, put

    M_j=(d-a)max(a^j,d^j), j=-1,0,1,2.

For the Haar change D_ab of the complete covariance, independent integration
of the current conservative source gives exactly

    delta P=int (Dtt-Drr)/R - epsilon partial_Z int Drz,
    delta Jz=int R[Dzz-(Drr+Dtt)/2]
                          +(epsilon/2)partial_Z int R^2 Drz. (1)

This keeps axial velocity, anisotropy and the radial-axial covariance.
Drt and Dtz stay in the angular target/scalar equations; their absence from
(1) is not their deletion from the full residual. The wave does not change
mean velocities, so no old-mean term is being set to zero here.

At each of the two wave operations let the actual complete old wave and
increment satisfy the parent bounds epsilon^(1/2)w_n(S) and
epsilon^(C-1/2)q_n(S). Their complete tensor identity gives for C>=1

    |D_ab|_n <= epsilon^C d_n(S), d_n=2w_n q_n+q_n^2.      (2)

The old wave in the second operation is the updated one. From (1),

    |delta P|_m <= epsilon^C[2M_-1 d_m+(m+1)M_0 d_(m+1)],
    |delta Jz|_m <= epsilon^C[2M_1 d_m
                                      +(m+1)M_2 d_(m+1)/2]. (3)

The extra input order m+1 in Drz must be retained even though the axial
derivative has a small epsilon coefficient. No radial derivative is lost
in these integrated identities. A wider wave shell uses its own M_j;
it cannot be assigned to the reserved slow hull.

## 2. Fast pressure and axial changes keep the actual radial correlations

Hold W fixed during the fast step. Write (beta,v,gamma) for old means and
(B,u,g) for the ACTUAL zero-Haar increment, including its curl remainder.
Define averaged products AFTER multiplication:

    S_f=bar(2beta B+B^2),
    H_f=bar(beta g+gamma B+Bg),
    U_f=bar(2v u+u^2), Z_f=bar(2gamma g+g^2).

Only linear Haar terms vanish. Exact changes are

    delta P=int(U_f-S_f)/R-epsilon partial_Z int H_f,
    delta Q=int R(S_f+U_f)-epsilon partial_Z int R^2 H_f,
    delta Jz=int R Z_f-delta Q/2.                         (4)

Let mu=1-kappa, H=C-kappa, C>=C0=2-2kappa. At input order n=m+1 assume
old tangential/radial norms <=epsilon^mu O_t, epsilon^(mu+1)O_r and fast
norms <=epsilon^H F_t, epsilon^(H+1)F_r. All four majorants are nonnegative
polynomials in S. Set

    s_f=2O_r F_r+F_r^2,
    h_f=O_r F_t+O_t F_r+F_r F_t,
    u_f=z_f=2O_t F_t+F_t^2.

Every original exponent is at least C, since mu+H-C=1-2kappa>0 and
2H-C=C-2kappa>0; the radial factors only strengthen it. Consequently

    pi_f=M_-1(s_f+u_f)+(m+1)M_0 h_f,
    q_f=M_1(s_f+u_f)+(m+1)M_2 h_f,
    zeta_f=M_1 z_f+q_f/2                                 (5)

bound |delta P|_m/epsilon^C and |delta Jz|_m/epsilon^C, respectively.
The input order is m+1 for all factors in this conservative budget.

Add (3) for both current wave operations and (5) to the incoming scalar
polynomials. This gives explicit preslow Pi_P,m and Pi_z,m with

    |P_postfast|_m<=epsilon^C Pi_P,m(S),
    |Jz_postfast|_m<=epsilon^C Pi_z,m(S).                  (6)

Together with the parent's Pi_theta,m, these are all three fresh scalar
inputs. They have been bounded through actual current differences, not
inferred from the projected signed target.

## 3. Fresh slow output: explicit postslow polynomials

On the reserved potential hull use G=0, keep b nonzero, W fixed and the
CURRENT three linear cancellation rows. Let N=m+2. Retain the sharper
Haar rates and current inverse estimates in the form

    |b|_N<=epsilon b_N, |bar(beta)|_N<=epsilon^2 r_N,
    |bar(v)|_N,|bar(gamma)|_N<=epsilon a_N,
    |B|_N<=epsilon^(C+1) X_N,
    |u|_N<=epsilon^C U_N, |g|_N<=epsilon^C Z_N.            (7)

The majorants on the right are polynomials in S; notation Z_N is a norm
majorant, not the coordinate. On this hull define

    s=2(b_N+r_N)X_N+X_N^2,
    h=(b_N+r_N)Z_N+a_N X_N+X_N Z_N,
    t=2a_N U_N+U_N^2, z=2a_N Z_N+Z_N^2,
    t_m=(m+1)+(m+1)(m+2).

All geometry constants below refer to the reserved hull, not the wave
shell. The exact postslow moments from the dynamic mean identities yield

    Theta_P,m=M_0 t_m X_N+M_-1(s+t)+(m+1)M_0 h,
    Theta_Q,m=M_2 t_m X_N+M_1(s+t)+(m+1)M_2 h,
    Theta_z,m=M_1 z+Theta_Q,m/2.                          (8)

For example the time term epsilon partial_T int B has exponent C+2,
axial viscosity epsilon^3 partial_ZZ int B has exponent C+4, old angular
crosses have exponent C+1, and the new tangential squares have exponent
2C>=C+1. The nonlinear radial/base products and their axial gradients are
no worse. Dropping their extra epsilon factors using epsilon<=1 gives

    |P_new|_m<=epsilon^(C+1)Theta_P,m(S),
    |Jz_new|_m<=epsilon^(C+1)Theta_z,m(S).                 (9)

Radial vector viscosity is retained in the pointwise source; its two
integrated boundary terms vanish by smooth compact extension. That fact
does not remove axial viscosity or prove a pointwise residual bound.

The inverse coefficients make (8) an explicit coupling of the three scalar
ledgers. At each fixed N, the retained tangential map has finite majorants

    U_N<=L_u,N(Pi_P,N+Pi_z,N),
    Z_N<=L_g,N Pi_theta,N,
    X_N<=L_B,N Pi_theta,N+1.                             (10)

The L majorants keep their profile, normalization, a^-1 and fixed 1/lambda
costs; they are inputs from the actual inverse estimates. The curl defining
B accounts for N+1. Thus scalar output order m can require incoming target
order m+3. Equations (8)--(10), with the parent's angular output, supply a
finite coupled recurrence at every fixed order. Since C+1 exceeds
Cnext=C+1/2-kappa, all three outputs provide next-cycle scalar inputs.
This does not replace the separate projected-target or regular-source proof.

## 4. Whole moving pressure and the four-extra-target-jet ledger

With c_rho=(1/2)int R^2 rho, the full axial residual moment is

    Mz_new=epsilon partial_Z(Jz_new+c_rho P_new).

For the retained c_(m+1)=|c_rho|_(m+1), (9) at order m+1 gives

    |Mz_new|_m<=epsilon^(C+2)(m+1)
                  [Theta_z,m+1+c_(m+1)Theta_P,m+1].       (11)

This differentiates the WHOLE moving product, including c_rho,Z P_new.
In field variables it uses (7) through m+3. After (10), the radial curl
can require angular target jets through m+4. An m+3-only TARGET budget
for (11) is not established; m+3 is the field-input order. Finite smooth
induction can read the necessary higher jets. A same-order contraction,
uniform derivative ladder or infinite-cycle estimate does not follow.

Both background pressure memories Pi1/Pi2, complete axial velocity and
covariance, and actual A0/A1 terms remain in the same base/current fields.
These scalar identities cancel the same complete base equation only.
They do not set an exceptional pressure tail or its physical gradient to
zero. All polynomial losses in (6)--(11) remain; the parent's positive-
margin rational tail test can check actual tolerances on an already
accepted tail, but is not itself certification of those inputs.

Next independently verify the actual canonical seed, inverse and current
field/source/edge hypotheses, then join the complete physical residual.
The scalar recurrence is now quantitatively coupled under those hypotheses;
no unconditional pulse, residual improvement or blow-up is proved.

## 5. Evidence

Exact complete-current covariance changes verify (1) before applying norm
bounds. Nonzero-Haar old fields and actual zero-Haar increment products
verify every fast term in (4); radial correlations are not replaced by
products of averages. Independent radial-source integration verifies the
fresh slow remainders and the compact moving-pressure moment. Rational
positive-polynomial budgets and exponent controls compare (8) with the
unreduced dynamic budget. Diagnostics have finite endpoint regularity and
are not actual global fields or continuum supremum certificates. The
conditional estimates follow analytically from sections 1--4.

```sh
python3 research/navier_stokes_cascade/results/core_coupled_defect_budget/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_coupled_defect_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_coupled_defect_budget/evidence.json \
  --output /tmp/core_coupled_defect_budget_fresh.json
```
