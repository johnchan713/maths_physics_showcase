# Explicit fast inverse, actual axial repair and current-residual jet costs

Parent: `2adeca6b61022d0a952989c031d46bb94f741c03`.
The slow-map constants are now explicit. This checkpoint supplies fast-map
constants in terms of the CURRENT zero-Haar residual, including its actual
A1 axial velocity and the exact fast-time error. Qualitative compact and
temporal inverse existence and flatness were already established in
`compact_mean_foundation`; those results are not claimed as new.

The new result is the quantitative substitution into the coupled scalar
budget and the explicit input-jet ledger. No canonical residual norm,
core floor, complete pointwise residual or infinite realization is certified.
Profiles, width, tilt, carriers, historical evidence and common tail remain
unchanged. No numerical evaluation of the selected enormous fields is used.

## 1. A temporal bound with an explicit Fourier constant

Use the common torus y in R^2/Z^2 and actual directions
vr=(1,-b), vt=(b,1), b=sqrt(2)-1. Work in U=R^d, d=d_r>0 fixed.
Let N_m(f) be the maximum supremum of all RAW coefficient derivatives
in (U,Z,T,y1,y2) through total order m. These derivatives hold M,c fixed.
This norm is not the parent's normalized slow sum norm; its conversion
factor to that norm is J_m=(m+1)(m+2)/2. Fixed changes between R and U
have finite profile-dependent jet constants on the positive radial hull.

For k!=0 and r=|k|_1, the integer quadratic norm gives

    |vr.k|^-1, |vt.k|^-1 <= (1+sqrt(2))r < 3r.

If f has zero auxiliary Haar mean, the real zero-Haar inverse satisfies

    N_m(c^-1 Nt^-1 f) <= c^-1 Gamma_m N_(m+4)(f),
    Gamma_m=24*2^(m+4).                                 (1)

To verify (1), for output torus order a choose |k_j|>=r/2 and integrate
the Fourier coefficient by parts a+4 times in y_j. Input total order is
at most m+4. Output powers and the inverse divisor leave r^-3, with
constant at most 3*2^(m+4) after dropping inverse powers of 2pi. There
are 4r lattice points of norm r, and sum r^-2<=2, so the lattice sum
is at most 8. This proves absolute convergence and differentiation.
Reality follows by conjugate symmetry. On the same accepted band
c^-1<Tg*S; (1) has one explicit S factor and no epsilon loss.

## 2. Compact operators: coefficient jets have no M power

On the fixed positive U integration hull [u0,u1] set L=u1-u0. The SAME
moving cutoff chi(U,Z,T) is independent of y. Define

    f_e(U)=R^(e-d+1)/d, g_e(U)=R^-e,
    k_e(U)=d R^(d-1-e), F=f_e f,
    I_-=int_(-infinity)^0 F(U+a,y+Ma vr)da,
    I_+=int_0^infinity F(U+a,y+Ma vr)da,
    Te f=g_e[(1-chi)I_--chi I_+],
    Ae f=k_e chi_U int_all F(U+a,y+Ma vr)da.             (2)

Extend the input smoothly by zero. At fixed a the shift is independent
of U,Z,T, so differentiation hits only F and the displayed multipliers.
No M factor is introduced. Each integral has effective length at most L.
For N_m(fg)<=2^m N_m(f)N_m(g), valid sufficient constants are

    KT_e,m=2^(3m)L N_m(g_e)[N_m(1-chi)+N_m(chi)]N_m(f_e),
    KA_e,m=2^(3m)L N_m(k_e)N_m(chi_U)N_m(f_e).           (3)

Thus N_m(Te f)<=KT_e,m N_m(f) and N_m(Ae f)<=KA_e,m N_m(f).
The cutoff cost includes order m+1 of chi in the second bound. All
these multiplier seminorms are finite and uniform in the band by the
fixed-profile composition and positive-hull estimates already proved.
Their dependence on the actual hull, d, profile and normalization stays.
This argument bounds coefficient derivatives, not graph or Cartesian
derivatives; those carry the separately retained physical losses.

Both operators preserve the slow projection of support and have smooth zero
extension on the full fixed radial potential hull. They can fill radial
gaps and the auxiliary torus; A1 can live on the cutoff's radial support.
There is no rectangular weight-floor assertion near moving shell edges.

## 3. Bound the actual fast velocity without deleting A1

Let ftheta=Etheta^circ and fz=Ez^circ be the CURRENT angular-mean
residuals after removing their auxiliary Haar averages. Keep all wave,
base pressure, covariance and existing-mean contributions inside these
inputs. Set

    u=-c^-1 Nt^-1 ftheta, gd=-c^-1 Nt^-1 fz,
    psi=T1 gd, B=-epsilon partial_Z psi,
    g=D1 psi=gd-A1 gd.                                  (4)

The divergence identity is exact, including derivatives of the moving
cutoff. All increments have zero Haar mean, so both velocity moments are
preserved. The actual g in (4) obeys

    N_m(u) <= Gamma_m Tg S N_(m+4)(ftheta),
    N_m(g) <= (1+KA_1,m)Gamma_m Tg S N_(m+4)(fz),
    N_m(B) <= epsilon KT_1,m+1 Gamma_(m+1) Tg S
                                                   N_(m+5)(fz). (5)

The KT order in the radial bound is m+1. Bounding D1 psi naively by a
graph derivative would introduce an unnecessary M loss; the exact
gd-A1 gd identity gives (5) while retaining the entire correction.

If N_n(ftheta),N_n(fz)<=epsilon^H Etheta_n(S),epsilon^H Ez_n(S),
the parent's slow-norm fast tangential/radial majorants can be taken as

    Ft_m=J_m Gamma_m Tg S
                    [Etheta_(m+4)+(1+KA_1,m)Ez_(m+4)],
    Fr_m=J_m KT_1,m+1 Gamma_(m+1) Tg S Ez_(m+5).         (6)

These bound the actual (u,g) by epsilon^H Ft_m and B by
epsilon^(H+1) Fr_m. Constants in (3) may be replaced by their positive
polynomial S majorants if that is the supplied input format.

Substitute (6) at order m+1 into the coupled pressure/axial fast budget.
It requires ftheta through m+5 and fz through m+6. The angular scalar
change int R^2 bar(gamma_old u+v_old g+ug) has no additional slow
derivative and needs fast residuals only through m+4 at output order m.
Consequently the subsequent slow scalar output of order m can use fast
residuals through m+8; its full axial residual moment can use m+9.
Indeed the tangential slow response reads P,Jz at m+2 (cost m+8),
and its radial response reads Jtheta at m+3 (cost m+7). Increasing output
order once for epsilon partial_Z(Jz+c_rho P) gives m+9. These counts
concern the prefast residual inputs only, holding the other current-wave
and old-field hypotheses; they are sufficient, not optimal. The earlier
m+4 INCOMING SCALAR target cost is a different ledger and stays valid.

## 4. The fast-time error is an operator on the residual itself

Since chi and f_e,g_e,k_e are independent of y, Nt commutes with Te,Ae.
The chart coefficients c,M are constant during these derivatives. Thus

    c Nt g=-fz+A1 fz.                                  (7)

No temporal-inverse derivative loss or c^-1 factor is required to bound
the extra term A1 fz. It generally is nonzero. For any integer p>=1,
write n=m+2p+3, B_m=N_m(k_e chi_U), A_n=N_n(f_e). Then

    N_m(Ae fz) <= C_(m,p) M^-p N_n(fz),
    C_(m,p)=2^(m+n)L B_m A_n [8*3^p*2^(m+p+3)].         (8)

For each nonzero mode, integrate the whole-line radial integral p times
by parts; smooth compact extension kills all boundary terms. Use p
radial derivatives, the divisor bound to power p, and a+p+3 torus
derivatives to leave the same summable r^-3. Total input order is at
most m+2p+3. Leibniz on the fixed multiplier gives the constant in (8).
The zero Fourier mode vanishes because bar(fz)=0. A nonzero zero-mode
radial moment instead leaves an M-independent cutoff term; it is not flat.

The retained floor gives M^-p<(Tg S)^(rho_g p)epsilon^(p kappa).
Replace its S exponent by ceil(rho_g p) for a polynomial bound. For any
fixed physical loss and requested order choose p with STRICT positive
exponent margin. This gives flatness on the same punctured domain with
finite constants, not a derivative-dependent new tail or a finite-order
certificate of all flatness orders. Finite regularity alone is insufficient.
All actual A0/A1 velocities and full pressure memory remain in subsequent
products and pointwise residuals, including when (8) is used.

Next establish the actual canonical seed and current sourced/fast residual
norms Etheta_n,Ez_n with their support/edge estimates, then verify the
complete viscous, pressure, cutoff and nonlinear residual. These explicit
operator estimates do not independently establish those inputs or blow-up.

## 5. Evidence

Exact integer quadratic norms test the actual divisor arithmetic. Separate
formal phase fixtures with rational nonzero frequencies verify compact
radial inversion, moving-cutoff divergence, (7), its sign and the nonzero
A1 term without approximating oscillatory integrals. Their finite endpoint
regularity and frequencies are not actual canonical fields. Independent
positive polynomial substitutions check (6) and the channel-specific jet
counts. Analytical arguments above, not passing fixtures, establish the
stated conditional continuum estimates.

```sh
python3 research/navier_stokes_cascade/results/core_fast_inverse_budget/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_fast_inverse_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_fast_inverse_budget/evidence.json \
  --output /tmp/core_fast_inverse_budget_fresh.json
```
