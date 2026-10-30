# Q32: physical stress normalization and the native root input

Parent: `a951d7094b44df5dacfcd04f30ee4b46f4aaa9f6` (Q31).
Q31 deliberately stopped before identifying its scalar Tsum with physical
stress. This checkpoint supplies that normalization and its finite-order
derivative costs, including the optional U=R^d chart used by compact maps.
It derives the epsilon-cubed higher-stress input used by the initial Haar
row estimate. All statements remain relative to Q31's exact-profile and
combined Borel premises and the retained chart/support hypotheses.

## 1. Distinguish the three stress quantities

Let A=1/2+h, alpha=A+1/2=1+h, Q=2^-ell, epsilon=Q^h and s=q/Q.
The physical stress and original leading box target, as specified in
`core_localized_stress_assembly` section 2, are

    Sigma_phys=q^-alpha Tsum,
    T_box=s^-alpha T0.

The normalized velocity is Q^A times physical velocity. Its covariance,
and hence the stress in its normalized equation, scale by Q^(2A), giving

    Sigma_norm=Q^(2A) Sigma_phys=epsilon s^-alpha Tsum,
    Sigma_norm,0=epsilon T_box.                         (1)

Indeed 2A-alpha=h. In particular the canonical equation is H y=T_box,
with both signs inside H, and the normalized principal flux is epsilon T_box.
Q31's T0/zeta bound must therefore be multiplied by s^-alpha before it
is a bound for this physical chart's root target. This clarifies notation
where older root recipes call the box target simply T; it does not change
any root or insert a new partition factor.

The complete correction is

    Sigma_norm-Sigma_norm,0=epsilon s^-alpha (Tsum-T0). (2)

Q31 bounds the parenthesized expression by q^(2h) zeta times a finite
derivative constant and delta^-3m-3. Since q^(2h)=epsilon^2 s^(2h),
(2) has normalized epsilon order THREE. It is not order two, and it is
not an absolute physical decay claim: removing Q^(2A) gives the value
power Q^(-1+h). No residual cancellation follows from (1)--(2).

## 2. Differentiate the full scale multiplier on the enlarged chart

Use the same closed chart s in [1/4,4], X=R^2/(2s), eta=Z s^(-(1/2-h)),
and fixed positive annular bounds Rmin,Rmax. Let Y_m and Comp_m be the
zero-constant logarithmic increment sum and finite Taylor composition
constant of `core_background_finite_head_budget` equation (7).
Define Exp_m(x)=sum_(k=0)^m x^k/k! and

    C_m=16 m! Exp_m(2Y_m) Comp_m.                       (3)

For the actual 0<h<1/100, alpha<2 and s_ref^-alpha<=16. Write
s^-alpha=s_ref^-alpha exp(-alpha log(s/s_ref)). Its normalized finite
jet sum is at most 16 Exp_m(2Y_m). Multiplying by the finite composition
jet and converting a degree at most m coefficient to a raw derivative
costs m!, proving (3). This is finite Taylor algebra with remainder,
not analyticity of the heat profile. Every derivative of s^-alpha is
retained. Setting that factor to its value would miss terms even for
constant T0.

With Q31's T_m,S_m and its pointwise delta at the base point, (3) gives

    |Dslow^I (T_box/zeta)| <= C_m T_m delta^-3m-3,
    |Dslow^I (Sigma_norm-Sigma_norm,0)|
        <=4 epsilon^3 C_m S_m zeta delta^-3m-3, |I|<=m. (4)

For the second line all profile/Euler derivatives in Q31 already
differentiate q^(2h). The only extra conversion is q^(2h)<=4epsilon^2;
there is no second independent differentiation of that extracted bound.
For the first line zeta is inside the original quotient before any
composition. No inverse edge floor or extra pole is introduced. The
same bound with T_m also controls the leading weighted stress itself,
because Q31's T_m dominates its original leading weighted raw bound A_m.

## 3. The optional radial change and factorial-normalized root input

The stress is a scalar component under the auxiliary reparameterization
U=R^d (d>0 fixed). This is composition of the cylindrical stress entries,
not a change of the physical vector basis and not an integration measure.
For U>=Umin>0 and R<=Rmax, define

    r_j=|(1/d)(1/d-1)...(1/d-j+1)| Rmax Umin^-j/j!,
    Z_m=sum_(j=1)^m r_j+4,  m>=1; Z_0=0,
    J_m=m! sum_(k=0)^m Z_m^k/k!.                       (5)

The r_j bound the normalized derivatives of R(U). The four identity
coordinates pad the rest of the five-variable coefficient chart; unused
auxiliary coordinates only enlarge this bound. For a raw through-m
bound M, multivariate finite Taylor composition has coefficient sum at
most M sum_(k<=m)Z_m^k/k!, by the multinomial theorem. This proves the
cost J_m in (5). There are no additional input derivatives, epsilon powers
or edge poles. For d=1 the padded cost is still valid, though not sharp.
For the R chart this conversion is omitted. No uniformity as Umin tends
to zero is claimed; edge/axis charts still have their separate obligations.

Thus put Cnative_m=C_m in the R chart or J_m C_m in this U chart. The
original root budget's total-degree Taylor norm in five variables is bounded by

    target_m = (5^m/m!) Cnative_m T_m,                 (6)

with the same pole 3m+3. This uses sum_(|I|=m)1/I!=5^m/m!, not an
identification of raw derivatives with Taylor coefficients. `transfer`
returns (6) as a positive constant polynomial compatible with the original
`RootBudget`, while preserving its independent column, wave, complete-curl,
loghalf and value-floor inputs. It also returns 4 Cnative_m S_m as the
raw higher-stress constant AFTER extracting epsilon^3 and the full weight.
Labels and Q are frozen; none is differentiated.

## 4. Join and remaining obligations

The exact physical leading assembly remains

    sum_beta eta_beta^2 Q_beta^(-2A) epsilon_beta T_box,beta
      = (sum_beta eta_beta^2) q^-alpha T0 = Sigma_phys,0.

This identity includes derivatives of all partitions and all scale factors
before cancellation. One must not count the two signs again or divide by
a partition. Complete curl squares, other tensor entries, axial velocity
and pressure are outside this principal identity and remain in the full
mean/source recipes. Higher-stress bounds (4) feed their M_3 input on each
accepted chart; existing overlap/common-atlas derivative costs must still
be applied when assembling norms. A bound on this scalar component does
not justify dropping those costs or propagating a root through arbitrary
new charts.

Established here: the missing physical prefactor, its all-order finite
majorant, the optional radial composition and the actual root-input norm
conversion. Next assemble these with the canonical entrance/column/curl
budgets and common-atlas estimates, and check every remaining lower
inequality on Dv before claiming the complete first source bound.
The background's full Pi1/Pi2 and later pressure memory is unchanged.
No complete residual improvement, all-scale pulse theorem or blow-up follows.

Tests use exact rational manufactured parameters, independent differentiated
coordinate fixtures, the original RootBudget interface and negative controls
for a missing scale factor or covariance square. They do not evaluate
continuum seminorms at the enormous fixed root parameter. Q31's evidence
and its transitive source hashes are replayed unchanged. On 30 October 2026
parent focused CI 41086336398 had passed; broad CI 41086336169 was running.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json` in this directory.
