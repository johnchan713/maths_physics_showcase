# Q124: all-order full-pressure moment and incoming repair discrepancy bounds

Parent: `921cfa63cb051c6be7cefe3f0fb42247c15451f6` (Q123).
On 3 November 2026 the freshly fetched branch agreed with this parent.
Q123 research CI 41631558113 and field CI 41631558011 were in progress;
Q122 field CI 41622757388 and Q121 field CI 41614701499 had passed.
Current notes, frozen evidence, background_summation, the coefficient
support foundation and higher_pressure_source_moments were read first.
Historical records, scales, supports, roots and pressure datum are preserved.
The core_growth_activation scope remains local leading core/first activation
collar and its width C_setup^-20 remains unchanged.

The full ordered pressure moment IDENTITIES already exist at every order.
Their explicit quantitative estimate in higher_pressure_source_moments was
only C_eta^1 and its repair discrepancy substitution was order two. Here we
extend those estimates to every fixed angular order and every positive
coefficient n. This supplies all five incoming discrepancy norm formulas
needed by Q122/Q123, conditional on actual whole-support lower-field and
cut-seed norms. It does not supply those native norms or a pointwise source
bound. In particular signed total control is NOT absolute-integral control.

## 1. Fixed angular algebra and a closed reciprocal estimate

Write |f|_m=sum_(j<=m) sup|partial_eta^j f|/j! on [-1,1]. This is a
product algebra; |fg|_m<=|f|_m|g|_m and
|f_eta|_m<=(m+1)|f|_(m+1). These norms are not raw derivative maxima.
Keep A=1/2+h, D=1/2-h, L=1-2h eta^2, d=1-eta^2, 0<h<1/100.
A sufficient all-order reciprocal bound is

    |1/L|_m <= R_m:=4^(m+1).                          (1)

Indeed the positive normalized Taylor majorant at every center is
[ell-4ht-2ht^2]^-1, ell=1-2h. At t=1/4 the subtracted variation is
9h/8<ell/2, so the series is at most 2/ell. Every normalized coefficient
is <=(2/ell)4^j. Summing through m gives less than
[2/(3ell)]4^(m+1)<4^(m+1). This uses analyticity of L's reciprocal only,
not analyticity of the profiles. The earlier sharper low-order bounds
remain unchanged.

For T_b f=(d f_eta+2b eta f)/L, |d|_m<=4 and |eta|_m<=2 imply

    |T_b f|_m <=4R_m(m+1+|b|)|f|_(m+1),
    |T_(b-D)T_b f|_m
      <=16R_m R_(m+1)(m+1+|b-D|)(m+2+|b|)|f|_(m+2). (2)

Both angular operators in the shifted axial viscosity remain. No division
by d is used, including at eta=+-1.

## 2. Full signed pressure moments at every fixed m

Let lambda_i=2hi, M_i=int_0^X U_i, v_i=V_i/X with its regular axis value.
The regular average M_i/X=int_0^1 U_i(tX,eta)dt and the divergence relation
from the inherited hierarchy give, for every r>=0,

    sup_X |v_i|_r
     <=R_r[4+4(D+lambda_i)+4(r+1)]sup_X|U_i|_(r+1). (3)

All U bounds here hold on the ENTIRE [0,Xv], containing every finalized
lower support. The exact zero total primitive and smooth endpoint joins
are hypotheses needed for the full-support integrated identities. Local
vanishing on I3 or I4 cannot replace these hypotheses or reset pressure.

For pressure order n>=1 put k=n-1. Supply u_i>=sup_X|U_i|_(m+2),
i=0,...,k. Define

    v_i^*=4R_(m+1)(m+3+D+lambda_i)u_i,
    B=sum_(i+j=k) u_i v_j^*, C=sum_(i+j=k) v_i^* v_j^*.

These bound the required |v_i|_(m+1), the ordered U_i v_j products and
the ordered v_i v_j products. Every orientation is retained. For k>=1
supply in addition u_prev>=sup_X|U_(k-1)|_(m+3) and put

    p=4R_(m+2)(m+4+D+lambda_(k-1))u_prev,
    beta_j=lambda_(k-1)+j,
    H_j=16 nu R_m R_(m+1)
          (m+1+|beta_j-D|)(m+2+|beta_j|)
          Xv^(j+1)/(j+1) *p, j=0,1.                 (4)

For k=0 set H_j=0; no fictitious predecessor is requested. By the existing
exact signed moment identities and (1)--(4), I_kj=int_0^Xv X^j Omega_k/X
obeys the following FACTORIAL angular algebra bounds:

    Q0=R_m Xv[(lambda_k+2D(m+1))v_k^*
                +(4(m+1)+4|lambda_k-A|)B]
          +Xv C/2 +2nu v_k^* +H_0,
    Q1=R_m Xv^2/2[(1+lambda_k+2D(m+1))v_k^*
                +(4(m+1)+4|lambda_k+D|)B]
          +Xv^2 C/4 +H_1.                           (5)

The unweighted axis contribution 2nu v_k(0) is bounded, not silently set
to zero. Actual positive-order zero axis traces can remove it only after
that identity is established. Weighted radial viscosity cancels, but H_1
usually does not. Shifted axial viscosity first appears at n=2 and needs
U_0 through m+3 there. At m=1 this recovers the inherited order demand:
U_i through 3 and the shifted predecessor through 4. The constants here
are conservative and need not match the sharper old C_eta^1 constants.
No radial derivative norm enters these TWO signed totals; this does not
remove radial derivatives from the local pointwise source.

## 3. All five ordinary incoming discrepancies at general n

Let 0<Xcut<=Xvel<=Xv cover the same tentatively cut current coefficient
and the true positive lower-field supports. Supply

    Jf>=sup_[0,Xcut] |F_n^cut|_m,
    Ju>=sup_[0,Xcut] |U_n^cut|_m,
    f_i>=sup_[0,Xvel] |F_i|_m, i<n.

The u_i above also dominate their lower-order angular norms. Set

    FF=sum_(i=1)^(n-1) f_i f_(n-i),
    UF=sum_(i=1)^(n-1) u_i f_(n-i),
    UU=sum_(i=1)^(n-1) u_i u_(n-i).

An empty sum is zero. In the ORIGINAL ordinary moment order, the complete
incoming discrepancy norms satisfy

    mu1=Xcut Ju,
    mu2=Xcut^2 Jf,
    mu3=2Xcut f0 Jf +Xvel FF +Q0/2,
    mu4=Xcut^2(u0 Jf+f0 Ju)+Xvel^2 UF,
    mu5=2Xcut u0 Ju+Xcut^2 f0 Jf
          +Xvel UU+Xvel^2 FF/2+Q1/2.                (6)

To derive (6), retain Pi_n,X=2F0 F_n+sum_(1<=i<n)F_i F_(n-i)
-Omega_(n-1)/(2X), E_i=sqrt(2X)F_i and the five ordinary total densities
from coefficient_support_foundation (2). Integrating 2X over the current
cut support yields Xcut^2; integrating X over the lower support yields
Xvel^2/2. All lower-order squares and cross terms remain fixed incoming
terms, not a nonlinear current-order repair. At n=2, with Jf=Ju, (6)
reduces exactly to the previously derived second-repair discrepancy formula.
The pressure terms Q0,Q1 are through Xv even when Xcut<Xvel<Xv.

Q122's discrepancy inputs use RAW maxima through m, so a sufficient
conversion is

    D_i=m!*mu_i, i=1,...,5.                           (7)

This factorial must not be discarded. The moderate-order API exposes (7)
with exactly the Q122 keys m1,...,m5. For a tuple request of order r,
Q122/Q123 need discrepancies through m=r+1; then (3) requires ordinary
lower U through r+3 and the shifted predecessor through r+4. Actual
angular input availability must be checked at those increased orders.
This is a bound on the complete incoming moment FUNCTIONS, not merely
pointwise values or pressure derivatives on a reserved mean interval.

## 4. Deferred binary implementation and checks

bounds(binary=True) accepts integer binary upper budgets for the SAME
native product-algebra inputs. Positive constants are padded to one before
ceil(log2); a zero term has a separate zero sentinel. Product budgets add,
and a sum of l nonzero terms costs max budget+ceil(log2 l). The geometric
R_m budget is exactly 2m+2, without evaluating 4^(m+1). Equation (7)
costs at most m ceil(log2 m), zero for m=0. Thus all five raw discrepancy
budgets at m=100005 can be calculated without factorials or high powers.
These outputs still require actual native input proofs. They do not select
any a_n or certify Q120's retained-scale inequality. The exact-number API
is deliberately capped at moderate angular orders.

Independent symbolic controls form the regular full order-two pressure
source directly, including both ordered lower interactions and the nested
axial operator. They integrate it and all five ordinary incoming densities
in X before checking angular jets at endpoints and interior eta. Polynomial
coefficient norms supply conservative whole-rectangle input bounds. These
are manufactured algebra checks, not actual coefficient solutions or a
continuum proof by sampling. Further checks recover the old n=2 formula,
compare binary bounds with moderate exact values, and distinguish zero
previous input from an omitted shifted input. Historical full-support and
finite-endpoint identity tests remain part of the inherited suite.
The audit replays Q123's complete source/evidence closure; a fresh record
must match byte-for-byte.

Next prove actual cut-seed and globally finalized lower F/U input bounds
on these domains at the stated orders, retaining the leading/core and
full-support identities. Then insert the resulting five raw budgets into
the current repair/tuple/schedule chain and verify its original scale gates.
The native order-100005 pointwise source/common-pole obligations and the
compact-modulation/exact-repair cone join are not closed by signed totals.
Global support/mean identities, complete-wave covariance means including
axial velocity and pressure, uniform actual physical pulses and the full
viscous/cutoff/interaction residual remain open. No unconditional pulse
bound, complete residual improvement or blow-up follows.
