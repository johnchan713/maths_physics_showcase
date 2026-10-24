# Conditional symbolic two-edge target and unchanged signed-root join

Parent: `4fd985d9936c16d06a39437340bf4a8c9cb1b0d3`.
On 24 October 2026 its logarithmic pullback CI step passed; the full
run remained in progress without failures. The previous projector
seminorm boundary remains OPEN. The recorded normalized/support
conditions do not alone establish the actual uniform family jets.
This checkpoint advances the other quantitative parts of the target
route while retaining that lower input explicitly. No projector,
physical width, support, first band, gauge or root is replaced.

## 1. Keep the physical logarithmic shell unchanged

The parent supplies the pulled-back projected-source coefficients
H_e,n at epsilon order alpha=2-3kappa, weight zeta and common pole
q=2p. The actual source has EXACT zero total moment from its original
projector normalization. Both endpoints are fixed in x=log(X/Xa),
with L=log(Xb/Xa). Let rho=Xb/Xa>1 and choose the rational enclosures

    Llo=2(rho-1)/(rho+1) <= L <= Lhi=rho-1.

The lower inequality is the existing log-gap estimate; the upper
follows by integrating 1/t<=1 from1 to rho. These are proof constants,
not a change of the physical interval. The code recovers the SAME rho
from the parent's retained geometry: with l=2Rminus, u=Rplus/2,
ra=(2+da)/(2-da), rb=(2+db)/(2-db), rho=ra rb (u/l)^2.

Set r=min(1,Llo/4), a=w^2, b=4, k=max(0,ceil((q-3)/2)), and

    D(c)=1/(2c)sum_(j=0)^k binom(k,j)j!(r^2/c)^j,
    C=D(a)exp(b/(Llo-r)^2)+D(b)exp(a/(Llo-r)^2)
                        +Lhi Eedge(q,w)exp((a+b)/r^2). (1)

The historical two-edge proof applies with this SMALLER proof collar.
On each edge use the correct directed primitive, with exact zero moment
for the backward form. On the middle interval use Lhi and the same
weight floor. L>=Llo bounds the opposite-edge denominators safely.
The sum dominates all three branches. C is a finite positive ExpSum;
the actual giant exponentials are never numerically evaluated.

For the normalized exact-degree jets of F_e=int_left^x g_e,

    F_n <= C H_n+H_(n-1), H_-1=0.                    (2)

Pure slow jets use the integral bound; jets with an x derivative use
the source one order lower, with normalized factor1/i<=1. Through-n
source bounds are safe overestimates for these exact-degree inputs.
No additional incoming source derivative is introduced. Zero extension
and support remain conditional on the original smooth compact source.

## 2. Preserve the full reciprocal weight and return to native U

Multiplication Xi_e=-R^-e F_e has fixed log-chart coefficient bound
Rmin^-e E_j(e/2), using the parent's E_j for exp(b[t+L_s(t)]).
Convolve it with (2) to obtain K_e,n at the SAME full zeta and pole q.

For j>=1 the normalized jets of log(zeta^-1) satisfy

    |partial_x^j log(zeta^-1)|/j!
        <=(w^2+4)(j+1)delta^(-j-2)
        <=(w^2+4)(j+1)delta^(-3j).

The exact exponential-jet recurrence therefore gives reciprocal
coefficients B_0=1,

    B_n=(1/n)sum_(j=1)^n j(w^2+4)(j+1) B_(n-j).

This bounds J_n(zeta^-1) by zeta^-1 B_n delta^-3n. No bound on the
zeroth value of log(zeta^-1) is assumed. Convolution cancels the full
zeta and bounds xi_e=Xi_e/(epsilon^alpha zeta) at pole q+3n.

For the ORIGINAL native U=R^d_r chart,

    x=(2/d_r)log(U)-log(s)-log(2Xa), 1/3<d_r<5/2.

Let Ufloor=min(1,Rmin^3), valid on the same positive radial hull.
The normalized U jets of (2/d_r)log U are <=6/(j Ufloor^j).
The sum of the three native-to-log coordinate increments is thus
bounded by coefficients

    A_j=2 delta_j1+6/(j Ufloor^j)+L_s,j, j>=1.

Finite Taylor composition gives

    V_e,n=sum_(k=0)^n Vx_e,k [t^n] A(t)^k,
    J_n(xi_e)<=V_e,n delta^(-3n-q).                   (3)

The pole is monotone, so padding lower-degree terms costs no new
edge power. Both nonlinear chart directions are now explicit. The
weight, width squared and outer constant4 are retained throughout.

## 3. Feed the ORIGINAL root and full curl operator

Use V_1,n+V_2,n as a sufficient two-component signed-target norm.
The code substitutes it into the existing positive recipe

    Z_n=(Cmass*S/2)sum_(i+j+k+l=n) E_i (K*V)_j R_k W_l.

It uses the same returned primary U-chart root, inverse, reciprocal,
wave and complete curl coefficients, copied into the symbolic ring.
It never recomputes a positive root at an updated signed target. The
historical mapping bounds the signed principal coefficient at order
alpha-1/2 and the complete curl remainder at order alpha-kappa, with
the same edge allowances d=q. The curl needs one more target jet.
The target route itself requests no derivative beyond n of its source;
the initial source's safe covariance requirement remains n+13.

The exact principal cross identity is linear in the signed target.
Its signed square and all later complete covariance entries remain;
there is no complete residual gain implied by this norm mapping.
Original root/value/cone, entrance, profile and support premises are
unchanged. The actual projector family's uniform seminorms are STILL
required: this is a conditional composition, not verification of the
actual signed wave from normalization alone.

## 4. Evidence and next obligation

Tests compare the two-edge constant to historical diagnostic branches,
independently differentiate the full reciprocal weight, check native
nonlinear chart jets and reproduce the ORIGINAL signed-root formula
and its curl rule. Marked manufactured target/covariance fixtures
keep finite tests moderate; they are not the giant continuum field.
All historical frozen records are preserved and rebound.

Next identify or bound the actual fixed projector family so that its
native seminorm input is justified, then carry the resulting real signed
wave through its complete covariance update and freshly recomputed
mean/pressure state before the first sourced inverse. Retain signed
squares, later complex cross terms, whole c_rho P, axial velocity,
viscosity, cutoffs and all interactions. No unconditional pulse bound,
complete residual improvement, infinite realization or blow-up is proved.
