# Q72: explicit mixed bounds on the short native and taper interval

Let t=C^(-20), Xa=4/Lambda, y=log(X/Xa), Y=4 exp(y), B=Xa exp(2t).
The interval [Xa,B] ends precisely where Q71 starts. We quantify the
finite-seminorm construction of [the prescribed continuation](../prescribed_continuation_foundation/README.md)
on this interval using the SAME infinite core and steps. The input constants
below involve frozen parameters and finite sums, not unspecified core
seminorms. No later matching root, pressure reset, or radial analyticity
of the cutoff is needed.

## Seed slopes in log radius and core traces

Use the polydiscs of the [all-order core](../core_all_order_foundation/README.md):
a=1/(100S), b=rho/(100S), Re Phi>1/5 and |Phi|<=2S there, since
the series majorant has margin >=5/8. Write D_s for the complex
s-neighborhood of the real eta interval. For real 0<=y<=2t, 4 exp(y)<4.1.
If |v|<=h0=a/16 then |Y exp(v)-Y|<=2Y h0<a. Thus on the product
of the complex v disk and D_b the principal log Phi has modulus at most
L0=log(10S)+2 (bound its real and imaginary parts separately).
U_nat has modulus <=U0=2(9+S/Lambda) by the same series majorant.
Cauchy bounds the normalized y derivatives of these functions of order
i by L0 h0^(-i) and U0 h0^(-i). Only the native seed uses these y disks.

For each requested radial order m use the norm sum_(i=0)^m
sup |partial_y^i v|/i! on the real interval and the indicated eta strip.
Products cost multiplication, exp costs exp of the norm. The primitive
from zero with trace v0 costs |v0|+(2t+1) times the integrand norm
(the additive factor multiplies only the integrand). A derivative in eta
from D_s to D_(s') costs 1/(s-s') with no radial order loss.
For psi=1-sigma((y-t)/t), kappa=1-(1-t)sigma(y/t), their norms cost
at most T_m=1+sum_(i=1)^m 162·3260^i(i!)²/t^i, by Q63 and the
real value bounds. In particular no complex y shift acts on a step.

At Xa the native F=g Phi and its inverse cost F0=2S and
V0=5C exp(2 Lambda Mcore) respectively, using the original core g and
the constant Mcore of the zero-free estimate. The same F0,U0 bound
the core down to X=0 on D_b. Hence the five trace bounds at Xa are

    m0=Xa U0,  i0=Xa² F0,  j0=Xa² F0 U0,
    s0=Xa U0²+Xa² F0²/2,  p0=20P²+Xa F0².                  (1)

These follow by integrating the original X densities with E=sqrt(2X)F.
The pressure datum costs <=20P² since the measure has mass <5 and
|1+eta²|>=3/4 on this tiny strip, with exponent in [0,2].
Thus the last row bounds the full Pi, not just its increment. The traces
in the two maps below are the same functions at Xa only.

## Reference fields and one explicit positive moment-map calculation

Put D_m(H)=H sum_(i=0)^m (i+1)h0^(-i-1), the norm bound for a first
native slope. Integrate the tapered slopes from the actual core trace:

    J_r=int psi (log Phi)_y dy,   F_r=F(Xa) exp(J_r),
    U_r=U(Xa)+int psi (U_nat)_y dy.

Here g is independent of y. On D_b put J=(2t+1) T_m D_m(L0).
Then F_r,1/F_r,U_r have norm bounds f=F0 exp(J), v=V0 exp(J),
u=U0+(2t+1)T_m D_m(U0). These are size bounds at order m, not
small perturbations. They include the first interval with psi=1.

Here is a reusable upper bound for the original state formulas, with
inputs f,v,u on D_s and outputs on D_(s/2). Let x=e B and
q=e/Xa bound the norms of X and 1/X (sum of normalized derivatives
of exp(±y)). Denote k=2/s, so each angular derivative costs at most k.
Integrating the five y densities gives, in the order M,I,Jmoment,S,Pi,

    (m1,i1,j1,s1,p1mem) = (m0,i0,j0,s0,p0) + (2t+1) times
      (x u, 2x² f, 2x² f u, x(u²+x f²), x f²).              (2)

All entries in (2) are magnitudes; signs are retained in the defining
fields. With |eta|<=2, |1-eta²|<=5 and |1/L|<=2 on the strips, set
A=1/2+h, D=1/2-h. Direct product bounds for the full forward map give

    w=1+q(4D+5k)m1,
    nQ=(1-h+2D k)i1+(5k+4|h-D|)j1,
    Q1=2[x w+q v nQ/2],
    N=2[w u+q(D(1+2k)m1+(8h+5k)s1)+(8A+5k)p1mem].         (3)

Indeed w bounds W, Q1 bounds p1=(-XW+N_Q/(2XF))/L, and N bounds
ns=N_s/L, whose last terms are (4A eta Pi-(1-eta²)Pi_eta)/L.
In particular the two reference source rows are Q1,N, not p1,p2.
This estimate uses all incoming row constants and does not assume a
three-power radial representation while F_r and U_r still vary.

## Actual slopes and independent memory

Use (2)--(3) first for the reference on s=b. On D_(b/2), before X=100
the actual equations are (log F_a)_y=-kappa p1r/2 and
(U_a)_y=-kappa X nsr/2. Here B<100. Set
J_a=(2t+1)T_m Q1/2. For the actual F,1/F,U one may therefore use

    fa=F0 exp(J_a), va=V0 exp(J_a),
    ua=U0+(2t+1) T_m x N/2.                                 (4)

Now repeat (2)--(3) on s=b/2 with THESE fields and the original traces
(1). This step cannot copy reference moments after the prescribed slopes
change. It bounds the actual integrated coordinates on D_(b/4): p1
by the new Q1 and p2=X ns/E by the new N times sqrt(B/2) exp(1/2) va,
since X/E=sqrt(X/2)/F. Also E=sqrt(2X)F costs sqrt(2B) exp(1/2) fa.
For b_s=2(U_a)_y/E use the prescribed slope in (4) before integration,
and 1/E bounded by exp(1/2)va/sqrt(2Xa); a=kappa p1r directly.
These are the pre-modulation state and fields, not the finalized pulse.

For any such order-m output bound V on D_(b/4), real mixed physical
derivatives satisfy, for 0<=i<=m and any j>=0,

    |partial_X^i partial_eta^j output|/(i!j!)
      <= (4/b)^j Xa^(-i) V/i! sum_(r=0)^i [i,r] r!,          (5)

by angular Cauchy and the signed Stirling formula for log X. Constants
are explicit at fixed order and may be enormous. At B the steps are
flat and the defining fields/moments join Q71's continuation without
altering a trace. Together with Q69--Q71 this supplies a common mixed
bound (six angular velocity, five angular state, every finite radial
order) on the native pre-modulation interval from Xa to Xrep under the
combined inherited premises. Higher angular orders on the later
continuation and annulus, the normalized repair defect and actual
physical-scale inequalities still require estimates. No uniform physical
pulse, improved pressure tail, complete residual improvement or blow-up
conclusion follows.

Run test_audit.py and audit.py --verify-record evidence.json --output a
fresh path. The analytical inequalities are distinct from the exact
manufactured algebra checks and frozen provenance replay.
