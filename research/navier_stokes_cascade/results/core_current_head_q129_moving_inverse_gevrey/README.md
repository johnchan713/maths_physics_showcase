# Q129: closed real angular bounds for the same moving phase inverse

Parent: `ef19ac81120289a8fff635938a6babbe45bcad68` (Q128).
Fresh fetch on 3 November 2026 agreed with this head. Q128 research CI
41681262022 and field CI 41681262055 were in progress; Q126 field CI
41659390704 had passed. Current notes, Q128 frozen evidence, Q64 and the
original loop bounds were read first. Historical evidence, accepted root,
phase origin, frequency, pressure datum, supports, scales and C_setup^-20
width remain unchanged. The setup activation result remains local.

Q128 closes the real cutoff/signed-root angular majorant conditional on
PREMOD complex caps. Q64 gives a finite moving-inverse procedure, but its
high slow jets were not yet a closed Gevrey bound. Here we prove a closed
real angular inverse/composition/zero-mean primitive estimate. The precise
mixed density and source inputs below still require binding to the same
actual loop from Q128. Neither slow analyticity nor a final repair bound
is inferred from this conditional estimate.

## 1. Original phase and the exact remaining mixed input

Write eta for the REAL slow variable and theta for the fast loop angle.
Use the original normalized phase

    phi(eta,theta)=integral_0^theta w(eta,s) ds,
    w=a(1+t_loop^2)/(2*pi*v), integral_0^(2*pi) w dtheta=1.

The same orientation-preserving inverse satisfies phi(eta,Theta)=zeta.
Its real density floor is w>=d>0 (the original d=G^-3 suffices).
The origin is fixed at theta=0; it cannot be reselected at each slow point.
For every finite requested order, require the uniform mixed input

    |partial_eta^i partial_theta^j w|/(i! j!)
                      <= W B^i C^j (i!)^2,          (1)
    W,B,C>=1.

This is a real Gevrey-3 angular bound with an analytic fast-phase cost.
For an all-order statement W,B,C are independent of i,j; a degree-m
application needs these bounds only through i+j<=m. It is NOT a complex
eta extension. Q64 supplies fast-phase analyticity;
Q128 supplies s,mu slow jets, and the original removable loop formula is

    t_loop=-b/a+d0*mu*integral_0^1 q_z(u*mu*p,theta) du.

Binding (1) requires quantitative mixed composition of this formula and
its positive denominators, with the SAME full p=p2 from Q128. In particular
incoming Pi and Pi_eta cannot be dropped. We do not declare that binding
proved by entering rational W,B,C into the helper.

At each real inverse center take theta0 in [0,2*pi]. Since 2*pi<8,
phi's normalized coefficients of positive total order satisfy

    |c_ij| <= A B^i C^j (i!)^2, A=8W.              (2)

For j=0 integrate the i-th eta derivative over the fixed interval from 0
to theta0. For j>=1 use c_ij=(partial_eta^i partial_theta^(j-1) w)/
(i! (j-1)! j), bounded by (2) because C>=1. The constant phi value is
irrelevant. Normalization makes the positive slow derivatives periodic
in the lifted angle; no periodic value bound for Theta itself is claimed.

## 2. A convergent majorant for finite REAL inverse jets

At fixed zeta put u=Theta(eta0+h,zeta)-theta0 and w0=w(eta0,theta0).
In the finite chain rule the only degree-m unknown is w0*u_m. All other
terms come from c_ij h^i u^j excluding (0,0),(0,1). Divide by w0>=d.
For any positive composition i+k1+...+kj=m,

    (i! product k_l!)^2 <= (m!)^2.                 (3)

Consequently the transformed coefficients |u_m|/(m!)^2 are dominated by
the coefficients of the NONNEGATIVE analytic majorant V determined by

    V=(A/d)[1/((1-Bh)(1-CV))-1-CV].               (4)

Here u_m=partial_eta^m Theta/m!. Equation (4) is a majorant only, not
an assertion that the actual slow Taylor series converges. Induction in
finite degree proves domination, since the excluded linear term is the
only one that could involve the unknown coefficient of that degree.

Set x=Bh, T=CV, k=AC/d, K=max(1,k), and

    t=min(1/8,1/(64K)), eps=min(1/8,t/(16K)).       (5)

On |x|<=eps, |T|<=t, the right side of
T=k[(1-x)^(-1)(1-T)^(-1)-1-T] has constant part <=t/14.
Its T derivative is bounded by 8k(eps+t)<=17/128<1/4.
Indeed, for real 0<=x,T<=1/8 the positive geometric series gives
(1-x)^(-1)(1-T)^(-2)-1<=8(x+T), and absolute values of the complex
series have the same bound. The image is within t/14+t/4<t/2.
Uniform contraction gives a holomorphic NONNEGATIVE majorant T on |x|<eps.
Picard iteration from zero identifies its positive formal coefficients.
Cauchy and (3) now give the closed estimate

    |partial_eta^m Theta| <= (t/C)(B/eps)^m(m!)^3,
    m>=1, at fixed zeta.                           (6)

The actual inverse is smooth by the real IFT and the original floor.
Only this auxiliary majorant is holomorphic in h. No division by pressure,
mu or the cutoff amplitude, and no change to the signed root, occurs.

## 3. Source composition and the SAME zero-mean primitives

For an original unpulled loop source f require the mixed real bound

    |partial_eta^i partial_theta^j f|/(i!j!)
                      <= F Bf^i Cf^j(i!)^2,
    F,Bf,Cf>=1.                                   (7)

Let R=eps/B and use the SAME majorant V in (4). On |h|<=r<=R/2,
V(r)<=2(t/C)r/R since V(0)=0 and its nonnegative coefficients have
Cauchy bound (t/C)R^-m. Choose

    r=min(R/2, R*C/(8*Cf*t), 1/(4*Bf)).            (8)

Then Bf*r<=1/4 and Cf*V(r)<=1/4. By the finite chain rule and (3),
the transformed jets of g(eta,zeta)=f(eta,Theta(eta,zeta)) are bounded by
F/[(1-Bf*h)(1-Cf*V(h))], whose value at r is <=16F/9<2F. Thus

    |partial_eta^m g| <=2F r^-m(m!)^3, m>=0.       (9)

For the original primitive source the required mean identity is precisely

    integral_0^(2*pi) f(eta,theta)*w(eta,theta) dtheta=0. (10)

For the ACTUAL original loop sources, (10) follows algebraically from
its accepted variance and mean identities. Write t_s=-b/a,
v_s=a(1+t_s^2), rho=zeta(v_s)^2*(2+delta/2-v_s), and v=v_s+rho.
The removable q_z integral has zero theta mean, hence average t_loop=t_s;
the signed variance equation gives a*average(t_loop-t_s)^2=rho.
Thus a*average(1+t_loop^2)=v, proving the phase normalization. With

    aL=v/(1+t_loop^2), bL=-v*t_loop/(1+t_loop^2),
    fA=-(aL-a)/2, fB=E*(bL-b)/2,

one has integral aL*w=a and integral bL*w=-a*t_s=b, so both weighted
source means are exactly zero. All integrals here are over theta in
[0,2*pi]. Crucially the denominator uses the AUGMENTED v, not v_s.
This proves the loop source mean identities conditional on the original
accepted mean/variance equations. It is not the full wave-covariance mean
correction with axial velocity and pressure, which remains open.

Equation (10) is equivalent by change of variable to integral_0^1 g dzeta=0.
An unweighted theta mean is insufficient. Differentiation is at FIXED
zeta; it includes phi_eta and motion of Theta. Retain the original
zero-mean primitive

    P(eta,zeta)=integral_0^zeta g(eta,s)ds
                    -integral_0^1 integral_0^v g(eta,s)ds dv.

It is periodic by (10), with exactly zero mean. Differentiating these
fixed finite integrals and using (9) gives

    |partial_eta^m P| <=4F r^-m(m!)^3, m>=0.       (11)

No mean correction or integration constant has been changed. Apply (7),
(10) separately to Q64's Acal_zeta and Bcal_zeta sources. At the original
fixed integer N, the actual axial modulation Umod=U+Bcal(eta,Ny)/N obeys

    |partial_eta^m Umod| <= Uraw_m+4F r^-m(m!)^3/N. (12)

This is before the LATER exact five-row repair. It is not the missing
final compact bound in Q127's global maximum. Pure angular differentiation
holds y fixed; (12) does not give radial smallness or alter Q64's N^(r-1)
radial cost. Exponentiated swirl still needs its own composition bound.

For (6),(11) the binary bound is ceilLog max(1,prefactor)+
m ceilLog max(1,base)+3m ceilLog m (last term zero at m=0).
The inverse lift rejects m=0. Order 100005 requires no factorial expansion.

## 4. Checks, limitations and next obligation

Exact rational tests check both contraction gates and the composition disk
for small and enormous inputs. Independent implicit polynomial roots check
(6); positive finite convolutions check (3),(4) and detect omission of the
Gevrey factorial. A periodic phase control detects the omitted moving-origin
term. A weighted-mean countercontrol shows why zero theta mean cannot be
used as (10); periodicity and the fixed zero-mean primitive are checked
exactly. Moderate exact bounds are enclosed by the high-order binary budget.
These are manufactured algebra diagnostics, not evaluation of the actual
selected parameters or continuum proof by testing. The audit replays Q128
and checks every source/evidence hash in its frozen dependency closure; fresh evidence
must agree byte-for-byte. Historical frozen records remain untouched.

Next bind (1),(7) and the accepted mean/variance premises to the actual removable loop and full incoming
pressure state, using Q128, then propagate through exponentiated swirl and
the SAME later five-row repair. Actual PREMOD state caps/floors, the final
compact norm and cone join, coefficient/reserved-mean supports, uniform
physical pulse sizes/derivatives/supports and complete wave-covariance
mean updates remain open. The viscosity/cutoff/all-interaction residual is
unproved. No unconditional pulse, full residual improvement or blow-up
conclusion follows.
