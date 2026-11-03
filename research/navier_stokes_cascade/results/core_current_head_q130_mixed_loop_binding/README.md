# Q130: bind mixed loop bounds to the six original scalar inputs

Parent: `41736482dd840aed1e802dfad2ea70495536cc40` (Q129).
Fresh fetch on 3 November 2026 agreed with this head. Q129 research CI
41689746649 and field CI 41689746383 were in progress; Q127 field CI
41671058659 had passed. Current notes, Q129 frozen evidence, Q128, Q64
and the original loop bounds were inspected first. Historical evidence,
accepted roots, phase origin, frequency, supports, pressure datum, physical
scales and C_setup^-20 width remain unchanged. The setup activation result
remains local leading core/first collar.

Q129 still requested independent mixed density/source bounds. Here we
remove those independent high mixed inputs: they follow from Q128's six
scalar angular jets and the original real loop caps. This is conditional
on the SAME actual PREMOD caps, positive floors and real value bounds;
entering numbers into a helper does not establish those premises. No
analytic continuation in eta through the flat cutoff is asserted.

## 1. The six scalar jets, including the augmented variance

Use Q128's complete state caps A,B,P for a,b,p=p2 on eta disks rho0,
original E cap e, positive real a floor a_*, and its radius ra where
|1/a|<=2/a_*. Keep its amplitude and signed-root bounds Ms,Ls and
Mmu=max(G_mu,a_Psi), Lmu. Set

    L=max(1,rho0^-1,ra^-1,Ls,Lmu),
    Mts=2B/a_*, Mv=A+2B^2/a_*+A*Ms^2,
    N=A+Mts+P+Mmu+e+Mv, L0=3L.                  (1)

The SAME six scalar parameters are lambda=(a,t_s,p,mu,E,v), where

    t_s=-b/a, v=a+b^2/a+a*s^2=v_s+rho.

They satisfy for every requested finite angular order m

    sum_l |lambda_l^(m)|/m! <= N L0^m (m!)^2.  (2)

For a,b,p,E use Cauchy on rho0 and pad by the Gevrey weight. For 1/a
use ra and 2/a_*. For mu retain its original value cap at m=0, and its
Q128 positive jets at m>=1. Products of k Gevrey-3 functions with common
base L have prefactor product M_l and base kL: each convolution satisfies
product(i_l!)^2<=(m!)^2, and the coefficient of the product of geometric
series is bounded by (kL)^m. Apply this to b/a, b^2/a and a*s^2.
The sum in v gives Mv in (1). This is the AUGMENTED v, not v_s; the latter
would fail the phase normalization and the original weighted means.
No pressure term is dropped: p is exactly Q128's state with incoming
Pi and Pi_eta, not a local pressure surrogate.

The original real loop bounds must also supply a common G>=2 for
|a|,|t_s|,|p|,|mu|,|E|,|v|,|mu*p| and |t_loop|, with a,v>=1/G.
These are sharper real value premises, not deductions from (2)'s huge
upper constants. Supply E1>=exp(G) as in Q128/Q65. Define E2=9E1^2;
since e<3 this gives E2>=exp(2(G+1)) without evaluating a huge exponential.

## 2. A joint analytic neighborhood in scalar inputs and fast angle

At any REAL lambda0 and real theta0, use exactly the removable loop

    z=mu*p, q(z,theta)=exp(z sin(theta))/I0(z),
    t_loop=t_s+d0*mu*integral_0^1 q_z(u*z,theta) du,
    0<d0<=1 fixed.                              (3)

Only the scalar inputs and fast theta are continued complexly. Put

    r0=min(1/4,1/(8G),1/(64G E1)),
    T=2G+24G E1 E2,
    r=min(r0/2,r0/[112T(G+1)]).                  (4)

On the seven-variable polydisk of radius r0 (six scalars plus theta),
|Delta z|<=2G r0+r0^2<=3G r0<=3/(64E1). Uniformly for u in [0,1],
|I0(u*z)-I0(u*z0)|<=E1(exp(|Delta z|)-1)<=2E1|Delta z|<1/2,
using the real integral bound |I0^(j)(u*z0)|<=E1 and I0(u*z0)>=1.
Thus |I0^-1|<=2 and |I0'|<=2E1. Since |sin(theta)|<2 and |u*z|<=G+1,

    |q|<=2E2, |q_z|<=12E1 E2, |t_loop|<=T.      (5)

The u integral remains holomorphic by uniform bounds. This representation
works at p=0 and mu=0 and never divides by either of them. Also
|v^-1|<=2G on this polydisk, and each scalar value is <=2G.

On the half polydisk Cauchy gives each of the SEVEN first partials of t
at most 2T/r0. Integrating along the segment to a point of radius r gives

    |t-t0|<=14Tr/r0<=1/[8(G+1)].                 (6)

Here the original real |t0|<=G is essential; using only a crude complex
T as that real cap would change the radius cost. Equation (6) implies
|t^2-t0^2|<1/4, so |1/(1+t^2)|<=2 on radius r. The loop denominators
are therefore uniformly zero-free in this scalar/angle neighborhood.
This is an estimating polydisk, not a new physical band or slow eta domain.

For the ORIGINAL density and two unpulled sources

    w=a(1+t^2)/(2*pi*v), aL=v/(1+t^2), bL=-v*t/(1+t^2),
    fA=-(aL-a)/2, fB=E*(bL+a*t_s)/2,

uniform complex value caps on radius r are

    W=4G^2(1+T^2), F=8G^3(T+1),
    |w|<=W, |fA|,|fB|<=F.                      (7)

Indeed |aL|<=4G, |bL|<=4GT, |a*t_s|<=4G^2, so the second source
costs at most 4G^2T+4G^3<=F; the first costs at most 3G<=F.
The coarser W discards only the harmless factor (2*pi)^-1<1.

## 3. Finite real composition supplies Q129's mixed inputs

Cauchy in the six scalar variables and theta gives normalized coefficients
bounded by (W,F,F)*r^(-|alpha|-j). For any positive eta index composition
k1+...+kk=i, product(k_l!)^2<=(i!)^2. Sum scalar increments using (2):
multinomial coefficients in (sum Delta lambda)^k dominate the sum of all
scalar monomials of total degree k. The finite chain rule then yields

    |partial_eta^i partial_theta^j (w,fA,fB)|/(i!j!)
              <= (W,F,F)*Bmix^i*Cmix^j*(i!)^2,
    Bmix=L0*(1+N/r), Cmix=1/r.                  (8)

For i>0 the remaining sum is
sum_(k=1)^i (N/r)^k binom(i-1,k-1)<=(1+N/r)^i; i=0 uses (7).
This keeps the factorial in the slow index i alone. It does not charge
an artificial (i+j)! Gevrey factor or request extra slow derivatives for
high fast j. The argument uses only finite real jets of the scalars;
it does not compose a convergent eta Taylor series through the flat cutoff.
For a finite total degree m, all premises need only the corresponding
indices through m. An all-order assertion uses common constants independent
of i,j. No radial regularity is gained by relabeling these angular bounds.

Equation (8) removes Q129's independent W,B,C,F,Bf,Cf inputs. Feed W,Bmix,
Cmix and the real density floor d=1/(8G^2) into its inverse recipe; feed
F,Bmix,Cmix into its composition/primitive recipe. This floor follows from
a>=1/G, 0<v<=G and 2*pi<8: w>=1/(2*pi*G^2)>1/(8G^2).
The generic assumption G>=2 alone does NOT imply the older G^-3 floor;
that floor requires its retained stronger hypotheses. This derivation
avoids importing it without those hypotheses. No physical scale changes.
Its fixed-origin,
density-weighted mean proof still applies to the SAME augmented v and
accepted variance root. This binds the Q128 scalar bounds to Q129's closed
inverse and zero-mean primitive bounds, including zero pressure/variance.
It does not perform the exponentiated swirl update or later exact repair.

## 4. Checks, scope and remaining inequality

Exact rational controls check the I0 disk margin, real-center variation,
zero-free quadratic denominator and all value-cap inequalities. Independent
polynomial scalar jets check (1),(2), including the augmented variance.
Finite positive mixed convolutions check (8) and detect omission of the
slow factorials. A removable zero-pressure loop control checks the actual
w,fA,fB mixed derivatives against (8). Binding tests feed the resulting
constants into Q129, keep its stage/mean guards, and evaluate its binary
budget at order 100005 without materializing factorials. These are
manufactured controls, not evaluation of the selected enormous profile.
The audit replays Q129 and checks the source/evidence hashes in its frozen
closure; fresh output must match byte-for-byte. Tests are not a continuum proof.

Next verify the actual Q128 PREMOD complex caps, real floors, G caps and
outward exp(G) input across every original active region. Those are still
premises; (8) eliminates an additional mixed high-jet obligation rather
than certifying them by finiteness. Then propagate through exponentiated
swirl and the SAME later five-row repair with complete pressure memory.
The final compact norm/cone join, coefficient/reserved-mean supports,
uniform physical pulse bounds, full wave-covariance axial/pressure mean
corrections and the viscosity/cutoff/all-interaction residual remain open.
No unconditional pulse bound, full residual improvement or blow-up follows.
