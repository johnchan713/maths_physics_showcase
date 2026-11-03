# Q128: closed real Gevrey-3 bounds through the variance cutoff and signed root

Parent: `5a0fcb117a907fa6ead1d802c27e604e478f640e` (Q127).
On 3 November 2026 the freshly fetched branch agreed with that parent.
Q127 research CI 41671058791 and field CI 41671058659 were in progress;
Q125 field CI 41649355674 and Q124 field CI 41641232330 had passed.
Current notes, frozen evidence and Q64--Q67 were inspected first.
Historical evidence, cutoff, accepted roots, frequency, supports, pressure
datum, scales and C_setup^-20 width are unchanged. The setup result remains
local leading core/first activation collar.

Q127 supplies complex PREMOD estimates, but the variance cutoff need not
preserve slow analyticity. Here is a closed REAL angular Gevrey-3 estimate
through that cutoff and the SAME signed variance inverse, using finite
complex value/first-radial-derivative caps and genuine incoming memory.
Q66's finite jet construction remains valid; the new bound avoids
constructing its order-100005 polynomials. This checkpoint does not yet
bound the moving phase, modulation primitives or final compact field.

## 1. Complete forward-state value caps before modulation

Fix an actual PREMOD interval [x0,x1], x0>0, with endpoints independent of
eta. On every complex eta disk of radius 0<rho<=1/2 about [-1,1], supply
uniform caps e,u,e_X,u_X,i_E for the SAME E,U,E_X,U_X,E^-1, together with
caps t_M,t_I,t_J,t_S,t_Pi for their five actual incoming primitive functions
at x0. These are complex analytic caps, not merely real value bounds or
sampled points. Q127 and the existing regional formulas provide the route
to these low radial orders; their identity/domain binding is still required.
In particular no complex Gamma heat profile is assumed.

Put ell=x1-x0 and r=max(1,2x1)>=sqrt(2x1). The five density caps are

    g_M=u, g_I=r e, g_J=r e u,
    g_S=u^2+e^2/2, g_Pi=e^2/(2x0).

Set Mbar=t_M+ell g_M, and similarly for all five primitives. Finite
integration with fixed endpoints preserves holomorphy. On the smaller
disk rho/2 their first eta derivatives are bounded by 2/rho times their
caps. Keep Aphys=1/2+h, D=1/2-h, L=1-2h eta^2, d=1-eta^2 and the
original Q67 identities

    w=X-2D eta M-d M_eta,
    a=1-2X E_X/E, b=2X U_X/E,
    p=p2=[-w U+D(M-eta M_eta)+4h eta S-d S_eta
                      +X(4Aphys eta Pi-d Pi_eta)]/(L E). (1)

Using |eta|<=2, |d|<=5 and |L^-1|<=(1-8h)^-1 for 0<h<1/100 gives
uniform analytic state caps on rho0=rho/2:

    A=1+2x1 e_X i_E, B=2x1 u_X i_E,
    W=x1+(4D+10/rho)Mbar,
    P=i_E/(1-8h) [W u+D(1+4/rho)Mbar
           +(8h+10/rho)Sbar+x1(8Aphys+10/rho)Pibar]. (2)

Thus |a|<=A, |b|<=B, |p|<=P. I,J are retained in the five-row ledger
although these three loop inputs use M,S,Pi. Both pressure terms and all
incoming traces remain. A real pressure value cap alone is insufficient:
replacing incoming Pi by Pi+c(eta) changes p by
X(4Aphys eta c-d c_eta)/(LE), equal to -X/E at eta=0 for c=eta.
Equation (2) charges the derivative through the complex trace bound;
it never sets that trace to zero on a velocity plateau.

## 2. Zero-free real radicand neighborhoods, only where needed

The retained real positive floor a>=a_*>0 is an additional original
loop hypothesis, not a consequence of the upper cap A. Define

    ra=min(rho0/2,rho0 a_*/(4A)),
    V=A+2B^2/a_*.

On ra, |a-a0|<=a_*/2, so |a^-1|<=2/a_* and the analytic function
v=a+b^2/a has modulus <=V. Keep the same variance window delta=delta_L,
Hstar=2+delta/2, g=delta/4. Where the REAL cutoff is nonzero,
H=Hstar-v>=g and H<=Hstar because real v>0. Put

    rh=min(ra/2,ra g/(4V)), C=2 max(1,Hstar/a_*).     (3)

On rh the change of H is <=g/2. The positive real branch of sqrt(H/a)
therefore extends holomorphically with modulus at most
2 sqrt(Hstar/a_*)<=C. This is a local estimating disk about real centers
in the closed nonzero-cutoff region. It is NOT a definition of sqrt(H)
where H may be negative. Outside the cutoff support the original amplitude
is zero; flatness gives the same zero jets at its boundary. These facts
are enough for the uniform REAL derivative estimate below.

## 3. Closed cutoff-amplitude Gevrey bound

Use the SAME Q66 variance shape K built from
phi(t)=exp(-1/t-1/(1-t)), normalized by its integral, with flat plateaus.
Its global raw derivative bound is

    ||K^(j)||<=729*8^(j-1)*((j-1)!)^3, j>=1,
    ||K||<=1.                                      (4)

This is not the squared-flat Borel cutoff, and no other smooth cutoff is
substituted. For a convenient common majorant, (4) is <=729*8^j*(j!)^3.
The amplitude is s=K(8(v-2)/delta-1) sqrt(H/a). The positive-order
normalized Taylor coefficients of analytic v are bounded in magnitude by V ra^-k. For a
single-variable analytic increment bounded by G R^k, its jth power has
order-m coefficient at most G^j R^m binom(m-1,j-1). Therefore finite
chain-rule composition and (j!)^2<=(m!)^2 give

    |partial_eta^m K(8(v-2)/delta-1)|
       <=729 Lz^m (m!)^3,
    Lz=(1+64V/delta)/ra.                            (5)

Indeed sum_(j=1)^m (64V/delta)^j binom(m-1,j-1)
=(64V/delta)(1+64V/delta)^(m-1). No analytic extension of K is used.
Multiplying its normalized Gevrey jets by the analytic radicand-factor
coefficients C rh^-k, and using
sum_(i=0)^m x^i y^(m-i)<=(x+y)^m, gives

    |s^(m)| <= Ms Ls^m (m!)^3,
    Ms=729C, Ls=Lz+rh^-1, m>=0.                     (6)

The factorials in (5)--(6) are real derivative costs. They are not evidence
of a convergent slow Taylor series. The derivation applies at both flat
cutoff boundaries by smooth extension and does not divide by s or H=0.
The only radicand centers used have the original positive gap g.

## 4. Same signed root, including zero pressure and variance

Use Q65's scalar input-space inverse mu=Psi(p,s). From the unchanged real
center caps |mu|, |p|, |mu*p|<=G with G>=2, positive 0<d0<=1,
and an outward E_exp>=exp(G), its explicit radii
give a_Psi,sigma_Psi>0 and, for i+j>=1,

    |partial_p^i partial_s^j Psi|/(i!j!)
                 <=a_Psi(2/sigma_Psi)^(i+j).        (7)

The accepted real value bound is |mu|<=G. The input-space inverse is signed
and removable at zero variance/pressure; no inverse of p, mu or s occurs.
The condition E_exp>=exp(G), and all original real value caps/floors, are
premises. Computing the rational radius formulas does not prove them for
the selected profile or evaluate exp of its enormous G.

Since Ls>=rho0^-1, the positive normalized eta jets of p and s together
are bounded by N Ls^k(k!)^2, N=P+Ms. In a product of j such increments,
all positive index compositions k1+...+kj=m satisfy
product (ki!)^2<=(m!)^2, and their count is binom(m-1,j-1).
Use Q65's finite multivariate chain-rule bound (7); summing over j gives

    |mu^(m)| <=a_Psi Lmu^m (m!)^3, m>=1,
    Lmu=Ls[1+2(P+Ms)/sigma_Psi].                    (8)

At m=0 retain the original G, not a_Psi. This uses analytic scalar
input-space coefficients at the accepted real centers, not convergence of
a composed slow series through a flat cutoff. It identifies the same
signed root as Q65; no new root, frequency, band or cutoff is selected.

Equations (6),(8) are the new closed all-order REAL angular estimates from
finite low-radial PREMOD complex caps and original value floors. They do
not assert slow analyticity of s, mu or the final modulated field. For raw
jets M L^m(m!)^3 an integer binary upper budget is

    ceil(log2 max(1,M))+m ceil(log2 max(1,L))
                         +3m ceil(log2 m),          (9)

with the last term zero at m=0. Thus order 100005 can be budgeted without
factorials, high powers or order-100005 polynomial expansions. Normalized
jets cost division by m!, leaving (m!)^2; that conversion cannot be omitted.
No retained-scale inequality is inferred merely from these finite bounds.

## 5. Checks, binding and remaining obligation

Independent forward-state polynomial-memory controls retain both pressure
terms, all incoming traces and their angular derivatives. Real plateau
radicand derivatives are checked against (6). Positive finite convolution
controls check both Gevrey composition inequalities and detect dropping
the normalized factorials. Exact moderate derivative bounds are enclosed
by the binary budgets, including m=0 and zero pressure-cap cases. These
are manufactured algebra controls, not a selected-parameter evaluation or
continuum proof by sampling. The audit replays Q127's immutable evidence
and source closure; fresh output must match byte-for-byte. Resource-heavy
replays are run separately after the concurrent kills observed in Q127.

Next supply and bind the finite complex caps in (2) on each actual active
PREMOD region, using Q127 and the original low-radial regional formulas.
Then propagate (6),(8) through the MOVING phase inverse and its zero-mean
primitives, with the original density floor, and through fixed-frequency
modulation and the SAME later five-row repair. The final compact leading
norm remains open. Full pressure memory and the complete nonlinear wave
covariance/axial/pressure mean update, original scale and high-pole source
contracts, uniform actual physical pulses and the full viscous/cutoff/all-
interaction residual remain unproved. No unconditional pulse bound,
complete residual improvement or blow-up follows.
