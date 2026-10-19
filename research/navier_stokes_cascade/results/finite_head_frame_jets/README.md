# Finite-head third-jet Borel bound and radial spectral-frame derivatives

Parent: `79e54e2be7c0b250516e1e0cc2eb9462ffb3bf98`.
Status: `finite-head-third-jet-Borel-and-conditional-radial-spectral-frame-C2`.

The preceding collar estimate supplies an explicit inner first-coefficient
J3. It does NOT by itself give the summed background's third derivatives.
This checkpoint makes the missing finite-head contribution explicit, then
propagates a valid third-jet background bound into first/second radial jets
of N,K,c and the positive growth rate. Physical radial scale losses remain.

This is a quantitative refinement of the ORIGINAL finite-head/tail proof,
not a claim that its all-jet prescription omitted a requirement. The earlier
first-coefficient-plus-half formula was proved for value/first radial inputs,
not for every derivative order. A separate second-coefficient third-jet
constant is needed unless its stronger cutoff budget was also verified.
No cutoff, diagonal schedule, field, frame convention, root, width, tilt,
pressure gauge or historical evidence is changed. These are conditional
analytical estimates, not selected-profile measurements or a full pulse proof.

## 1. Use the actual diagonal without upgrading its quantifiers

At fixed q,eta the native tangential increments are

\[
 e=\sum_{n\ge1}\chi(c_nq)q^{2nh}E_n,\qquad
 u=\sum_{n\ge1}\chi(c_nq)q^{2nh}U_n.                   \tag{1}
\]

The SAME chi satisfies `0<=chi<=1`; it is independent of X,eta. Hence
ordinary native radial derivatives do not differentiate chi. This is NOT
permission to omit its Euler/axial/time or potential-curl derivatives.
The retained diagonal guarantees, for each fixed m<=n, on the WHOLE support,

\[
 C_{n,m}q^{nh}\le2^{-n},                              \tag{2}
\]

where Cn,m includes the required coefficient seminorms. For the statement
here it must bound EACH X derivative of En,Un through order m on the same
closed positive annulus. The original stronger mixed-jet conditions remain.

Set `z=q^h<=1`, `J=max(2,m)`. Supply finite HEAD bounds Hn,m for
n=1,...,J-1; these are norms of the actual globally finalized coefficients.
The actual sums are locally finite at q>0. A positive geometric majorant
for the rest gives, for each displayed jet,

\[
 |\partial_X^j(e,u)|\le
 \sum_{n=1}^{J-1}H_{n,m}z^{2n}
       +{(z/2)^J\over1-z/2},\qquad 0\le j\le m.       \tag{3}
\]

The combined norm means a COMMON component bound, not their Euclidean sum.
Every support point satisfies (2), not merely the cutoff endpoint. In
particular at m=3 the second coefficient belongs to the finite HEAD:

\[
 \boxed{|\partial_X^j(e,u)|\le q^{2h}B_3(z),\quad
 B_3(z)=H_{1,3}+z^2H_{2,3}+{z\over8(1-z/2)},\quad j\le3.} \tag{4}
\]

B3 is increasing in z. Replacing z by a valid upper zstar supplies one
whole-small-q constant; on q<=1 it is at most `H1,3+H2,3+1/4`.
It is not just H1,3+1/2. The n=2, m=3 norm does not follow from the
generic diagonal m<=n. If the SAME retained schedule has ALSO been proved
to satisfy `C2,3*q^(2h)<=1/4` on its support, the stronger tail starting at
two recovers `H1,3+1/[4(1-z/2)]<=H1,3+1/2`. That is an EXTRA verified
input, not a proposal here to alter c2. Additional inherited physical
requirements might already imply it; their mere existence is not checked
as that specific quantitative bound.

For m=1,2 the head has only n=1, so the old half-tail formula is recovered.
For arbitrary fixed m, (3) gives the honest finite-head constant
`sum_(n<J)Hn,m+2^(1-J)` after extracting q^(2h). No bound uniform over m,
or common first-coefficient-only constant for all m, follows.

## 2. Why the second input genuinely cannot be discarded

This is an ABSTRACT scalar profile-budget control, NOT an actual coefficient
of the Navier--Stokes construction. Use the earlier smooth plateau psi,
zero outside (1/4,3/4) and one near X=1/2, with derivative bounds
`P1=256,P2=77824,P3=112640000`. Let

`U2,N(X)=a N^-2 psi(X) sin(N(X-1/2))`, `a=(1+2P1+P2)^-1`.

Every value/first/second X norm is <=1 for N>=1. Its third derivative at
the center is -aN. Put all other coefficients zero, take chi=1 at the
evaluation point and fix a z<=1/2. The order-two PROFILE diagonal budget
`1*z^2<=1/4` passes uniformly; the tail n>=3 is zero. Nevertheless the
third background jet has magnitude `aN z^4` and exceeds `z^2/2` for all
sufficiently large N. A fixed stage's additional cutoff-derivative constants
can be absorbed into another finite low-order constant and a smaller fixed
z; the missing third norm still grows with N. No claim that this control
meets every full vector-field or moment equation is made.

A valid H2,3 is at least the actual third norm; the explicit enclosure
`max(1,a[N+3P1+3P2/N+P3/N^2])` restores (4). Four rational controls check
this algebra. The continuum obstruction is the exact family, not the scan.
It proves only that the low-order diagonal assumptions cannot imply the
proposed third-jet bound. It is not a failure of the actual Borel field.

## 3. Join the first coefficient and retain inverse-radius shear terms

On the SAME native annulus let `0<Rlo<=Rprof=sqrt(2X)<=Rhi`. Suppose
S0,...,S3 bound the finalized first F1,U1 jets there, in particular with the
preceding J3, the unchanged cutoff's L3 and the FULL global repair inverse.
Since E1=Rprof F1, use the product rule with radius derivative bounds

`Rj=(Rhi,Rlo^-1,Rlo^-3,3Rlo^-5)`.

Then `AE1,j=sum_(i<=j)binom(j,i)Ri*S_(j-i)`, `AU1,j=Sj` supplies H1,3
as their maximum (padded by one). No radial factor is frozen.
H2,3 still needs its OWN finalized order-two coefficient bounds. Its pressure
and moments depend on the complete lower-order source, including future
memory; coefficient velocities vanishing on I4 do not erase Pi1 or Pi2.

If (4) bounds each e,u derivative by `B*q^(2h)`, put

`Ij=(Rlo^-1,Rlo^-3,3Rlo^-5,15Rlo^-7)`,
`Fj=B*sum_(i<=j)binom(j,i)Ii`, j<=3.

These bound derivatives of `delta F=e/Rprof`. The exact native shear is
`g=(2X F_X,Rprof U_X)`; hence for r<=2,

\[
 |\partial_X^r\delta g_\theta|
 \le q^{2h}[Rhi^2F_{r+1}+2rF_r],
\quad |\partial_X^r\delta g_z|
 \le q^{2h}B\sum_{j\le r}{r\choose j}R_j.             \tag{5}
\]

The component SUM bounds the Euclidean error. The `-e/Rprof` term in
delta gtheta and ALL its differentiated Jacobians are included through
Fj; they are not replaced by a derivative of e alone. These are native
radial estimates at fixed q,eta. Add valid leading jet bounds to (5) to
obtain TOTAL current jet bounds. A value-only cone test does not supply
those derivative inputs.

## 4. Explicit first and second jets of the CURRENT spectral map

Suppose the previous value-domain gate has independently established

`0<f<=F<=Fmax`, `-gtheta>=t>0`, `|g|<=G`,
`D=|g|^2+2F*gtheta>=d>0`.

Let Fj,gj now denote TOTAL current derivative bounds for j=1,2, not the
error-only factors above. Put `P=-2F*gtheta`, `Pmin=2ft`, `Pmax=2FmaxG`,
`Dmax=G^2+2FmaxG`. The frame is EXACTLY `N=g/|g|`, `K=rot(N)`,
`c=-sqrt(D/P)`. Rational valid magnitude enclosures are
`cmin=min(1,d/Pmax)`, `cmax=max(1,Dmax/Pmin)`.

The normal projection identity gives

\[
 |N_X|,|K_X|\le g_1/t=:N_1,\qquad
 |N_{XX}|,|K_{XX}|\le g_2/t+3g_1^2/t^2=:N_2.          \tag{6}
\]

Indeed `N''=(I-NN^T)g''/|g|-2N'(|g|)'/|g|-N|N'|^2`.
Define the positive derivative caps

`D1=2Gg1+2F1G+2Fmax*g1`,
`D2=2g1^2+2Gg2+2F2G+4F1g1+2Fmax*g2`,
`P1=2(F1G+Fmax*g1)`, `P2=2(F2G+2F1g1+Fmax*g2)`.

For Qc=c^2=D/P, quotient differentiation yields

\[
 Q_1=D_1/Pmin+Dmax P_1/Pmin^2,
\]
\[
 Q_2=D_2/Pmin+(2D_1P_1+Dmax P_2)/Pmin^2
                         +2Dmax P_1^2/Pmin^3.
\]
\[
 \boxed{|c_X|\le Q_1/(2cmin),\qquad
 |c_{XX}|\le Q_2/(2cmin)+Q_1^2/(4cmin^3).}             \tag{7}
\]

The last inverse-cmin cubed cost is essential. The growth rate is
`lambda=2F*(-Ntheta)*|c|`; its first two derivatives are bounded by the
full three-factor product rule in the calculator. No sampled square root
or pressure-free alternative discriminant is substituted.

## 5. Physical radial losses and remaining scope

At fixed physical z,t, q and eta do not depend on r, `X=r^2/(2q)`.
For any of these scalar/vector radial frame components H,

\[
 |\partial_r H|\le q^{-1/2}Rhi H_1,\qquad
 |\partial_r^2 H|\le q^{-1}(H_1+Rhi^2H_2).             \tag{8}
\]

These identities hold at every positive physical scale, given the common
input bounds. They are not scale-loss-free estimates. Band normalization
`s^(-1-h)` leaves N,K,c unchanged; it MULTIPLIES lambda and its radial
derivatives, since s is constant in a radial derivative. Frozen representative
data remain frozen: this does not replace the original primary covariance
frame by a different current one. Phase-normal frames, axial/time/angular
Cartesian derivatives, complete potential curls and generated nonlinear
frames still require their own hypotheses. The derivative bounds here add
no derivative-dependent change to the original first-band choice.

Symbolic checks verify radius/shear, quotient, normal and physical-chain
identities. Five manufactured spectral profiles test the estimates; four
second-coefficient controls test the finite-head issue. They are not the
selected enormous fields or continuum proofs. Complete pressure/stress
budgets remain separate inputs; bounding tangential velocity alone does
not bound the full relative stress or covariance.

Next obtain H2,3 (or verify the extra SAME-schedule n=2 third-jet budget),
the unchanged cutoff shape jets, and actual value/scale acceptance. Then
join the remaining phase/frame/potential operators and full covariance,
axial-pressure mean and viscous-cutoff-interaction residual obligations.
No unconditional pulse, full residual improvement or blow-up is proved.

Run `python3 test_audit.py -v` and
`python3 audit.py --verify-record evidence.json --output /fresh/path.json`.
