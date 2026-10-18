# Explicit phase-to-frame and primary-column error constants

Parent: `255ae2e44292419235dd9b7e1a340f20e690dc39`.
Status: `explicit-current-phase-frame-and-primary-column-comparison-budgets`.
The parent hosted CI was still progressing without a failed step when
inspected on 18 October 2026. All historical source/evidence bytes remain.

This supplies the next input of the Gaussian/Haar checkpoint: explicit
frame errors CE,Cd and primary-column constants D0,D1,Cdir. The proof uses
the actual phase, specified frame left inverse and complete moving-normal
pressure projection. It replaces an appeal to smooth dependence by a
two-coordinate identity and displayed inequalities. No roots, tilt,
carriers, radius, width C_setup^-20, pressure gauge or supports are changed.

The result is relative to the reviewed leading/core and background jet
bounds. It makes their finite seminorm dependence explicit; it does not
numerically materialize the enormous selected seminorms or first band,
close every nonlinear covariance/mean/residual input, or prove blow-up.
Exact algebra and manufactured phase controls are not continuum proofs.

**Frame-convention correction.** Sections 2--4 below apply to the inherited
repository's actual-normal mixing convention, NOT verbatim to manuscript
(7.8). That equation uses the reference s(v) inside its mixing matrix,
while both conventions use the actual sa inside U. Section 5 supplies the
exact bridge, a separate manuscript-frame estimate, and the entrance-data
obstruction. No historical pulse or evidence is silently replaced.

## 1. Primitive errors from the original phase and background jets

Write L=Ls, S=ell^2, epsilon=Q^h, k=ceil(epsilon^-1/2), and retain
`1<=L<=S`, `S^2(epsilon+epsilon^2+k^-1)<=1`, `1<=epsilon*k^2<=4`.
Let R,R0 belong to [Rmin,Rmax], Rmin>0, and |R-R0|<=Dbox*S^-3.
The current frozen leading representative gives F0, g0=(R0 F_R0,G_R0),
N=g0/|g0|, K=(-Nz,Ntheta), lambda0 and c0. Use the current bounds

`0<F0<=Fmax`, `gmin<=|g0|<=gmax`, `cmin<=|c0|<=cmax`,
`lambda_minus<=lambda0<=lambda_plus`, `bmin<=Bs<=bmax`.

The spectral floors provide these. Bs is the original calibrated amplitude:
`epsilon*k^2*Bs^2=lambda0/(1+u^2)^(3/2)`. The background foundation supplies
one finite Db bounding each of F-F0, F_R-F_R0, G_R-G_R0 by
`Db*(S^-3+epsilon^2)`, and one finite Mz bounding |F_Z|,|G_Z|.
These are norms of the actual fixed background, not manufactured values.
Db can be chosen from its C2 leading bound, positive-order C1 comparison
and fixed box constant; all are uniform on the accepted closed charts.
Use `Frmax=gmax/Rmin`, and define

\[
 Pmax=Rmax\,bmax(1+u/gmin)+1,\quad
 Zmax=bmax(1+u/gmin),\quad Zc=(Pmax+Zmax)Mz,
\]
\[
 De=2(Pmax+Zmax)Db+Frmax,
\quad Dr=De,\quad Dr'=bmax\,u+De,
\]
\[
 Dt=bmax\,Dbox/Rmin+bmax\,uRmax/(Rmin\,gmin)
       +1/Rmin+bmax\,u/gmin+Zc,\quad Dt'=Zc,
\]
\[
 DF=2Db,\qquad Dg=2(Rmax+1)Db+Dbox\,Frmax.             \tag{1}
\]

For the SAME rounded p and pz as the parent, set W=pF+pzG and
`n=(x0-v W_R,p/R,pz-epsilon*v W_Z)`.
Its reference is `nref=Bs(sref,K)`, `sref=sigma*u*(1/2+v/L)`.
The exact cancellation `ptilde F_R0+pz G_R0=-sigma Bs u/L` gives

\[
 |n_r-Bs\,sref|\le Dr/L,\quad |n_t-Bs K|\le Dt/L,
 \quad |n_r'|\le Dr'/L,\quad |n_t'|\le Dt'/L,         \tag{2}
\]
\[
 |F-F0|\le DF/L,\qquad |g-g0|\le Dg/L,
 \quad g=(RF_R,G_R).                                 \tag{3}
\]

To check (2), the radial-derivative error is at most De/S^2; multiplication
by v<=L gives De/L. Angular normal errors are the R/R0 displacement,
the deliberate tilt and the full nonzero-integer rounding error 1/k.
The axial error includes epsilon*L*Zc<=Zc/L. The pulse derivative includes
the SAME axial term. No L/k term, epsilon*v term or angular inverse is lost.
The two S^-3/epsilon^2 contributions in (3) each fit S^-2<=L^-2.

## 2. The exact two-coordinate projected generator

Put t=|n_t|, a=n_t/t, b=(az,-atheta), s=n_r/t, w=sqrt(1+s^2),
and omega=a'.b. Write

`U=[er-s a,b]`, `M=[[1,1],[c0*w,-c0*w]]`, `B=U M`,
`Bleft=M^-1[er^T;b^T]`.

This is the specified ambient left inverse, not an assumed Moore--Penrose
inverse. With
`Kbg=[[0,-2F,0],[2F+RF_R,0,0],[G_R,0,0]]` and
`An=-Kbg+n(n^T Kbg-nprime^T)/|n|^2`, direct multiplication gives

\[
 H=[er^T;b^T](An U-U')
 =\begin{pmatrix}
 {s(a\cdot g-s')\over1+s^2}&{2F b_\theta-s\omega\over1+s^2}\\
 -(2F b_\theta+b\cdot g)+s\omega&0
 \end{pmatrix},
\quad Bleft(AnB-B')=M^{-1}HM-M^{-1}M'.                \tag{4}
\]

In particular the pressure contribution nprime and frame derivative B'
are both present. Also

\[
 M^{-1}M'={ss'\over2(1+s^2)}
        \begin{pmatrix}1&-1\\-1&1\end{pmatrix}.       \tag{5}
\]

For a=K, b=N, g=g0, omitting only the O(1/L) connection entries defines
`H0(s)=[[0,2F0*Ntheta/(1+s^2)],[-(2F0*Ntheta+|g0|),0]]`.
The ORIGINAL c0 exactly diagonalizes it:
`M^-1 H0(s) M=diag(lambda0/w,-lambda0/w)`.

## 3. Explicit Euclidean operator-norm CE and damping Cd

Assume the original finite value-acceptance inequalities
`L>=max(1,2Dt/bmin,Ds)` below. Define

\[
 beta=bmin/2,\quad Sref=3u/2,\quad S0=Sref+1,
 \quad Da=2Dt/bmin,\quad Ds=(Dr+Sref\,Dt)/beta,
\]
\[
 Omega=Dt'/beta,\qquad S1=(Dr'+S0\,Dt')/beta.         \tag{6}
\]

Normalization gives t>=beta, |a-K|<=Da/L, |s-sref|<=Ds/L,
|s|<=S0, |omega|<=Omega/L and |s'|<=S1/L. Let

\[
 D11=S0(Da\,gmax+Dg+S1),
\quad D12=2DF+2Fmax\,Da+S0\,Omega,
\]
\[
 D21=2DF+2Fmax\,Da+Dg+gmax\,Da+S0\,Omega,
\quad DH=D11+D12+D21,
\]
\[
 kM={\max(1,cmax(1+S0))\over\min(1,cmin)},\quad
 CE=kM\,DH+S0\,S1+lambda_+\,Ds.                     \tag{7}
\]

The entries of H-H0(s) are bounded by Dij/L. The two singular values of
M are sqrt(2) and sqrt(2)|c0|w, so its condition number is at most kM.
Equation (5) costs at most S0*S1/L. Finally lambda0/sqrt(1+s^2) changes
by at most lambda_+*Ds/L when s is replaced by sref. Consequently

\[
 \|Bleft(AnB-B')-diag(\lambda,-\lambda)\|_2\le CE/L,
 \qquad \lambda=\lambda0/\sqrt{1+sref^2}.             \tag{8}
\]

This is an OPERATOR-norm bound; an entrywise CE is not substituted for it.
For damping put Dn=Dr+Dt and Nref=bmax(1+Sref). Then

\[
 Cd=4Dn(2Nref+Dn),\qquad
 |epsilon\,k^2|n|^2-dref|\le Cd/L,\quad
 dref={\lambda0(1+sref^2)\over(1+u^2)^{3/2}}.         \tag{9}
\]

The factor 4 is the original rounding/calibration bound, not a new viscosity.
All constants above are explicit finite functions of reviewed background
seminorms and current floors. They are not evaluated at the giant actual
selected parameters. Once these value floors hold, every fixed derivative
of (4)--(9) follows by product/inverse rules with finite polynomial S and
edge losses. No derivative-dependent first band is introduced.

## 4. The complete principal-column ratio and its remaining connection error

Use the original Riccati acceptance L>=8CE/mu, where
`mu=lambda_minus/sqrt(1+Sref^2)>0` (or its conservative lower bound
`lambda_minus/(1+Sref)`). With q=zminus/zplus, |q|<=CE/(mu L)<=1/8.
For x=zplus+zminus, direct use of BOTH frame coordinates gives

\[
 {t_{tan}\over x}=-s a+c0\,w\,{1-q\over1+q}\,b.      \tag{10}
\]

The factor (1-q)/(1+q) cannot be replaced by one. It differs from one
by at most 16CE/(7mu L). Let Wref=1+Sref, W0=1+S0. Then

\[
 D0=Ds+Sref\,Da+cmax\{Ds+Wref\,Da+16W0\,CE/(7mu)\},
 \qquad D1=u(1+cmax).                                \tag{11}
\]

Equations (6),(10) imply a D0/L error from the ideal instantaneous
direction `-sref K+c0 sqrt(1+sref^2)N`. That direction differs from its
midpoint value by at most D1*|v-L/2|/L. The previous Gaussian/Haar result
therefore gives the explicit normalized Haar-column bound

\[
 |e_\sigma|\le Cdir/\sqrt L,\quad
 Cdir=D0+D1\,{x_+^2e^{C_G}\over2x_-^2c_G}.            \tag{12}
\]

Here x_-,x_+,c_G,C_G are the SAME supplied Gaussian constants. The original
Haar mass, angular factor 1/2 and physical area r0^2/sqrt(L) stay intact.
Only principal columns define the original positive roots; complete curl
covariance remainders still enter the later full tensor/mean system.

Insert (12) into the EXISTING cone budget, for example
`L >= [4*Cdir*max(Amin^-1,u^-1)/gamma]^2`, with the existing freezing
condition and L>=25. This is a finite sufficient value inequality, not a
claim that an arbitrary first numeric band satisfies it. Its potentially
enormous constants are best retained logarithmically. It exposes the
remaining acceptance check without changing the construction or root.

## 5. Correction: manuscript mixing and inherited replay are not identical

The manuscript with the SHA256 bound in protocol.json, equation (7.8), uses
`Br=U Mr`, `Mr=[[1,1],[c0*wr,-c0*wr]]`, `wr=sqrt(1+sref^2)`.
The inherited `physical_pulse_foundation/README.md` and the actual basis
code in `geometric_pulse_replay/diagnostics.py` instead use
`Ba=U Ma`, `wa=sqrt(1+sa^2)`. Here sa is s in (4)--(11).
Both are valid frames for the SAME plane. They are not the same frame.
Statements describing the inherited actual-mixing replay as literally
the original manuscript frame must be read with this correction.

For either frame (4) remains exact, with its own M and Mprime. In the
manuscript convention, comparison to H0(sref) gives

\[
 kMref={\max(1,cmax(1+Sref))\over\min(1,cmin)},\qquad
 CEref=kMref(DH+2Fmax\,Ds)+Sref\,u,                  \tag{13}
\]
\[
 \|Brleft(AnBr-Br')-diag(\lambda,-\lambda)\|_2
 \le CEref/L.                                      \tag{14}
\]

Indeed H-H0(sa) costs DH/L; the derivative of 1/(1+s^2) has magnitude
at most one, so H0(sa)-H0(sref) costs at most 2Fmax*Ds/L.
Mr diagonalizes H0(sref) exactly and its connection costs
`|sref*sref'|/(1+sref^2)<=Sref*u/L`. Damping (9) is unchanged.
This proves an explicit version of the manuscript's operator comparison,
not an equality of its primary pulse and the inherited one.

The exact change of coordinates is

\[
 T=Ma^{-1}Mr={1\over2}
 \begin{pmatrix}1+wr/wa&1-wr/wa\\1-wr/wa&1+wr/wa\end{pmatrix},
 \quad Br=BaT,\quad za=Tzr,                         \tag{15}
\]
\[
 Gr=T^{-1}GaT-T^{-1}T',\qquad
 fr=T^{-1}fa,\qquad zr(0)=T(0)^{-1}za(0).           \tag{16}
\]

These formulas include scalar viscous damping (which commutes with T),
the complete forcing coordinates and the frame derivative. For matching
initial data, uniqueness gives the same physical t, and therefore the same
normal pressure and covariance. In particular zero-data forced inverses
are invariant under this FULL transformation. No pressure memory is lost.
For the inherited homogeneous primary with za(0)=P0*(1,0), however,

\[
 zr(0)={P0\over2}(1+wa(0)/wr(0),\ 1-wa(0)/wr(0)).   \tag{17}
\]

Resetting zr(0)=P0*(1,0) generally changes the physical primary. A precise
obstruction appears already before rounding, at v=0, R=R0, epsilon*v=0:
the original deliberate tangential tilt is
`nt=Bs*(K-sigma*u*g0/(L*|g0|^2))`, with K.g0=0, nr=Bs*sigma*u/2.
Hence, for finite L and u>0,

\[
 wr(0)^2-wa(0)^2={u^2\over4}
 {u^2/(L^2|g0|^2)\over1+u^2/(L^2|g0|^2)}>0.         \tag{18}
\]

Thus the required equality wa(0)=wr(0) is NOT an identity of the phase
construction. Rounding may alter the difference; it does not provide a
universal equality. Even arbitrarily small nonzero entrance differences
cannot be discarded in an EXACT identification of pulses or their roots.
From (6), |wa-wr|<=Ds/L, so the transported unwanted-to-growing entrance
ratio has magnitude at most Ds/(2L); this is small, not zero.

This checkpoint preserves the existing repository pulse and its roots.
Its Gaussian/column bounds in sections 3--4 concern that pulse. Equation
(14) separately supports the manuscript frame. Before claiming their
canonical growing-eplus primaries or root values coincide, one must either
explicitly carry (17) through every relevant manuscript identity, or
construct and validate the manuscript-canonical pulse/root family anew.
An O(1/L) comparison alone does not supply that exact identification.
Neither alternative is asserted completed here. This is an additional
identified dependency, not an unconditional pulse/residual theorem.

## 6. Evidence and next obligation

The analytical proof is (1)--(12) on the whole accepted domain.
`frame_budget.py` implements exact rational budgets. `frame_identity.py`
independently forms An, B, Bprime and the left inverse with Fractions,
checks (4)--(5), and detects separate omissions of nprime, Bprime and the
column ratio. It also checks the time-dependent transition (15)--(17),
transformed forcing, and the entrance-reset counterexample with exact
Fractions. Equation (18) has a positive exact squared-scale fixture.
A separate high-precision phase control uses the original
integer rounding, full epsilon*v axial normal and all pulse endpoints.
Its finite parameters are manufactured controls, not the selected giant
background or a continuum supremum certificate.

The complete pressure `i(n.Kbg*t-nprime.t+n.f)/(km|n|^2)`, Pi1/Pi2 memory,
axial mean velocity and moving pressure moment remain. Gaussian linear
tails remain additive; actual cutoff/A0/A1 velocities stay in products.
First resolve the canonical manuscript entrance versus inherited mixing
dependency above; then join these explicit frame/column budgets to the primary
roots, generated complete covariance, fresh mean maps and the complete
viscous/cutoff/interaction residual on one current domain. Retained core/
background analytical dependencies must still be discharged before any
unconditional construction or blow-up claim.

```sh
python3 research/navier_stokes_cascade/results/frame_column_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/frame_column_foundation/audit.py \
 --verify-record research/navier_stokes_cascade/results/frame_column_foundation/evidence.json \
 --output /tmp/frame_column_fresh.json
```
