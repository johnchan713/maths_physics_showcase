# Explicit leading-reference frame freezing with a pointwise target

Parent: `b0f179b085aa0a64f281fd9ada47951c28e0e196`.
Status: `explicit-leading-frame-freezing-from-global-jets-with-pointwise-full-target`.

This checkpoint bounds the previously symbolic representative-freezing
constant for the LEADING reference frame, from the retained global F/U jets
and endpoint value floors. Crucially, the target stress remains at the
evaluation point. Neither a derivative of its unit direction nor a positive
stress-amplitude floor is needed. The global pressure-sensitive orientation
and cone assumptions ARE still needed. The bounds below do not independently
prove those inherited assumptions, or cover an arbitrary generated nonlinear
frame. No cutoff, carrier, root, width, support, gauge or mean is changed.

The defining comparison is the supplied manuscript (7.24), (7.28)-(7.29),
printed page 83; SHA256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The [primary covariance foundation](../primary_covariance_foundation/README.md)
already keeps this pointwise target and a separate freezing gate. Its words
"target-freezing" must not be read as replacing the stress by its value at a
representative. The historical evidence is unchanged.

## 1. Relative coordinate comparison, not stress-amplitude freezing

Write g=(gtheta,gz), N=g/|g|, K=(-Nz,Ntheta), c=|c0|>0, and
w=sqrt(1+u^2) for the SAME retained u>0. For the complete target
T(p)=m(p)n(p), |n|=1, m>0 in the active interior, set

\[
 H=-N(p)\cdot n(p)\ge\omega>0,\qquad
 V_0={mH\over cw},\quad W_0={mK(p)\cdot n(p)\over u}.
\]

Freeze ONLY N,K,c at p*, and put Hf=-N(p*) dot n(p).
The same m(p), n(p) appear in Vf=mHf/(cf w), Wf=m K(p*) dot n(p)/u.
Thus, exactly,

\[
 {V_f-V_0\over V_0}
 ={c\over c_f}{H_f-H\over H}+{c-c_f\over c_f},\qquad
 {W_f-W_0\over V_0}
 ={cw\over uH}(K(p_*)-K(p))\cdot n(p).                 \tag{1}
\]

There is no difference n(p*)-n(p): the target is NOT frozen. Both target
components are retained; in particular no axial/pressure contribution is
declared zero because the local axial profile happens to vanish. The complete
pressure memory still enters the global orientation, cone and Cstar bounds.
At zero-stress edges use the inherited limiting direction and weighted
extension; (1) is proved on the open interior, without division by zero.

Suppose BOTH endpoints have

\[
 f\le F\le\bar F,\quad -g_\theta\ge t>0,\quad |g|\le\bar G,
 \quad c_-\le c\le c_+,\qquad
 |\Delta F|\le J_F d,\quad |\Delta g|\le J_g d.           \tag{2}
\]

Set P=F(-gtheta), D=|g|^2+2F gtheta, so c^2=D/(2P). Endpoint
product/quotient differences, not derivatives of a square root, give

\[
\begin{split}
 \bar D&=\bar G^2+2\bar F\bar G,\\
 L_D&=2(\bar G+\bar F)J_g+2\bar GJ_F,&
 L_P&=\bar GJ_F+\bar FJ_g,\\
 L_C&={L_D\over2ft}+{\bar D L_P\over2f^2t^2},&
 L_N&=2J_g/t,\\
 |\Delta(c^2)|&\le L_Cd,&|\Delta N|=|\Delta K|&\le L_Nd.
\end{split}                                                    \tag{3}
\]

For example, use
`D/P-D*/P*=(D-D*)/P+D*(P*-P)/(P P*)`, and P,P*>=ft.
Normalization of g uses |g|,|g*|>=t. Therefore
|Delta c|<=LC d/(2c_-). Equation (1), w<=1+u, implies

\[
\boxed{
 L_V={c_+L_N\over c_-\omega}+{L_C\over2c_-^2},\quad
 L_W={c_+(1+u)L_N\over u\omega},\quad
 \max(|V_f-V_0|,|W_f-W_0|)\le\max(L_V,L_W)dV_0.}        \tag{4}
\]

The use of 1+u is solely an upper bound, not a physical replacement of w.
No norm or derivative of T/zeta is assumed bounded. In particular the known
outer inverse-distance losses in primary roots are not removed.

## 2. Explicit jets on a controlled slow-coordinate neighborhood

Use Xa=a^2, Xb=b^2, a>0. Leading F and g have the same positive band
factor s^(-1-h), which cancels in N,K,c. Hence work with the native profiles
F(X,eta), g=(2X F_X,sqrt(2X) U_X). The leading-frame foundation supplies
endpoint floors: `a_log=-2 partial_y log F>=alpha=1/(2G^3)`,
`v=a_log+b_log^2/a_log>=2+nu`, `nu=1/(2G^4)`, and a finite upper Vb.
Thus t=alpha f* and c_- =1/(2G^2) are valid, with any proved c_+ dominating
sqrt((Vb-2)/2), e.g. Vb. These are profile floors, not automatically the
band-factor-inclusive perturbation floors in a different comparison.

For an active representative on 1/4<=s<=4, Xa<=X<=Xb, control the larger
rectangle 1/8<=s<=8, Xa/2<=X<=2Xb. In slow coordinates (R,Z,T),

\[
 s-Z^2s^{2h}=T,\quad X=R^2/(2s),\quad
 \eta=Zs^{-D_h},\quad D_h=1/2-h,\quad A_h=1/2+h.
\]

On T>=0, |eta|<=1 and L=1-2h eta^2>0.98 for 0<h<.01.
Implicit differentiation gives s_T=1/L, s_Z=2 eta s^A_h/L, and, INCLUDING
the cancellation at eta=+/-1,

\[
 \eta_Z=s^{-D_h}(1-\eta^2)/L,\qquad
 \eta_T=-D_h\eta/(sL).
\]

Consequently |s_T|<=2, |s_Z|<=16, |eta_T|<=8, |eta_Z|<=6.
On the enlarged rectangle sqrt(2X) lies in [a,2b]. With Xhi=2b^2,

\[
 K_X=8b+144Xhi,\quad K_\eta=14,\qquad
 |\nabla_{slow}X|_1\le K_X,\quad |\nabla_{slow}\eta|_1\le K_\eta.
\]

The positive-chart continuation bootstrap is valid whenever the straight
segment distance d satisfies

\[
 d\le d_{chart}=\min(1/288,Xa/(4K_X)).                    \tag{5}
\]

Indeed |Delta s|<=18d<=1/16 and |Delta X|<=KX d<=Xa/4 before any first
bootstrap exit. Starting in the active rectangle these remain STRICTLY
inside the enlarged rectangle. The implicit derivative is nonsingular, so
there is no first exit. Both endpoints and the segment must have T>=0;
this argument does not invent a negative-time heat extension.

Let the global swirl input bounds be f3,fx2,fxx1 and the axial ones ux2,uxx1.
The swirl checkpoint's native mixed envelope holds for ALL X>=Xrep, not
only up to its chosen Xc. Its conversion plus compact/axis estimates
therefore applies up to the enlarged finite endpoint 2Xb as well: call its
`swirl_inputs(Xa,2Xb,...)`. This is an explicit extension of that interval,
not an assumption that the old Xc equals 2Xb. The axial profile and its jets
extend by zero beyond the original Xv after the exact moment closure;
this does NOT extend the pointwise pressure source by zero.

Native differentiation gives

\[
\begin{split}
 G_X&=2fx2+2Xhi fxx1+ux2/a+2b\,uxx1,\\
 G_\eta&=2Xhi fx2+2b\,ux2,\\
 J_F&=fx2 K_X+f3 K_\eta,\qquad
 J_g=G_X K_X+G_\eta K_\eta.                              \tag{6}
\end{split}
\]

The inverse-radius term ux2/a in g_z,X is essential. The rectangle in
(X,eta) is convex, so the global raw-jet bounds prove the endpoint
differences in (2); positivity of c need not be assumed at intermediate
points of an arbitrary slow segment. Its value floors are used at the two
endpoints only. No derivative of the representative or a discrete label is
taken.

## 3. Retained boxes, explicit freezing constant and joined gates

The auxiliary-geometry audit fixes delta=ell^-6 and closed slow-box
halfwidth delta. A representative is in its original box. A point in the
k-fold box is at distance at most (k+1)sqrt(3)delta. Thus safe rational
constants are Dbox=4 for the original and Dbox=6 for the doubled box.
For larger fixed enlargements use Dbox=2(k+1), not an untracked constant.
Intersect with T>=0 when applying the positive chart argument above.

Combine (4)-(6):

\[
\boxed{C_{freeze}=Dbox\max(L_V,L_W),\qquad
 \tau_{freeze}\le C_{freeze}\ell^{-6}.}                  \tag{7}
\]

This removes the need to estimate a derivative of the flat limiting stress
direction for THIS frame-only step. Its pointwise orientation floor remains
indispensable. With a valid unfrozen same-u margin gamma0>0, the existing
primary comparison gives a frozen margin at least gamma0/2 if

\[
 \boxed{\ell^6\ge\max(Dbox/d_{chart},8C_{freeze}/\gamma_0).} \tag{8}
\]

The exact half-margin threshold is gamma0/(4-gamma0); the padded gamma0/8
is retained. Original strict leading admissibility with raw deficit d0
also supplies a leading margin at least d0/(2 kappa_u), hence the smaller
`gamma_B=omega^2 d0/(9 kappa_u)` from the previous checkpoint is safe for
the leading target. Composing its complete Cstar background requirement
with (8) yields the independently necessary-to-this-certificate gates

\[
 E=C_\star(M,N)2^{-2h(\ell-2)}\le r_b,\quad
 R_{act}+\max(2L_\Psi E/\omega^2,
       72\kappa_u C_{freeze}/(\omega^2\ell^6))\le\kappa_u,
 \quad\ell^6\ge Dbox/d_{chart}.                          \tag{9}
\]

This is a CONDITIONAL composition, not acceptance of the retained first
band. In (9), Cfreeze must correspond to the frame/target actually being
frozen. Our derived constant applies to the leading-reference frame. A
perturbed/nonlinear frame needs its own complete jets and endpoint floors;
its target orientation may be only omega/3, not omega. Passing the leading
test cannot silently stand in for that further obligation. Pulse length,
column-direction error, carrier normal, actual inverse and support gates
are still separate. No old claim that another lower-band selector dominates
Cfreeze is restored by merely proving this constant finite.

## 4. A flat-amplitude countercontrol and reproducibility

It would be invalid to freeze T itself. For m(y)=exp(-1/y^2), choose points
y and y+y^3, with identical constant frame and direction. Then

\[
 {m(y+y^3)\over m(y)}
 =\exp\left({2+y^2\over(1+y^2)^2}\right)\longrightarrow e^2,
 \qquad |(y+y^3)-y|=y^3\longrightarrow0.
\]

Thus relative stress-amplitude freezing is not O(distance), although the
true frame-freezing error is exactly zero. This is a manufactured analytic
countercontrol, not a Navier-Stokes field and not an allegation that the
original source froze the amplitude.

The evidence includes 216 exact rational endpoint comparisons with both
target components nonzero, three unchanged tilts, eight perturbation sign
choices, three distances, and magnitudes 1,10^-40,10^-200. It checks the
mixing and relative-coordinate bounds and exact magnitude cancellation.
Symbolic controls check pole cancellation, native shear derivatives, (1)
and the flat countercontrol. A rational auxiliary composition uses the
complete global jets, doubled-box cost, full pressure/viscosity rows and
same-deficit witness. These fixtures do not identify the actual huge profile,
cutoff norms, raw deficit, retained first integer or selected physical pulse.

Run `python3 test_audit.py -v`, then
`python3 audit.py --verify-record evidence.json --output /fresh/path.json`.
Evidence binds new sources and recursively verifies historical bytes with
the original comparators. Tests and symbolic checks are not a continuum
proof. Remaining: evaluate valid unchanged-cutoff/raw-deficit inputs and
all retained-scale gates, extend freezing to any required current frame,
then complete uniform actual pulses, full covariance/mean corrections and
the complete viscous cutoff/interaction residual.
