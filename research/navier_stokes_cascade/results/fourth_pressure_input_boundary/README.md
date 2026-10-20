# Fourth moving-phase transfer: the fifth-moment input boundary

Parent: `187149d4cdd64da1d82d405dc3cfbfd3c92519f5`.
Status: `exact-fourth-pressure-principal-rows-and-moving-phase-input-boundary`.

The parent reduced the compact C4 certificate to the actual primitive jets
a4,b4. This checkpoint resolves WHICH higher native inputs those jets need.
It extracts the exact fifth-moment combinations in the fourth pressure
coordinates, retains their lower-jet remainders and all incoming pressure,
and gives the full fourth moving-phase transfer into the primitives. An
exact compatible-field control proves that a C4 native/moment bound ALONE
cannot be promoted to a C4 pressure-coordinate bound. This is a limitation
of that inference, NOT evidence that the selected actual profile fails.

All roots, cutoffs, radial widths C^-20, supports, N=1+floor(H^32), pressure
gauges and historical evidence are unchanged. The actual higher inputs and
the compact gate are NOT evaluated. No global H2,3, physical-pulse/covariance/
residual or blow-up result follows. Symbolic and manufactured tests below
check algebra, not a continuum existence proof.

## 1. Exact pressure rows, with the original normalizations

At a fixed X>0 let d=1-eta^2, D=1/2-h, Aphys=1/2+h,
L=1-2h eta^2, Hx=sqrt(2X)E. The defining forward map is the unedited
`compact_pressure_bridge/bridge.py:forward`, also transcribed in
`compact_jet_envelope/README.md` section 5:

\[
\begin{split}
 W&=1-(2D\eta M+dM_\eta)/X,\\
 Q&=-W+[(1-h)I-D\eta I_\eta-dJ_\eta+2(h-D)\eta J]/(XH_x),\\
 N_s&=-WU+[D(M-\eta M_\eta)+4h\eta S-dS_\eta]/X
                  +4Aphys\eta\Pi-d\Pi_\eta,\\
 p_1&=XQ/L,\qquad p_2=XN_s/(LE),\qquad\Pi=\Pi_{axis}+C_p.
\end{split}                                                     \tag{1}
\]

Write f_j=partial_eta^j f/j! at any base point eta0, and m5,i5,j5,s5,
c5,k5 for the NORMALIZED fifth jets of M,I,J,S,Cp,Pi_axis.
In this section E,U and every denominator retain their full jets through
four. Define R1,R2 by computing the fourth coefficient of (1) with ONLY
those six fifth moment/datum coefficients set to zero. Then EXACTLY

\[
\begin{split}
 (p_1)_4&=R_1+{5\over L}\left[d m_5-
                    {D\eta i_5+d j_5\over H_x}\right],\\
 (p_2)_4&=R_2+{5\over LE}\left[(dU-D\eta)m_5-ds_5
                                  -Xd(c_5+k_5)\right].
\end{split}                                                     \tag{2}
\]

Every coefficient outside a fifth jet in (2) is its VALUE at eta0.
The factor 5 is 5!/4!, not an optional convention. To prove (2), the
only fifth inputs in (1) arise when all four derivatives hit M_eta,
I_eta,J_eta,S_eta or Pi_eta. All other Leibniz terms use moment jets at
most four and belong to R. In particular -WU contributes +dU m5 in
the second row; omitting axial U changes the highest-order combination.
Differentiated L,E,Hx and d are still present in R. The program uses an
exact rational Taylor algebra and calls the ORIGINAL forward implementation,
then compares with the independently written principal rows in (2).

Consequently, define the two explicit combinations in brackets as T1,T2.
A sufficient actual input certificate is

\[
 \sup|(p_1)_4|\le\sup|R_1|+5\sup|T_1/L|,\qquad
 \sup|(p_2)_4|\le\sup|R_2|+5\sup|T_2/(LE)|.          \tag{3}
\]

These are combined rows: individual fifth moments need not be bounded
separately if their ACTUAL combinations can be controlled. But neither
combination may be declared zero without an identity. At eta=+-1,
d=0 does NOT remove the I5 term of p1 or the M5 term of p2. The same
axis datum has k5; it does not vanish merely because a modulation density
has stopped. The earlier same-gauge cancellation applies to DIFFERENCES,
not to the original pressure data needed to build the loop.

## 2. A genuine derivative loss, not a missing quadrature refinement

On 1<=X<=2 take E=1 and

`U=epsilon*n^(-4)*sin(n*(eta-1/2))`, n>=1.

With fixed zero incoming moment data at X=1, the EXACT moments at X=2 are

`M=U, I=(8-2sqrt(2))/3, J=I*U, S=U^2-1/2, Cp=log(2)/2`.

Keep the same constant axis datum. E,U and ALL five moments have uniform
Ceta4 bounds: |U|4<=(65/24)epsilon, and the algebra bound controls U^2.
Their fixed-order radial jets on this rectangle are uniform too. The
fifth moment m5 at eta0=1/2 is epsilon*n/120. Direct expansion of the
COMPLETE (1), subtracting the U=0 baseline, gives

\[
 \lim_{n\to\infty}{\Delta(p_2)_4\over n}
       =-{D\epsilon\over48(1-h/2)}\ne0.                \tag{4}
\]

Thus no bound depending only on fixed C4 input norms, geometry and these
positive E,L floors controls p2 in C4. Lower terms, the U-square in S,
the changing W and denominator derivatives are all included in (4);
they do not cancel its linear growth. This is a manufactured patch, not
the selected profile, a global axis construction or an accepted cone family.
It refutes only an invalid norm-transfer inference. An actual structural
identity or additional quantitative derivative input can close the gap.

Four exact rational comparisons use h=1/200, epsilon=1/10 and
n=8,32,128,512. In the p2 DIFFERENCE comparator only, I,J can be omitted
because p2 does not use them, and the unchanged constant Cp/axis terms
cancel. This does not assert zero physical pressure or reset any actual
incoming moment. The symbolic limit (4), not these four samples, proves
the non-implication.

The inherited higher_jet_budget supplies original fields/moments J4 and
loop inputs J3. Its native core derivatives through six and reference
angular C5 estimates are useful seeds, NOT a global pre-modulation C5
bound through all joins and exterior intervals. Also
`a=1-2E_y/E, bs=2U_y/E` at angular order four use E_(y eta^4),
U_(y eta^4), of total order five. A global J4 estimate alone does not
supply those mixed derivatives either. Prescribed-slope identities may
bound them without asking for every J5 derivative; that must be shown.

## 3. Full fourth moving inverse, with no phase freezing

The SAME signed variance root has F(mu,p)=s, p=p2,
`F=d0*mu*sqrt(R2(mu*p))`, and s is the retained flat-cutoff target.
At a base point use normalized jets. If R_F is the fourth Taylor
coefficient of F(mu_<4,p_<4), then

\[
 \mu_4={s_4-F_p p_4-R_F\over F_\mu}.                 \tag{5}
\]

The accepted positive F_mu is used, including at mu=0 where
F_mu=d0/sqrt(2); never divide by the vanishing V_mu. Partial derivatives
of F through four, hence derivatives of q through FIVE in its removable
integral formula, need quantitative estimates. The old q-through-four
ledger was for loop order three and is not relabelled. p1 fixes cone/cap
data, but is NOT a separate varying argument of this variance equation;
only its p2 fourth jet occurs in (5).

Let `Phi(eta,theta)=int_0^theta w(eta,r)dr`,
`w=a(1+t^2)/(2*pi*v)>0`, with the SAME mean-one lift, and
`Phi(eta,Theta(eta,varphi))=varphi`. At a fixed external phase let
qj=partial_eta^j Theta/j!, and write Fij for raw mixed derivatives
partial_eta^i partial_theta^j F at the base point. Define

\[
\begin{split}
 {\cal R}_4(F)={}&F_{11}q_3+F_{02}(q_1q_3+q_2^2/2)
 +F_{21}q_2/2+F_{12}q_1q_2+F_{03}q_1^2q_2/2\\
 &+F_{31}q_1/6+F_{22}q_1^2/4+F_{13}q_1^3/6+F_{04}q_1^4/24.
\end{split}                                                     \tag{6}
\]

Taylor multiplication proves `(F(eta,Theta))_4=F40/24+F01*q4+R4(F)`.
Therefore for either loop source
`gA=-(aL-a)/2`, `gB=E*(bL-bs)/2`, with all t-dependence retained,

\[
 q_4=-{\Phi_{40}/24+{\cal R}_4(\Phi)\over w},\qquad
 (g(eta,Theta))_4={g_{40}\over24}-{g_{01}\over w}{\Phi_{40}\over24}
            +{\cal R}_4(g)-{g_{01}\over w}{\cal R}_4(\Phi).
                                                               \tag{7}
\]

The fixed-angle coefficient Phi40/24 is the FULL integral of w40/24
from 0 to theta. The remainders in (6) use at most three slow derivatives,
but include FOUR phase derivatives. They can be bounded from lower slow
jets and the scalar loop functions, not by freezing Theta. Lift values are
bounded on one fundamental period only; all resulting periodic sources
and positive-order slow jets extend to every phase.

If L_A,L_B denote the suprema of the absolute right side of (7) for the
two sources on the compact rectangle and one phase period, the zero-mean
periodic antiderivative has sup-operator norm <=2 on period one: integrate
from 0 and subtract its phase mean. Slow differentiation commutes with
this operation. Thus the EXACT transfer gives

\[
                         a_4\le2L_A,\qquad b_4\le2L_B.          \tag{8}
\]

This does not evaluate L_A,L_B. Unlike simply naming unknown primitive
jets again, (2),(5)--(7) specify their missing native derivatives and all
mixed lower terms. For the explicit lift Phi=theta+eta*sin(theta)/2,
source cos(theta), eta=0 and external phase pi/3, (7) gives normalized
fourth coefficient 135/4096 (raw 405/512), while freezing Theta gives zero.
Even Phi40=0 does not eliminate the mixed remainder.
An independent Lagrange-inversion calculation agrees. This is NOT the
parent's different explicit map Theta=varphi+eta*sin(varphi)/2, whose
normalized fourth coefficient is 3/4096 (raw 9/512). An initial new check
incorrectly reused that different control value; it failed and was corrected
before freezing evidence. Historical results for the explicit map are intact.

## 4. Precise remaining gate; no frequency or pressure reselection

Use the parent's EXPLICIT positive coefficients
`alpha=3200000*A*Za`, `beta=1+3200000*A*Zb`, and
`R=N*(K-Kknown)-(1+3200000*A*ZH)*H` with the SAME definitions and N.
Equations (2),(5)--(8) give the sufficient inequality

\[
                             2\alpha L_A+2\beta L_B\le R.       \tag{9}
\]

It is sufficient for that positive compact C4 majorant, NOT necessary for
the actual field. No acceptance is asserted. The next obligation is to
quantify the actual T2 fifth-moment/axis-pressure combination and mixed
fourth-shear inputs, propagate them with the fourth signed-root/scalar and
full moving-phase estimates, and verify (9) at the fixed parameters.
T1 is additionally needed if a fourth full pressure-coordinate jet is
used downstream. Existing all-order smoothness alone does not evaluate
these constants. Finalized first U1 C3, second H2,3/pressure repair,
all physical scales/cutoffs/supports, complete covariance means including
axial velocity/pressure, viscosity and every residual interaction remain
separate open obligations. Historical conditional premises remain conditional.

## Reproduction

Run `python test_audit.py -v`, then
`python audit.py --verify-record evidence.json --output <fresh-path>`.
The audit verifies historical source bindings, exact principal splits,
the symbolic compatible-field limit and every fourth composition partition.
No historical frozen record or tolerance is modified.
