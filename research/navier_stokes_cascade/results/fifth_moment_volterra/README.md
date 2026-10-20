# Fifth forward moments: exact Volterra kernels and a quantitative input gate

Parent: `e85d339f19d696a915a047a923cb0303f4824acd`.
Status: `fifth-forward-moment-volterra-kernels-and-fourth-input-majorants`.

The parent identified the exact fifth-moment combinations needed by fourth
pressure-coordinate jets. Here these cumulative unknowns are ELIMINATED into
pure fifth field jets and the original incoming data. Combining the energy
and pressure rows before taking absolute values gives a local swirl kernel
which vanishes at the current radius and has a quadratic length bound.
We also quantify the parent's lower-jet remainders and isolate the two mixed
shear inputs. None of the remaining actual fifth field/mixed inputs is
silently assigned its old C4 bound. This is a conditional estimate, not
acceptance of the compact gate or a failure of the actual profile.

Keep all roots, supports, gauges, radial widths C^-20 and the selected
N=1+floor(H^32) unchanged. Historical evidence is frozen. Global leading
C4, finalized first C3/second H2,3, full physical pulses, complete covariance
means and the viscous/cutoff/interaction residual remain open. Tests check
exact algebra and manufactured fields, not the continuum construction.

## 1. Differentiate genuine forward moments, including the entering traces

Fix an eta-independent left endpoint a>0 and a current X>=a on the SAME
compact pre-modulation rectangle. All jth angular jets below are normalized
by j!. Write e5=E_eta5/120, u5=U_eta5/120 and

`Lu=sum_(j=1..4) U_j U_(5-j)`,
`Le=sum_(j=1..4) E_j E_(5-j)`,
`Lue=sum_(j=1..4) U_j E_(5-j)`.

The exact fifth forward identities are

\[
\begin{split}
 M_5&=M_{in,5}+\int_a^X u_5\,dx,\\
 I_5&=I_{in,5}+\int_a^X\sqrt{2x}\,e_5\,dx,\\
 J_5&=J_{in,5}+\int_a^X\sqrt{2x}(Ue_5+Eu_5+Lue)\,dx,\\
 S_5&=S_{in,5}+\int_a^X(2Uu_5+Lu-Ee_5-Le/2)\,dx,\\
 Cp_5&=Cp_{in,5}+\int_a^X(Ee_5+Le/2)/x\,dx.
\end{split}                                                     \tag{1}
\]

These follow from normalized Leibniz convolution, not a linearization.
Both squares are retained. Smoothness on this closed compact rectangle
justifies differentiation under finite integrals. The incoming functions
are those of the SAME profile; no entry value or derivative is reset.

## 2. Combine the pressure and energy memories BEFORE estimating them

Let d=1-eta^2, D=1/2-h, Hx=sqrt(2X)E(X), and k5=(Pi_axis)_5.
The parent's brackets are

`T1=d M5-(D eta I5+d J5)/Hx`,
`T2=(d U(X)-D eta)M5-d S5-Xd(Cp5+k5)`.

Substituting (1) gives the exact Volterra formulas

\[
\begin{split}
 T_1={}&T_{1,in}+\int_a^X d\left[1-{\sqrt{2x}E(x)\over H_x}\right]u_5\,dx\\
 &-\int_a^X{\sqrt{2x}(D\eta+dU(x))\over H_x}e_5\,dx
             -{d\over H_x}\int_a^X\sqrt{2x}Lue\,dx,\\
 T_2={}&T_{2,in}+\int_a^X[d(U(X)-2U(x))-D\eta]u_5\,dx\\
 &+d\int_a^X(1-X/x)E(x)e_5\,dx
       -d\int_a^X Lu\,dx+{d\over2}\int_a^X(1-X/x)Le\,dx,
\end{split}                                                     \tag{2}
\]

where

`T1,in=d Min,5-(D eta Iin,5+d Jin,5)/Hx`,
`T2,in=(d U(X)-D eta)Min,5-d Sin,5-Xd(Cpin,5+k5)`.

The factor 1-X/x comes from the COMPLETE Cp row together with the E-square
in S. Dropping Cp would replace it by 1, giving a wrong leading coefficient.
The axial factor U(X)-2U(x) similarly retains -WU and the U-square in S.
Neither is generically zero. At eta=+-1, T2=-D eta M5 still survives.
At X=a the integrals vanish but the entire entering memory remains.

For ell=X-a, the exact nonnegative kernel mass is

\[
 \kappa(X,a)=\int_a^X(X/x-1)\,dx=X\log(X/a)-\ell
       \le {\ell^2\over2a}=:\widehat\kappa.            \tag{3}
\]

Indeed (X-x)/x<=(X-x)/a. This is a LOCAL quadratic length gain for the
swirl part of the fifth combination, not removal of accumulated pressure
or a uniform smallness claim on the long compact interval. It does not
justify dividing by d or extending compact X floors to the axis.

## 3. Explicit positive bounds, with only two pure fifth field inputs

Use the inherited A^-1<=a<=X<=A, A>=1, native |E|4,|U|4<=B=Ahi,
the value bounds |E|,|U|,E^-1,Hx^-1<=A and all five moment/axis C4 norms <=B.
The Hx reciprocal is the inherited original coordinate value bound, not
inferred from E^-1<=A alone; geometry and that latter bound would only
give a larger power of A.
Keep the original lower C3 pressure/shear envelopes <=B. Put

`FE=sup |e5|`, `FU=sup |u5|` on [a,X] x [-1,1].

These actual higher inputs remain unquantified. The accepted C4 product
algebra bounds |Lu|,|Le|,|Lue|<=B^2. Let nonnegative
`m,i,j,s,c,k` bound the entering NORMALIZED fifth jets of the five
moments and axis datum. Then, at each eta,

\[
\begin{split}
 B_{in,1}&=dm+A(D|\eta|i+dj),\\
 B_{in,2}&=(dA+D|\eta|)m+ds+Xd(c+k),\\
 |T_1|\le B_1&:=B_{in,1}+\ell d(1+2A^3)FU
       +2A^2\ell(D|\eta|+dA)FE+2dA^2\ell B^2,\\
 |T_2|\le B_2&:=B_{in,2}+\ell(3dA+D|\eta|)FU
       +dA\widehat\kappa FE+dB^2(\ell+\widehat\kappa/2).
\end{split}                                                     \tag{4}
\]

Here sqrt(2x)<=2A and Hx^-1<=A suffice for B1. B2 keeps the pressure
kernel gain instead of separately paying for S5 and Cp5. All highest field
and entry jets occur LINEARLY. Taking sup_eta of (4), or conservatively
replacing both d and |eta| by one in these nonnegative majorants, gives
uniform bounds. Setting these incoming bounds to zero without matching
identities is invalid, even on an interval where U=0 or modulation stopped.

### The parent's lower pressure remainders now have a quantitative bound

The parent defines R1,R2 by setting ONLY fifth moment/datum coefficients
to zero in the fourth Taylor coefficient of the original forward map.
They are formal jet remainders, not a different physical profile. From
|eta|4<=2, |d|4<=4 and the factorial product algebra,

`|L^-1|4<=1/(1-2h*4)<=25/23<2`,
`|E^-1|4<=A sum_(j=0..4)(AB)^j<=5 A^5 B^4`,
`|Hx^-1|4<=5 A^6 B^4`.

Truncating a fifth coefficient makes the derivative-jet norm <=4B.
The Pi=Cp+axis jet norm is <=2B and its derivative <=8B. Applying the
original forward formula WITHOUT freezing any denominator gives the safe
ledger

| Formal truncated quantity | C4 jet bound |
|---|---:|
| W | 21 A B |
| numerator in Q's moment fraction | 27 B |
| Q | 156 A^7 B^5 |
| Ns | 91 A B^2 |
| R1 | 2496 A^8 B^5 |
| R2 | 7280 A^7 B^6 |

For example Ns costs 21AB^2 from -WU, 26AB from its moment numerator
and 44B from both pressure terms. Using the looser |L^-1|4<=16 in the
last multiplication yields the displayed constants. Thus both remainders
are <=R0:=10000 A^8 B^6. Combining with the parent's EXACT factor-5
principal rows gives

\[
 |p_1|_4\le B+R_0+10\sup_\eta B_1,\qquad
 |p_2|_4\le B+R_0+10A\sup_\eta B_2.                 \tag{5}
\]

The extra B is the already bounded lower C3 norm. Equations (4)--(5)
eliminate unknown cumulative fifth moments over the compact rectangle,
but not the fixed entering fifth data or actual pure fifth field norms.

## 4. Mixed shear derivatives: charge only the missing highest jets

Write YE=sup|E_(y eta4)|/24 and YU=sup|U_(y eta4)|/24.
The original J4 budget supplies |Ey|C3,|Uy|C3<=B, NOT YE,YU.
The exact fourth coefficients of a=1-2Ey/E and bs=2Uy/E split as

`a_4=-2 Ey_4/E0-2 sum_(j=0..3) Ey_j (E^-1)_(4-j)`,
`bs_4=2 Uy_4/E0+2 sum_(j=0..3) Uy_j (E^-1)_(4-j)`.

Consequently

\[
 |a_4|\le2A\,YE+10A^5B^5,\qquad
 |(b_s)_4|\le2A\,YU+10A^5B^5.                      \tag{6}
\]

These are SHEAR jets, not the parent's similarly named primitive a4,b4.
Prescribed-slope formulas can supply YE,YU on individual regions; no
full J5 or radial-frequency gain has been invented here.

## 5. Remaining actual obligation and reproduction

The next quantitative native inputs are now precisely FE,FU,YE,YU and
the six entering fifth jets. Their actual same-root/region estimates must
be propagated from the native core through the prescribed continuation,
shape and eta-dependent earlier annular repair. For the loop itself only
p2, not p1 as a new varying variance argument, enters the signed-root
equation. Apply the parent's fourth scalar/implicit/moving-phase formulas,
including q derivatives through five and zero-variance limits. Then check
its unchanged sufficient gate 2 alpha LA+2 beta LB<=R. No such acceptance
is asserted by (5) or by all-order smoothness alone.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output <fresh-path>`.
Twelve manufactured cases compare BOTH identities (2) with independently
differentiated genuine five-moment integrals, at three radii and four
angles including endpoints, with nonzero incoming data and axis pressure.
Exact original-forward remainder and mixed-shear checks accompany them.
None is a sample, bound or PDE solution of the enormous selected profile.
