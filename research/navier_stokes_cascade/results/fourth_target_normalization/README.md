# Tame fourth-angular five-moment target from the same loop primitives

Parent: `4f093b36c10abea2fb4c8e6b31bb522418df7963`.
Status: `tame-fourth-angular-density-moments-and-normalized-target-primitive-gate`.

The preceding checkpoint reduces a global leading axial C4 estimate to the
actual compact primitive B4 and fully normalized repair-target Z4. This
checkpoint ELIMINATES Z4 as an independent unknown: all five raw moment
differences and their normalizations have explicit bounds affine in the
two primitives' fourth angular jets. It then substitutes that estimate into
the previous sufficient compact certificate. The actual moving-phase
fourth primitive jets are NOT evaluated, so the compact/global C4 gate is
still open. No global H2,3, complete physical pulse/covariance/residual or
blow-up conclusion is claimed.

All parameters, roots, widths C^-20, N=1+floor(H^32), schedules, gauges,
supports and historical frozen evidence stay unchanged. The derivation is
conditional on the SAME original compact C4 envelope and accepted joint
C3 loop bounds. Smoothness at fixed parameters is not their quantitative
fourth-order bound, and finite manufactured checks are not continuum proof.

## 1. Use value/lower-jet smallness, not an exponential of the unknown C4 norm

Use factorial angular norms |.|k on [-1,1], including endpoints. Keep the
actual compact geometry A^-1<X<A, radial length <A, and the ORIGINAL
native |E|4,|U|4<=Ahi, with A>=1 and Ahi=A^256>=1. The accepted joint
primitive J3 norm is <=H. Define the actual highest jets

`a4=sup_(y,eta,phi) |partial_eta^4 Acal|/4!`,
`b4=sup_(y,eta,phi) |partial_eta^4 Bcal|/4!`.

These are partial derivatives at fixed AUXILIARY phase of the actual
pulled-back primitives. All slow derivatives of the moving inverse phase
inside Acal,Bcal belong in a4,b4. Only the subsequent substitution phi=Ny
is eta-independent; freezing the internal inverse phase is incorrect.

Put a=Acal(y,eta,Ny)/N, b=Bcal(y,eta,Ny)/N and t=H/N<=1. The selected
N>H^32 already ensures this; no frequency is reselected. We have
|a|3,|b|3<=t, fourth normalized jets <=a4/N,b4/N, and

`E_N=E exp(a)`, `U_N=U+b`.

For normalized Taylor jets a_j, the fourth coefficient of exp(a)/exp(a0) is

`a4jet+a1*a3+a2^2/2+a1^2*a2/2+a1^4/24`.

The known lower-jet sum costs at most (49/24)t^2<3t^2. Since
exp(|a|3)<=3, the C3 exponential difference is <=3t, while its fourth
coefficient is <=3a4/N+9t^2. This proves

`|exp(a)-1|4 <= (12H+3a4)/N`.

For exp(2a), its fourth coefficient contains
`2a4jet+4a1a3+2a2^2+4a1^2a2+(2/3)a1^4`.
Using exp(2t)<9 gives C3 difference <=18t and lower fourth cost <=96t^2.
Hence `|exp(2a)-1|4 <=(114H+18a4)/N`.
This is TAME: an arbitrarily large a4/N is allowed. It is not put into an
exponential, assumed small, or used to justify a new root/frequency choice.

For b^2 the C3 norm is <=t^2 and its fourth jet is at most
`2t b4/N+3t^2`; thus `|b^2|4<=(4H+2b4)/N`.
For b exp(a), the fourth Leibniz terms give
`3(a4+b4)/N+23t^2`, and its C3 norm is <=3t. Therefore
`|b exp(a)|4<=(26H+3a4+3b4)/N`.

## 2. Every leading modulation density, with both squares retained

The following are coefficients D such that the COMPLETE density difference
has C4 norm <=D/N. They use the actual E_N,U_N, not a linearization.
The entries are affine in H,a4,b4; Ahi is a fixed original envelope.

| Density difference | D |
|---|---|
| U_N-U | H+b4 |
| E_N-E | Ahi(12H+3a4) |
| U_N^2-U^2 | Ahi(6H+4b4) |
| E_N^2-E^2 | Ahi^2(114H+18a4) |
| U_N E_N-U E | Ahi^2(38H+6a4+3b4) |

For U^2 keep 2Ub+b^2, using the preceding tame b^2 estimate. For UE use
`UE(exp(a)-1)+E b exp(a)`. Neither square is removed, and the two different
exponentials in the E and E^2 rows are not confused. All high jets remain
linear because only lower jets multiply a fourth derivative. Applying a
crude C4 product bound to b^2 would unnecessarily introduce b4^2; doing so
for the swirl exponential would hide the available value/lower-jet control.

These are the five densities for the LEADING modulation/first exact
five-moment repair. They do not replace the higher-coefficient pressure
source w1 or its shifted viscosity, nor the F1^2 term of the full P2 source.
All those terms remain in the historical higher-pressure/second-source
checkpoints and must still be used for the later second-coefficient repair.

## 3. Integrate the whole incoming difference before normalizing

Let D_U,D_E,D_U2,D_E2,D_UE denote the table entries. The original moments
are M=int U, I=int sqrt(2X)E, J=int sqrt(2X)UE,
S=int(U^2-E^2/2), Cp=int E^2/(2X). The compact bounds on X, inverse X
and radial length give N times their C4 differences bounded by

`mM=A D_U`, `mJ=2A^2 D_UE`, `mI=2A^2 D_E`,
`mS=A(D_U2+D_E2/2)`, `mCp=A^2 D_E2/2`.                 (1)

All original forward moment data before modulation are the SAME in the two
profiles. Only their DIFFERENCES cancel those identical incoming data.
Neither actual incoming moments nor the eta-dependent pressure axis datum
are assigned zero. Full pressure is Pi_axis+Cp, so deltaPi=deltaCp with
the SAME axis datum, including its fourth derivative. A changed gauge would
give an additional term. The differences in (1) persist to repair entry;
they are not erased where the local modulation density has stopped.

At the unchanged first reserved repair interval let its left scale be rhop
and e(eta)=e_star/(1+eta^2). Supply
`s_rho >=1/(sqrt(2)*rhop^(3/2))`, with positive rhop,e_star and
0<lambda<=.01. The following are exact C4 normalizer bounds (or the stated
square-root enclosure), retaining ALL normalizer derivatives:

| Row | Ordinary dividing factor | C4 norm of reciprocal, bounded by |
|---|---|---|
| M | rhop e | wM=5/(rhop e_star) |
| J | sqrt(2) rhop^(3/2) e^2 | wJ=25 s_rho/e_star^2 |
| I | sqrt(2) rhop^(3/2) e | wI=5 s_rho/e_star |
| S | rhop e^2 | wS=25/(rhop e_star^2) |
| Cp | e^2 | wCp=25/e_star^2 |

Indeed |1+eta^2|4=5 and |(1+eta^2)^2|4=25. The C3 value 24 is not valid
for the latter C4 normalizer. The executable guard checks that s_rho is a
genuine upper enclosure, not an arbitrary positive number.

Only AFTER normalization transform to `(m,(m-j)/lambda,i,-s,c)`.
Thus each fourth normalized target jet is bounded by its full C4 row:

\[
 Z_4\le\max\{w_Mm_M,(w_Mm_M+w_Jm_J)/\lambda,
                         w_Im_I,w_Sm_S,w_{Cp}m_{Cp}\}.          \tag{2}
\]

Here Z4=N sup_eta ||partial_eta^4 z_N||inf/4!, exactly the input of the
parent checkpoint. The sign of the repair target is immaterial for this
bound. The inverse-lambda loss in the M-J row remains, even though its
exact limiting matrix is well conditioned. No cancellation between the
two differently normalized incoming rows is invented.

Each row in (2) has nonnegative affine coefficients of H,a4,b4. Taking
the componentwise maxima of those coefficients gives explicit rational
constants ZH,Za,Zb, implemented by `target_affines`, such that

\[
                              Z_4\le Z_H H+Z_a a_4+Z_b b_4.   \tag{3}
\]

Unlike a new unspecified target supremum, these constants are calculated
from A,Ahi,rhop,s_rho,e_star,lambda by (1),(2). No actual high jet is
silently replaced by the old H or H16 budget. This supplies the complete
target-bound reduction, not the remaining moving-phase primitive estimates.

## 4. Join to the SAME accepted compact repair and its global axial interface

The parent's actual B4 (full C4 primitive norm) is <=H+b4. Keep its SAME
low-order root cap r=2000H^16/N and known part

`Kknown=Ahi+6400Ar+960000000000Ar^2`.

Substituting (3) into its sufficient compact C4 majorant leaves only TWO
actual high-jet inputs:

\[
 \boxed{(3200000 A Z_a)a_4+(1+3200000 A Z_b)b_4
 \le N(K-K_{known})-(1+3200000 A Z_H)H.}                       \tag{4}
\]

This is equivalent to acceptance of the SUBSTITUTED positive majorant for
trial compact bound K. It is sufficient for the parent's actual-data
certificate, not necessary for the true field and not an assertion that
the selected data fail. Negative headroom blocks this particular majorant;
it does not prove nonexistence of a sharper estimate or of the construction.
No C, N, root, tilt or schedule is changed to make it pass.

The same exact five-row repair restores all moments as ETA FUNCTIONS, so
its future pressure memory is preserved. It does not reset the separate
positive-order Pi1/Pi2 constants on I4. If (4) is verified with actual a4,b4,
the parent's axis/compact/reference C4 cover and full-support signed w1
moment interface apply, with all viscosities and the true Xv unchanged.
Actual finalized first-coefficient angular bounds, coefficient cutoff jets
and full second-repair data are still required for H2,3.

Next quantify the actual fourth slow derivatives of both primitives,
including the variance root, moving inverse phase and periodic integration,
then verify (4) and the retained-scale inequalities. Full physical pulse
derivatives/supports, complete-wave covariance with axial velocity and
pressure means, complete viscous/cutoff/interaction residual estimates and
independent inherited foundations remain open. The result is not an
unconditional global profile or blow-up theorem.

## Checks

Exact algebra checks all fourth exponential partitions, normalizer factors,
inverse-lambda row, affine row maximum and substituted-gate identity. A
manufactured positive moving lift has nonzero fourth pulled-back derivative
9/512 where freezing the lift gives zero. Fifteen manufactured samples
independently differentiate all five density and target rows, including
eta=+/-1, and compare to normalized Leibniz convolutions. They are not
actual enormous loops, PDE solutions or interval certificates. The
analytical argument supplies the bounds, not passage of a sample grid.

```
python test_audit.py -v
python audit.py --verify-record evidence.json --output fresh.json
```

Fresh output only. Historical evidence is byte-bound and unchanged.
