# Native mixed fourth pressure and the remaining slow-jet boundary

Parent: `dca2a2cdbadd9f3e0a56f3a49fc0d4abf6012441`.
Status: `native-pressure-joint-four-and-explicit-fourth-slow-frequency-boundary`.

This checkpoint bounds the complete PREMOD pressure coordinates in joint J4
by A^1547. Radial source identities supply all mixed rows with a radial
derivative; the earlier pure-angular fourth rows retain their full incoming
pressure data. This needs no new fifth mixed field estimate or variance
cutoff assumption beyond the retained native/core hypotheses.

For the FINAL compact fields E,U we derive an exact fourth-order frequency
transfer, conditional on the existing explicit shape gate. One pure-slow
fourth primitive seminorm Gamma remains. A compatible functional control
shows that joint third primitive/source bounds and pure angular fourth
bounds cannot remove it. This is NOT a counterexample to the actual C.9
loop or a proof that the retained profile fails. The next native shear
inputs are explicitly identified as fifth mixed field derivatives.

All parameters, widths C^-20, cutoffs, supports, exact roots, frequency,
pressure gauges and historical evidence are unchanged. The actual shape
gate and Gamma are unverified. No global FINAL leading C4, first C3,
second H2,3, all-scale pulse, complete covariance/means/residual or blow-up
claim follows. Finite checks are not a continuum proof or an independent
foundation review.

Sources: retained [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(4.9), (4.15)--(4.16), (C.11)--(C.16); `higher_jet_budget`,
`actual_fifth_native`, `mixed_shape_cone_bridge`, `fourth_loop_shape_gate`,
`fourth_target_normalization` and `cutoff_higher_jet_interface`.
PDF SHA256: `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The manuscript's final theorem is not a premise.

## 1. Complete native radial pressure equations

Use the same factorial mixed norms Jk in (y=log X,eta),
A=C^4096>=2^16, Bn=A^256, G=exp(A^32), H=exp(G^4),
N=1+floor(H^32). All statements in this section concern PREMOD inputs.
On |eta|<=1, 0<h<=.01, set d=1-eta^2, D=1/2-h, Ap=1/2+h, L=1-2h eta^2,
Hc=D eta+dU, ell=1/2+E_y/E, Pi=Pi_axis+Cp. Exactly,

`p1_y = X*Sq/L-ell*p1`,
`p2_y = X*Sn/(L*E)-(E_y/E)*p2`,                         (1)

where

`Sq=-W*ell-h*(1-2eta U)-Hc*(E_eta/E)`,
`Sn=-W*U_y-Ap*(1-2eta U)*U-Hc*U_eta`
`   -d*Pi_eta+4Ap*eta*Pi+eta*E^2`.

These are identities for the ORIGINAL forward map, not just the reference
profile. Differentiate that map with the five actual density identities
M_y=XU, I_y=X Hx, J_y=XUHx, S_y=X(U^2-E^2/2), Cp_y=E^2/2,
Hx=sqrt(2X)E, to obtain (1). Both squares, axial U and the full axis
pressure survive. The equation for p2 follows from Ns_y+Ns=Sn;
it is NOT p2_y+p2=Sn. No incoming constant is set to zero.

The native ledger bounds E,U,Pi in J4, E^-1 in J4, W,p1,p2 in J3
by Bn. Consequently E_y,E_eta,U_y,U_eta,Pi_eta in J3 are <=4Bn,
Hc<=6Bn and ell<=5Bn^2. Product estimates give

`|Sq|J3<=34Bn^3`, `|Sn|J3<=60Bn^2`.

Also |X|J3<3A and |L^-1|J3<2<A. Equation (1) therefore gives

`|p1_y|J3<=107A^2 Bn^3`,
`|p2_y|J3<=184A^2 Bn^3 <200A^2 Bn^3 <A^771`.           (2)

For r>=1, normalized p_(r,s)=(p_y)_(r-1,s)/r. Summing all these rows
through total order four costs at most (2). The earlier actual fifth-moment
audit gives the PURE angular C4 norm <=A^1546, including all axis and
incoming data. Splitting the angular and radial rows proves

`|p1|J4, |p2|J4 <=A^1546+A^771 <A^1547 <G`.            (3)

Radial differentiation cannot establish the pure angular row or erase its
memory. Bound (3) is not a FINAL modulated-pressure bound.

## 2. Four native shear inputs still use fifth mixed field jets

For each r+s=4 use normalized native coefficients E_(r,s), U_(r,s).
The actual shear definitions a=1-2E_y/E, bs=2U_y/E give

`a_(r,s)=Ra_(r,s)-2(r+1)*E_(r+1,s)/E0`,
`bs_(r,s)=Rb_(r,s)+2(r+1)*U_(r+1,s)/E0`.               (4)

Ra,Rb use only field jets through total four. The r=0, s=4 input was
bounded earlier by the defining slope identities. The remaining field
indices are (2,3), (3,2), (4,1), (5,0), for EACH of E and U. Existing native
J4, pure angular C5 and the E_(1,4),U_(1,4) bounds do not supply these
eight quantitative estimates. This is an input boundary, not evidence
that the actual prescribed slopes cannot bound them. The corresponding
mixed fourth target/loop calculation must retain them and the actual cutoff.

## 3. Exact fourth frequency transfer and a linear missing seminorm

For a zero-mean primitive P(y,eta,varphi) write g=P_varphi and
`P_N=P(y,eta,Ny)/N`. The exact normalized chain rule is

`(P_N)_(r,s)=sum_(c=0..r) N^(c-1)*P_(r-c,s,c)`.        (5)

There is NO extra binomial factor after factorial normalization. In this
section assume the retained sufficient shape condition (Z) of
`fourth_loop_shape_gate`. The parent then gives |P|J3<=H, |g|J3<=H
in all three variables, and the earlier pure angular fourth coefficient
is <=H. Define

`Gamma_P=sum_(r+s=4,r>=1) sup |partial_y^r partial_eta^s P|/(r!s!)`,
`Gamma=Gamma_Acal+Gamma_Bcal`.

For c>=1 in (5), P_(a,b,c)=g_(a,b,c-1)/c. Hence all top rows with at
least one phase derivative are controlled by g J3, although P J3 alone
would not do this. The c=0 slow fourth row is not removed. If b_j is the
sum of suprema of the degree-j jet of P_N, then

`b0<=H/N`, `bj<=N^(j-1)*H` (1<=j<=3),
`b4<=2N^3*H+Gamma_P/N`.                               (6)

For beta=Acal_N, exp(beta0)<2. The normalized fourth exponential jet is

`exp(beta0)*(beta4+beta1*beta3+beta2^2/2`
`                    +beta1^2*beta2/2+beta1^4/24)`

after scalar degree collapse; the same positive polynomial majorizes
the colored mixed jets. Sum degrees zero through four, using N>=H, to get

`|exp(beta)|J4 <=20N^3 H+2Gamma_Acal/N`.                 (7)

The unknown fourth slow jet occurs LINEARLY, not inside exp(Gamma).
It also cannot be deleted by retaining only the large N^3 terms.

## 4. SAME repair, a field-only gate, and no replacement of K

Under (Z) the same leading repair coefficient Ceta4 norm is
<=2000H^16/N. The original width-.1 unit-mass bumps have radial J4 norm
<1e106<Bn, using the established sigma derivatives through five. Their
original e(eta) factor has Ceta4<=Bn. Summing five supports safely bounds
each repair field J4 by

`10000Bn^2 H^16/N <1`.

No radial derivative of a tiny field change is declared small without its
bump-width cost. Apply the native field J4 bound Bn, (6)--(7), and the
same exact repair root to obtain the FINAL COMPACT FIELD estimate

`|Etilde|J4, |Utilde|J4 <=32Bn N^3 H+2Bn*Gamma/N+2`.     (8)

This is for E,U only. It does not bound fourth flux/pressure jets, which
can require additional fifth angular moment data, or all global regions.
The retained K=H^64 was a lower-order envelope and is NOT reselected.
For a separate fourth FIELD budget K^2=H^128, N<2H^32 gives

`32Bn N^3 H+2 <H^128/2`.

Thus the precise sufficient remaining field inequality is

`Gamma <=N*H^128/(4Bn)`.                               (F)

It is NOT evaluated at the actual primitive. A negative value of an
auxiliary gate is not failure of the actual construction. Even satisfying
(F) and (Z) would not prove the full leading C4 or physical-pulse theorem.

## 5. A protected-support countercontrol to deleting Gamma

Fix ANY finite integer N>=1. Choose a fixed smooth chi supported in
(.1,.9), equal to one on [.2,.8], and use the squared-flat radial step,
for example `chi=sigma(10(y-.1))*(1-sigma(10(y-.8)))`.
Set y0=1/4 if N=1, and y0=(floor(N/2)+1/4)/N otherwise. Then chi=1
near y0 and Ny0 is an integer plus 1/4. Define AUXILIARY primitives

`Acal_n=0`, `Bcal_n=n^-3*chi(y)*cos(n*(y-y0))*sin(2pi*varphi)`.

They are zero-mean, vanish in protected collars, and have zero pure
eta-fourth jets. Both Bcal_n and gB_n=partial_varphi Bcal_n have joint J3
norm <=1e211 independently of n,N: the existing fixed-step C3 bound gives
|chi|J3<1e208, and the trigonometric jets cost <1000. This fits the same
H, since log H=G^4>211log10. The fourth derivative after substitution is

`partial_y^4 (Bcal_n(y,eta,Ny)/N)|y0`
` =[n^4+6n^2*(2pi N)^2+(2pi N)^4]/(N*n^3)`.

Its normalized fourth coefficient divided by n tends to 1/(24N)>0.
No bound from those lower-order seminorms alone controls the fourth field
jet at fixed N. A later disjoint exact moment repair cannot change this
jet at y0. This is a continuum non-implication for the functional transfer,
NOT an actual C.9 loop, pressure-compatible global field or actual-profile
counterexample. The source bounds and mixed/native obligations must still
be proved for the selected profile.

## 6. Checks and next obligation

Exact symbolic checks differentiate the ORIGINAL forward pressure map with
the five density flows, extract (4), and use arbitrary independent Taylor
coefficients through degree four to retain every frequency/exponential
partition. Manufactured fixed-N controls independently differentiate the
trigonometric composition. Auxiliary exact-rational checks call the same
field majorant and include an explicitly rejected Gamma budget. Outward
80/110-digit comparisons check universal constants without materializing
the actual huge hierarchy. Historical records remain byte-identical.

Run `python test_audit.py -v`, then `python audit.py --verify-record
evidence.json --output FRESH.json`. Next bound the actual fifth mixed shear
inputs and full fourth slow loop seminorms, retaining the unchanged cutoff
budget. FINAL global jets, first/second coefficient estimates, all-scale
pulses, complete covariance/mean corrections and the full residual remain open.
