# Actual pre-modulation mixed shears and the variance-cutoff shape gate

Parent: `04a9058da056821a2f4833c10550ec958a9424c1`.
Status: `actual-premodulation-mixed-inputs-and-retained-variance-cutoff-budget-gap`.

This checkpoint DISCHARGES the parent's mixed inputs YE,YU for the actual
PRE-MODULATION compact profile, relative to the retained reference C4,
earlier annular-root and native-field estimates. It also bounds the native
a,bs,ts,vs in angular C4. This is not a global J5 estimate and does not
bound the finalized modulated field. Pure fifth field jets and entering
fifth moment data remain separate inputs.

A source audit finds an additional quantitative premise: the retained
variance cutoff is specified by smoothness, range and plateau locations,
not by a shape or derivative cap. Its fourth-order bound cannot be borrowed
from the explicitly chosen radial squared-flat step. We give the precise
shape-dependent variance-target estimate and a compressed-cutoff obstruction.
Earlier loop H estimates must retain the corresponding cutoff-budget premise
at their own orders too. Their frozen evidence is preserved, not rewritten.
Passing those tests does not supply the missing continuum premise.

No root, cutoff, width C^-20, support, pressure datum, N=1+floor(H^32) or
other construction parameter is changed. This remains conditional on the
inherited regional/core estimates; it is not independent foundation review,
global H2,3, unconditional pulse/covariance/residual closure or blow-up.

## 1. Use a prescribed first radial derivative, not a missing J5 norm

Distinguish the attachment parameter B_old from the native derivative
envelope Bn=Ahi=A^256, A=C^4096. On the same compact rectangle
`[Xa exp(tc/2),rhop exp(5)] x [-1,1]`, higher_jet_budget gives
`|E|Ceta4,|U|Ceta4<=Bn` and the original value/floor estimates. The
following bounds refer to that profile BEFORE compact modulation.

On the actual prescribed continuation through Xi=110,

`(log E)_y=(1-a)/2`, `U_y=-kappa*beta*X*ns,r/2`.

Here a is kappa*p1,r followed by its convex taper to 4/5. All radial
cutoffs are eta-independent. The retained reference estimates give

`|p1,r|Ceta4<=1e18 X Lambda K_old^12<B_old/100`,
`|ns,r|Ceta4<=1e8 K_old^8`.

Thus `|(log E)_y|Ceta4<B_old`, `|U_y|Ceta4<55e8 K_old^8<B_old`.
The last comparisons are checked at the ACTUAL retained Md64 parameters
in logarithmic interval arithmetic at 80 and 110 digits. There is no
inverse kappa and no differentiated kappa here: one radial derivative is
the DEFINING slope. Higher radial orders would still cost C^20 per narrow
cutoff derivative. The radially flattened reference's physical a_r=0 is
NOT substituted for its integrated p1,r in this actual prescription.

On the shape segment, with z=log(X/Xi),

`log(CE)=z/10+(1-sigma(z/Tsh))*ell_i+sigma(z/Tsh)*log f`,
`Tsh=2000 B_old`, `f=(1+eta^2)^-1`, `U=Gi`.

Use the retained |ell_i|C4<B_old, |log f|C4<=682 and sigma'<=64.
Since B_old>=682,

`|(log E)_y|C4<=1/10+64(B_old+682)/(2000 B_old)<1`, `U_y=0`.

The following ideal power has log slope 1/10. During the width-one axial
restoration, `U=4eta+g(eta)*(1-sigma(y+8))` with the SAME
`|g|C4<=3e-18`; consequently |U_y|C4<=64*3e-18<1. This uses the
retained C4 transfer, not a relabelled C1 estimate.

## 2. The earlier annular bumps are relative E edits, of width .08

On -6<log(X/XR)<-5 the actual earlier repair has

`E=E0*(1+h)`, `h=sum_(j=1..3) e_j(eta)*beta_j(y)`,
`U=4eta+sum_(j=1..2) u_j(eta)*beta_j(y)`,
`E0=P*f*exp(y/10)`, `|e_j|C4,|u_j|C4<=d4=5.72e-14`.

The bumps HERE are `beta=sigma'((t-c)/.08+.5)`, without a unit-mass
1/.08 prefactor. They are NOT the later width-.1 normalized leading
repair bumps or the positive-order I3 coefficient bumps. Their retained
bounds are beta<=64 and |beta_y|<=250000, including the inverse width.
Summing three E/two U coefficients is safe for the WHOLE angular norm;
disjoint radial supports do not justify replacing a sum of global suprema
of different angular orders by a single selected component norm.

Put h0=3*d4*64 and h1=3*d4*250000. Then

\[
 |h|_{C^4}\le h_0<1/2,\quad |h_y|_{C^4}\le h_1,
\quad |E_y|_{C^4}\le5P[(1+h_0)/10+h_1]<P,\quad
 |U_y|_{C^4}\le2d_4(250000)<1.                       \tag{1}
\]

Indeed the two poles of f give |f^(j)|/j!<=1 on real eta, hence |f|C4<=5,
and exp(y/10)<=1 here. Differentiating the RELATIVE edit gives
`E_y=E0*((1+h)/10+h_y)`. Its physical shear satisfies the exact identity

`a=4/5-2h_y/(1+h)`,

so |a|C4<=4/5+2h1/(1-h0)<1. No derivative of the annular root beyond
the already bounded fourth order is needed for this ONE radial derivative.
We do not replace the eta-dependent earlier quadratic map by the later
constant-Q map or change the accepted root.

## 3. Cover the whole remaining compact prefix and join uniform norms

Before the first reserved repair ends, every remaining swirl stage has
the form `log E=s(y)+theta(y)*log f`, including the constant log(2) term
in the interpolation. The exact schedule gives |slope_base|<=1,
|theta_y|<=64/Tf, Tf>=200, and |log(2/(1+eta^2))|C4<=683. Therefore
`|(log E)_y|C4<=1+(64/200)*683<1000`.
This includes eta-independent radial ramps and the terminal radial factor;
it does not evaluate the later heat endpoint at negative argument.
For `U=k(y)eta`, the prescribed taper has |k'|<=4, giving |U_y|C4<=8.

Use the following FIVE-region cover. The ideal/restoration row includes
the intervening ideal gaps; the actual parameter comparison Xsep/XR<exp(-8)
ensures it precedes the earlier repair. The later prefix ends at the
same first reserved-patch endpoint; no support is moved.

| Region | bound for |(log E)_y|C4 | bound for |U_y|C4 |
|---|---:|---:|
| Actual prescribed continuation | B_old | B_old |
| Shape | 1 | 0 |
| Ideal and axial restoration | 1/10 | 1 |
| Earlier annular repair | 1 | 1 |
| Remaining compact prefix | 1000 | 8 |

SUM the regional bounds, not merely their maxima, to control the norm
which sums the suprema of all angular derivatives over the whole rectangle:

`|(log E)_y|C4<=B_old+1003<C`, `|U_y|C4<=B_old+10<C`,
`|a|C4<=2B_old+2007<C`.

The inequalities follow from the SAME `log C=10000B_old+100log(1/epsilon)`.
Together with the global native |E|C4<=Bn and C<=Bn they give the actual
mixed-input estimates

\[
 YE:=\sup {|E_{y\eta^4}|\over4!}\le B_n^2,\qquad
 YU:=\sup {|U_{y\eta^4}|\over4!}\le B_n.             \tag{2}
\]

This discharges precisely two of the parent's high inputs. It does not
bound E_eta5 or U_eta5, entering fifth moment data, or the FINAL modulated
E_(y eta4), whose loop/pressure derivatives must still be included.

Using the inherited a>=1/A and native E^-1 value floor, finite reciprocal
Taylor algebra gives |a^-1|C4<=5A^9<=A^10 and
`|bs|C4=|2U_y/E|C4<=10A^5 Bn^5<=A^1286`.
Consequently

`|a|C4<=A`, `|ts|C4<=A^1296`, `|vs|C4<=A^2583`,

where ts=-bs/a, vs=a+bs^2/a. Both axial-square terms and all reciprocal
derivatives are retained. The harmless constants are absorbed using
A>=2^16; A,G,H,N themselves are never materialized.

## 4. Correction: a variance cutoff window is not its derivative budget

The retained manuscript, printed page 160, (C.9), gives the range and the
two plateau thresholds of zeta_L but no formula or norm cap. The retained
PDF was read and its SHA256 independently checked as
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
This agrees with the repository's historical manuscript binding.

Let w=delta_L/8 and let
`K(t)=zeta_L(2+delta_L/8+w*t)`, 0<=t<=1. Define the ACTUAL shape input

`Z=max(1,max_(j=1..4) sup |K^(j)|)`.

It is finite for a fixed smooth cutoff. But support, range and even
monotonicity do not give a parameter-only upper bound. For the auxiliary
family `K_n(t)=1-sigma(n*(t-1/2)+1/2)`, all permitted plateau constraints
hold. Since sigma'(1/2)=8 and sigma'''(1/2)=-3712, their first shape norm
is at least 8n. The third derivative has magnitude 3712n^3 at the center
and vanishes at distance 1/(2n); the mean value theorem gives

`sup |K_n''''|>=7424 n^4`.

This is a counterexample to an automatic uniform derivative budget, not
a replacement for the retained cutoff or proof that it fails a gate.
In particular the explicit RADIAL sigma bounds cannot be assigned to
this separate variance cutoff merely because both are smooth transitions.

The first/third-order loop estimates in loop_modulation_bounds and
higher_jet_budget used quantitative derivatives of the variance target.
Their existing parameter-only H derivations therefore require actual shape
enclosures, or a separate justified estimate of the actual composition that
avoids them. All-order smooth dependence still follows for a
fixed smooth cutoff; the claimed FIXED numerical jet envelopes do not
follow from the manuscript's qualitative cutoff clause alone. Future use
of those envelopes must retain this additional premise. Historical sources,
records, tolerances and roots are not edited to hide this distinction.

## 5. The explicit fourth variance-target shape gate at the SAME parameters

For `s=zeta_L(vs)*sqrt((2+delta_L/2-vs)/a)`, work on the cutoff support.
Here `delta_L/4<=w0:=2+delta_L/2-vs<=3`, a>=1/A, and V=A^2583 bounds
|vs|C4. Fixed-order scalar composition obeys
`|g(f)|C4<=3 M max(1,|f|C4)^4` if M bounds raw derivatives 0..4.
Applying this finite Taylor inequality, not an exponential of a jet norm,
gives

`|zeta_L(vs)|C4<=2e4 Z delta_L^-4 V^4`,
`|sqrt(w0)|C4<=4e5 delta_L^-4 V^4`,
`|a^(-1/2)|C4<=21 A^9`.

For the middle bound |w0|C4<=4V and raw derivatives of sqrt on the
stated positive range are <=512 delta_L^-4. For the last bound the
raw derivatives through four are <=7A^5. Flatness extends the product
by zero; no square root of a negative radicand is used outside support.
Thus

\[
              |s|_{C^4}\le10^{12}Z\delta_L^{-8}A^{20673}.       \tag{3}
\]

An explicit sufficient condition for |s|C4<=G=exp(A^32) is

\[
 \log Z+12\log10+20673\log A+8A^{20}\le A^{32}.                \tag{4}
\]

For EXAMPLE, a genuine bound Z<=1e100 would suffice at every A>=2^16:
`120log10+30000log A+8A^20<A^32`. Dividing by A^20, using log A<=A,
and checking the decreasing/increasing comparison at 2^16 proves it.
The interval code checks that comparison without forming G. A canonical
squared-flat specialization could provide such a budget, but this run
does NOT identify it with or substitute it for the retained cutoff.
Actual Z and (4) are NOT verified. Even accepting (4) would not supply
the missing pure fifth field/pressure inputs or the complete compact gate.

## Reproduction and next obligation

Run `python test_audit.py -v`, then
`python audit.py --verify-record evidence.json --output <fresh-path>`.
Rational annular bounds, symbolic prescribed/relative-bump identities,
two-precision actual parameter comparisons and the cutoff countercontrol
are separate from numerical diagnostics. The auxiliary logarithmic gate
can accept or reject supplied shapes; it is not an actual cutoff enclosure.

Next: identify or bound the unchanged variance-cutoff shape at all needed
orders; quantify actual pure fifth field and entering fifth moment jets
(the retained axis fifth datum already has its earlier bound); propagate
the complete signed-root/moving-phase estimates into the unchanged compact
inequality. Finalized global C4/H2,3, cutoff/scale acceptance, all physical
pulse derivatives/supports, complete covariance means and the viscous,
cutoff and interaction residual remain unproved.
