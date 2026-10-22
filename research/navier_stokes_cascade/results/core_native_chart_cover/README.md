# Existing-box native coverage and support-aware full velocity norms

Parent: `c69476a089fe3d894179e577b90b5f15ded3c4af`.
The latest checkpoints derive native velocity/primitive jets and certify the
implemented Borel cutoff. They leave actual native-to-band coverage and the
unchanged value/schedule gates. This checkpoint quantifies the ORIGINAL
selected doubled boxes, extends the existing leading VELOCITY rows to their
whole norm hull, and improves primitive bounds using the original axial
support. Pressure/moment source endpoints remain unchanged.

These are analytical implications of the retained regional, analytic-core,
cutoff and exact-repair estimates. The elementary chart bounds are proved
below; finite controls are not continuum proofs of the inherited profile.
No root, profile, frequency, width C^-20, support, pressure gauge, cutoff,
existing first band or native value neighborhood is selected again. No
unconditional actual pulse, complete residual improvement or blow-up follows.

## 1. Actual doubled boxes, not an unrelated Cartesian value domain

The original mesh has halfwidth delta=ell^-6 in (R,Z,T). A selected original
box meets the physical active set. Choose any such intersection point y0,
with native X0 in [Xa,Xb], s0=q/Q in [1/2,2], T0>=0. Every point y in
that box's DOUBLED box has

    |R-R0|, |Z-Z0|, |T-T0|<=3delta.                       (1)

For value bounds here intersect the doubled box with physical T>=0.
Choose fixed outward enclosures 0<rA<=sqrt(Xa), rB>=sqrt(Xb).
The implicit chart obeys

    s-Z^2 s^(2h)=T, eta=Z s^(-D), X=R^2/(2s),
    D=1/2-h, 0<h<1/100, Lprof=1-2h eta^2.                 (2)

On the physical branch, T>=0 implies |eta|<=1. As long as 1/4<=s<=4,
Lprof>=.98 and the exact derivatives satisfy

    |partial_Z log s|=|2eta s^(-D)/Lprof|<5,
    |partial_T log s|=1/(s Lprof)<5.                       (3)

Along the line segment from y0 to y (still physical T>=0), (1)--(3) give
|log(s/s0)|<=30delta. If delta<=1/100, this is <=.3<1/2<log 2.
A first-exit/continuation argument therefore keeps the whole segment in
1/4<s<4: at a putative boundary its logarithmic displacement from
s0 in [1/2,2] would be at least log 2, a contradiction. The implicit
derivative stays positive, so the physical branch continues throughout.
This is a whole-box argument, not a finite implicit-root scan.

If also 6delta<=rA, (1) and R0>=rA give R>=R0/2. Hence

    |log(R/R0)|<=2|R-R0|/R0<=6delta/rA,
    |log(X/X0)|<=delta(12/rA+30)=:w_native.                (4)

For an EXISTING prescribed log-native enlargement of halfwidth mu>0,
the precise sufficient additional gate is

    6 ell_existing^-6<=rA,
    100 ell_existing^-6<=1,
    ell_existing^-6(12/rA+30)<=mu.                        (5)

The code checks the supplied existing floor; it does not return a larger
one or increase mu. All three left sides decrease on every later band.
The second condition is already automatic for ell_existing>=16. The
first and third still require the actual original radius/neighborhood
data. If the original mesh gate is Dbox ell^-6<=dchart, sufficient
compatibility of its UNCHANGED constants is

    6dchart<=Dbox rA, 100dchart<=Dbox,
    dchart(12/rA+30)<=Dbox mu.                             (6)

Failure of (5) or (6) reports a sufficient-certificate gap; it does not
prove failure of the actual implicit chart or justify selecting another band.
All other original selector and q-domain gates remain necessary inputs.

## 2. A concrete norm hull for the same selected boxes

From the first two gates in (5), R0<=2rB and rA<=rB, derive

    Rmin=rA/2 <= R <= (5/2)rB=Rmax,
    1/4<=s<=4,
    Rlo=rA/4 <= Rprof=R/sqrt(s) <=5rB=Rhi,
    Xlo=rA^2/32 <= X <=(25/2)rB^2=Xhi.                   (7)

These are enclosures of the SAME supports, not new physical annuli or
cutoffs. They exactly meet the preceding radial checkpoint's coverage test:
Rlo=Rmin/2, Xhi=2Rmax^2 and Rhi^2=2Xhi.

The Cartesian norm hull (7) is much wider than the actual coupled box image;
it can contain X outside [Xa,Xb]. It is suitable for UPPER JET NORMS only.
Positive growth, nonzero shear and the stress cone are not inferred there.
Their application to the actual physical boxes instead uses (4)--(5) and
the independently justified original value neighborhood and floors.

The proof covers physical T>=0, including one-sided jets at T=0. It does
not analytically continue the Gamma integral to negative argument or
certify a negative-T background extension. A frozen representative outside
the physical covered set requires the original extension to be separately
bounded; it cannot be certified by these one-sided norms. No representatives
are moved to remedy that missing premise.

## 3. Extend the existing leading velocity rows to the norm hull

The old `global_swirl_envelope`, section 3, proves the native logarithmic
J3 bound |E0|J3<H^70 on the compact post-axis stages AND on ALL later
reference radii, not just up to the heat interface Xc. The latter proof
includes the actual normalized terminal factor, heat compensation and the
whole Gamma exterior. For the Gamma factor it uses nonnegative
Z=2(1-eta^2)/X; at eta=+-1 only one-sided limits are taken.
The exact tail estimates in `annular_exterior_foundation`, equations
(6)--(8), also apply for every X>=max(1,Xb). There is no finite-radius
truncation of that velocity argument.

Together with the retained axis-side common envelope Kaxis, the same
ordinary-derivative conversion therefore bounds VELOCITY on all X>=0:

    MF=2(1+Xa^-1) ME, ME=H^70,
    F3=max(Kaxis,MF),
    FX2=max(Kaxis,MF/Xa),
    FXX1=max(Kaxis,3MF/Xa^2).                              (8)

These are factorial C_eta^3,C_eta^2,C_eta^1 norms respectively. Before
the positive annulus use the regular axis-side estimate; only afterwards
use F=E0/sqrt(2X), X>=Xa and
F_XX=X^-2(F_yy-F_y). No division by zero or analyticity across flat joins
is introduced. This extends the stated velocity-row domain in the old
first-stress interface from [0,Xc] to any finite Xhi (indeed all positive
radii). It does NOT extend that interface's positive stress floor tau,
its cone/growth gaps, pressure source integrals or any finite source support.

Likewise `global_axial_source_envelope`, section 3, already bounds

    KU=max(Kaxis,2Ahi,352),
    R1=max(Kaxis,H^2/Xa),
    R2=max(Kaxis,H^3 N/Xa^2)                               (9)

on the whole original leading axial support [0,Xv]. The fixed smooth
zero extension applies to U0 and its jets after Xv. Thus (9) covers
[0,Xhi] for arbitrary Xhi, without replacing H^3 N by a smaller pure
angular norm. The globally finalized positive-order F1,U1 rows from the
same cutoff and full-pressure exact first repair are smoothly compact
on their original common support within Xv. Their old row bounds also
extend to the norm hull by zero. Pi1/Pi2 pressure memory is not inferred
to vanish on I4 merely because these velocities do.

`leading_global_rows` calls the historical formulas directly. Its first
three swirl outputs do not depend on the chosen finite Xhi. This is a
domain consequence of the existing complete velocity proof, not a new
selected-profile measurement or a change to the exact heat exterior.

## 4. A support-aware primitive improvement without lost pressure memory

On the whole hull, every fixed derivative of U_n is supported in the
ORIGINAL common [0,Xv]. Therefore, for any X<=Xhi,

    |partial_eta^k M_n(X)|
       <=min(Xhi,Xv) sup_(0<=Y<=Xv)|partial_eta^k U_n(Y)|.  (10)

For positive radial orders retain M_X=U and its differentiated versions.
At points where any U derivative is nonzero, X<=Xv, so every term in
any derivative of XU can use X<=min(Xhi,Xv) too. Outside that support,
all such derivatives vanish by smooth zero extension. The derivative of
the X coefficient is still ONE, not zero.

Consequently the preceding positive reconstruction recipe remains valid
on the ENTIRE hull with the replacement

    Xhi -> Vsupport=min(Xhi,Xv)                            (11)

ONLY in the XU coefficient bounds and zero-radial-order primitive volume.
Keep the full-hull Rlo,Rhi bounds in E=Rprof F, W=V/Rprof and
J=2eta M/(Lprof Rprof); keep every inverse-radius derivative and
factorial-to-raw angular conversion. `supported_reconstruct` reuses the
old finite polynomial calculator with precisely this NEW support argument.
Its domain remains X<=Xhi even when Vsupport<Xhi. For the actual cover,
Xv>=Xa>=rA^2 implies Vsupport>=Xlo, so the old calculator's positive
annulus guard is also satisfied.

This does NOT assume M=0 on a local U=0 plateau or even at infinity;
(10) is valid with nonzero accumulated primitive memory. It never changes
M's axis datum. The global first five-moment repair still integrates its
complete signed pressure sources Q03,Q13 to the full original Xv,
including the axis-viscosity boundary, not to Vsupport. If Xhi<Xv, future
pressure beyond the chart still affects that repair. The implementation
requires its source endpoint and the support used in (10) to be the SAME
original Xv; shortened pressure data are rejected.

All reconstructed raw velocity and primitive bounds weakly improve over
the old whole-hull volume bound. The E and U entries are unchanged;
W and J improve strictly for positive inputs when Xv<Xhi. This costs no
physical Q power and requires no changed cutoff or repair.

## 5. Join and exact remaining gate

Apply (8)--(11) to the leading tuple through order two and the finalized
first tuple through one. Use the derived band radii (7) and the SAME
implemented Borel cutoff bound K2=4928 in the seven-term Cphi formula.
This removes independent native-radius and hull-domain inputs relative
to (5) and the retained global velocity proofs. The supplied existing
native value-neighborhood gate is recorded separately; its validity does
not manufacture its growth/cone floors.

The unchanged phase/selector requirement remains

    beta Lreq^2>=144 r0^2 Cphi^2,
    beta=lambda_min/[4(1+U^2)^2],                          (12)

along with the rounding/slow gate and every original selector requirement.
The actual infinite Borel schedule and cutoff identity, frozen representative
coverage or extension, leading value/edge hypotheses, tangent propagator
and full supports are not certified merely by (7)--(11). No new first
band, length budget, log neighborhood or physical field is installed.

Checks use exact rational physical chart points, including both eta poles,
the original selector fixture, monotone support-aware majorants and a full
future-pressure negative control. They are manufactured controls, not
selected huge-field sampling or continuum proof of inherited premises.
The audit binds seven historical evidence records and their declared
source hashes. Run `python test_audit.py -v` and `python audit.py
--verify-record evidence.json --output /path/to/fresh.json`.

Next verify (5) and the actual frozen-representative neighborhood/extension,
identify the same infinite Borel prescription with its proved bounds, and
substitute the complete gap and norm data into (12) on the unchanged selector.
Then discharge the homogeneous propagator and remaining support/fixed-map
premises before the first sourced complete covariance and fresh full-pressure/
axial means. Unconditional actual pulse bounds, complete viscous/cutoff/
interaction residual improvement and blow-up remain unproved.
