# Q97: changed complex wave covariance and global pressure moments

Parent: `103b7e3a2433301edac90f0d04c128b6653ca9e9` (Q96).
This checkpoint forms the first sourced covariance change and its exact
radial/global pressure-moment identities WHILE the incoming final mean
velocities are held fixed. It does not perform the subsequent fresh mean
updates. Bounds remain conditional on Q96's actual native/source/inverse/
support contracts. No frozen evidence, original root, field, width or
accepted domain is changed. The general complex covariance and moment
identities were already established in core_sourced_covariance_moment_join.
Q97 replays that calculus with the CURRENT Q96 exponential-polynomial
response, restoration/evaluation costs and derived poles; these identities
are not presented as new discoveries.

## 1. Complex carriers: all six tensor entries and both orientations

Use real Fourier fields with a_-m=conjugate(a_m). The old first-real wave
has carriers +/-1; Q94's nonzero source and Q96's response have +/-1,+/-2.
Let w be the COMPLETE old cut wave and v the COMPLETE sourced addition.
Exactly

    Delta W=<w tensor v+v tensor w+v tensor v>_theta,
    W_changed=W_old+Delta W.                            (1)

For Fourier coefficients (not cosine amplitudes), averaging gives the sum
of a_m tensor b_-m, b_m tensor a_-m and b_m tensor b_-m over the appropriate
NONZERO input modes. Different output frequencies do not survive angular
averaging. This is physical angular mean, not auxiliary Haar mean.

Equivalently, write the old first carrier A cos(theta)+B sin(theta), and
the new first carrier C cos(theta)+D sin(theta), plus the new second carrier
E cos(2theta)+F sin(2theta). The exact change is

    [A tensor C+C tensor A+B tensor D+D tensor B
       +C tensor C+D tensor D+E tensor E+F tensor F]/2. (2)

Here C,D may BOTH contain a sourced principal component. The Q37 real-pair
formula must not be reused with D treated as only a small curl remainder.
No conjugate, second carrier, selected off-diagonal entry, old/new
orientation or new square is omitted. Formula (2) preserves all SIX
symmetric tensor entries. It asserts no cancellation or sign of a flux.
The old covariance is retained, not replaced by its increment.

## 2. Raw-jet bounds and the enlarged covariance domain

At fixed maximum order N, Q93 supplies old raw jet bounds
A_j=j!*Uold_j. Q96 supplies evaluated, restored and completely cut raw
response bounds V_(m,j) for positive carriers m=1,2. Their negative carriers
are conjugates by uniqueness of the real-background ODE. Thus a sufficient
new Fourier coefficient mass is V_j=2[V_(1,j)+V_(2,j)]. Use the SAME
incoming state, radial hull and label Gaussian P_i for both solves.
No extra copy or label averaging is inserted.

For raw mixed-jet convolution *, the complete covariance change has bound

    D_j=2(A*V)_j+(V*V)_j.                              (3)

The old and new wave orders are 1/2 and B=1-kappa. Crosses have order
C=B+1/2=3/2-kappa and the square has order 2B. Since 2B-C=1/2-kappa>0,
the new square embeds into order C for epsilon<=1. It is still retained
in (3); it is not subtracted from the coefficient polynomial.

Choose q as the maximum OLD and NEW wave pole at ALL requested jets.
The covariance pole is 2q. Its full pointwise weight is zeta P_i^2
(delta^-2q), from two half-weight waves. For subsequent mean/moment bounds,
use the original P_i<=1 contract and retain full zeta. The two waves must
actually have the same Gaussian at the mapped point; a radial weight or
physical trace does not establish that condition. Different overlapping
phases require the original closed support theorem, not this one-label
calculation. No unconditional provider inequality follows from (3).

## 3. Exact changed radial source and moments, with means held fixed

Put D=Haar_aux(Delta W) for this paragraph, using normalized Haar measure.
The pointwise radial source changes by

    delta gr=-D1 Delta Wrr-Dz Delta Wrz+Delta Wtt/R,
    Haar(delta gr)=-partial_R Drr-Drr/R
                              -epsilon partial_Z Drz+Dtt/R.       (4)

The second identity uses the original fixed common torus on each full
radial fiber; Haar removes the graph derivative. It does NOT remove the
pointwise graph term in the first equation. The pointwise source has the
sufficient order C-kappa; the Haar source has order C. All coefficients
are complete covariances, with no transverse-only replacement.

Integrate over the full FIXED positive R hull [a,b] containing the actual
wave supports. Under the retained compact endpoint and mixed-jet support
contracts, integration by parts gives

    delta P=int (Dtt-Drr)/R dR
                           -epsilon partial_Z int Drz dR,
    delta Jtheta=int R^2 Dtheta,z dR,
    delta Jz=int R[Dzz-(Drr+Dtt)/2] dR
                       +(epsilon/2) partial_Z int R^2 Drz dR.     (5)

Before compact endpoint cancellation, the first line has the additional
-[Drr]_a^b, and the last has +[R^2 Drr]_a^b/2. If the hull or its endpoints
move in a different chart, retain the Leibniz boundary terms. Compactness
on individual mesh boxes or reserved bumps is insufficient; the full hull
and original fixed-q fiber are required. These are analytical identities,
not numerical quadrature observations.

For the ORIGINAL unit-mass fixed compensator rho and c_rho=(int R^2 rho)/2,
Haar(delta pm) satisfies partial_R Haar(delta pm)=Haar(delta gr)-rho delta P.
Its compact primitive obeys the global pressure-moment identity

    int R Haar(delta pm) dR
             =-1/2 int R^2 Haar(delta gr) dR+c_rho delta P.       (6)

It follows by integration by parts with the same compact pressure/support
conditions. Thus the wave contribution to the full axial scalar remains

    epsilon partial_Z(delta Jz+c_rho delta P).           (7)

The derivative hits c_rho as well. Its slow variation cannot be frozen.
P_changed=P_old+delta P, J_changed=J_old+delta J and pm_changed=pm_old+delta pm;
old mean/base pressure, Pi1/Pi2 and preceding moving moments are not reset.
The sourced wave pressure from Q96 is a nonzero harmonic and is not a
replacement for this compact angular-zero pressure memory.

## 4. Explicit moment coefficient bounds and their limits

Let Eedge(2q,width) bound zeta delta^-2q on the WHOLE hull, with the
unchanged width, and put L=(b-a)Eedge. For covariance raw bounds D_j,
(5) gives sufficient coefficients

    P_j=2L D_j/a+L D_(j+1),
    Jtheta_j=L b^2 D_j,
    Jz_j=2L b D_j+(L b^2/2)D_(j+1).                    (8)

They all have coarse order C. The first and last formulas drop only the
beneficial epsilon of an axial derivative, not a radial graph loss.
For the pressure moment (6), the core has coefficient
L b D_j+(L b^2/2)D_(j+1), then add the binomial product c_rho*P.
For (7), differentiate the whole Jz+c_rho P product: output through n
requests covariance through n+2 and c_rho through n+1. Its sufficient
order is C+1=5/2-kappa. Every c_rho,Z P term remains in this envelope.
The code derives the edge absorption from the new covariance pole; it
does not reuse an old pressure-domain bound by assertion.

These identities and moment estimates do NOT supply a full POINTWISE
pressure bound. The precise remaining map inequality is

    |D^I T0(delta gr-rho delta P)|
        <= epsilon^(C-kappa) zeta delta^-p_out K_I(S)    (9)

on the SAME whole shifted-torus/radial domain, with the current compensated
source at pole 2q and all operator derivative requests. T0 must keep its
original moving cutoff, shifted radial integrand and compensator. An old
operator proof whose pole domain does not contain this new request cannot
be reused. A small pressure MOMENT, or Haar graph cancellation in (4),
does not prove (9) or erase the full pointwise graph source. The fresh
fast/slow means and final radial compact source require (9) and their
own changed-state rows before they can be claimed.

## 5. Requests and frontier

SourcedCovariance accepts only Q96 SourcedResponse and constructs BOTH
positive response carriers, retaining their conjugates. Its covariance
bundle records the old state, common wave pole and new covariance pole.
To bound the whole axial scalar at output n, covariance order n+2 requires
Q96 source jets through n+3 and Q94 mean-factory reserve n+5, hence its
conservative current-target/profile reserves n+19,n+32 before full primitive
DAG expansion. These are finite analytical requests, not evaluated giant
physical constants or assertions uniform in order/stage.

Established for the current Q96 state: a typed composition of the frozen
complex covariance/global moment calculus with explicit conditional
coefficients, current packet costs and newly derived proof-domain poles. Next supply (9) with the new source pole and
whole pressure memory, then recompute changed scalars, incoming rows,
complete fast/axial fields and refreshed slow means. Subsequent signed
updates must use current target rows and the original root. First-real
compact tails cannot be transplanted. Actual native continuum contracts,
unconditional pulse bounds, complete residual improvement, infinite
realization and blow-up remain unproved. Exact finite tests and frozen
provenance do not constitute a continuum proof. Run `python test_audit.py -v`
and `python audit.py --verify-record evidence.json --output FRESH.json`.
