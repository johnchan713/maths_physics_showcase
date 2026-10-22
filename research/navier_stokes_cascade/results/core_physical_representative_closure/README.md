# Original physical representatives close two extra floor requirements

Parent: `134d9a06a9c2e96ecb1b4debf7c1e7e8888bb1cc`.
This checkpoint verifies where the ORIGINAL representatives are frozen.
It removes an outside-physical-set representative concern and a separate
native value-neighborhood requirement from the tangential NORMAL FLOOR.
No representative, band, support, carrier, tilt or width C^-20 is changed.
The actual unchanged-band floor and complete residual remain unproved.

## 1. The source chooses an active point, not an arbitrary box point

The primary manuscript's section 6.2, equation (6.8), printed page 64,
defines an active chart set A_ell by 0<q<qbig, 1/2<=q/Q<=2, tau>=0,
and Xa<=r^2/(2q)<=Xb. Its representative satisfies

    x0_(ell,a) in B_(ell,a) intersect A_ell.                    (1)

Thus every representative has T0>=0, s0=q/Q in [1/2,2], |eta0|<=1,
and X0 in [Xa,Xb]. This includes boundary representatives at T0=0.
The choice is made before derivatives and used for both signs. A mesh
center or an arbitrary original-box point need not satisfy (1).

Section 7.1, printed page 74, uses physical relative neighborhoods and
one-sided boundary derivatives. This does not provide quantitative
negative-time extension bounds. Sections 6.2 and 5 retain a separate
smooth-extension obligation for assembled supports/fields.

See [the primary source](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
`source_contract.json` records these definitions, locations and the PDF
hash in historical metadata. The PDF was read directly for this checkpoint;
its binary hash was not independently recomputed. No theorem claimed by
that manuscript is used as a premise; only its defining point-selection
and physical-domain convention are identified here.

The older audit's statement that a representative is “in its original
box” is true but weaker. The source's intersection clause resolves the
conditional concern in `core_native_chart_cover`: an ORIGINAL representative
cannot lie at T<0 or outside the native active annulus. Historical records
remain frozen; arbitrary representatives from a different construction
would still require their own coverage proof.

## 2. The existing physical doubled-box norm hull applies from that point

Set delta=ell^-6. By (1), every y in the same doubled box has displacement
at most 3delta per slow coordinate from its actual frozen representative.
For this floor use only the physical part T>=0. The parent proves on the
same physical branch s-Z^2 s^(2h)=T that

    |partial_Z log s|<5, |partial_T log s|<5
       while 1/4<=s<=4 and |eta|<=1.                         (2)

The straight segment from x0 to any such y stays physical. A first-exit
argument gives |log(s/s0)|<=30delta<log 2 when 100delta<=1;
therefore 1/4<=s<=4 on the whole segment. This proves a whole-box bound,
not a test on finitely many chart points.

With fixed outward enclosures rA<=sqrt(Xa), rB>=sqrt(Xb), if

    6delta<=rA, 100delta<=1,                                 (3)

the same original norm hull is

    R in [rA/2,5rB/2], Rnative in [rA/4,5rB],
    Xnative <=25rB^2/2.                                     (4)

The barrier gate in (3) is automatic for ell>=16. The original geometry's
delta<=sqrt(Xa)/6 implies the clearance gate with rA=sqrt(Xa) as an
EXACT REAL enclosure. A stricter rational outward enclosure must still
pass its supplied gate; no band is increased here to repair failure.

All reference points now lie on the closed active annulus. The existing
leading/frame floor theorem therefore supplies their reference denominators,
relative to its retained core/regional proofs, with the old compact s
scaling factors. The much larger hull (4) needs UPPER VELOCITY NORMS only.
No growth, stress, orientation or cone floor is asserted on that hull.
All norm and gap bounds remain relative to their actual profile premises.

## 3. A native value-neighborhood gate is unnecessary for this floor

The complete tangential normal is

    n_tan=(p/R,pz-epsilon*v*H_Z), H=pF+pzG.                    (5)

Its reference is Bs Kunit at the ORIGINAL representative. The preceding
floor proof uses positive reference growth/shear bounds ONLY THERE.
Elsewhere in the physical doubled box it uses the radius bounds (4),
the carrier rounding cost, |R-R0|<=3delta, and a global upper bound for
the FULL H_Z. It does not use a pointwise growth/shear/cone lower bound
at an evaluation point outside the active annulus.

By the retained Borel prescription and its existing coefficient-regularity
premises, the complete normalized correction tuple has first jets <=q^h.
With the SAME leading first-tuple norm on (4), the previous construction
therefore gives Jfull1=(CF1+CG1)(Lead1+1) and

    |n_tan-Bs Kunit| <= ct/S+ca/k+cb/S^3+cz*S*epsilon.          (6)

The four coefficients are unchanged from `core_tangential_normal_closure`;
cz bounds the complete axial transport. No local U=0 pressure reset,
cutoff-curl omission or leading-only substitution occurs.

Consequently this floor needs NO additional prescribed log-native width
mu, and NO check delta(12/rA+30)<=mu. That check still matters for OTHER
estimates that actually require value/growth/cone bounds on a larger
native neighborhood. Likewise no negative-T representative or extension
is needed to apply (6) on T>=0. This does not discharge the global
smooth-extension, full cone, covariance or residual construction.

The actual remaining floor inequalities on the SAME band are precisely

    S=ell_existing^2, S^2(epsilon+epsilon^2+k^-1)<=1,
    4 Ctan^2<=beta*S^2, beta=lambda_min/[4(1+U^2)^2],           (7)

along with (3), the original q-domain and justified input bounds.
`floor_join` accepts only the old band, derives its norm radii, and reports
the supplied normal-margin inequality. It never chooses a value width or
changes an old scale. Dbox>=6 is required for the retained doubled-box
Euclidean displacement bound. Failure is a failure of this sufficient
certificate, not a fluid counterexample or proof of nonexistence.

Removed requirements FOR THIS FLOOR: separate physical/native membership
of the original representatives, a native value-neighborhood width, and
negative-T representative coverage. Added analytical assumptions: none.
Retained: the original representative rule, actual leading/coefficient
regularity and reference gaps, full norm hull estimates, original slow/
normal/clearance/q-domain gates. No new representative is substituted.

Next: substitute the actual leading first-tuple and reference gap bounds
into (7) at the existing band, or document its exact remaining margin.
Do not continue requiring a separate mu or negative-time representative
certificate for this floor. Complete cone, sourced covariance, pressure/
axial mean and physical residual estimates remain separate obligations.

## 4. Reproduction and exact scope

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output fresh.json`. Exact controls include both physical
poles, native shell edges, a selected boundary box with a negative-time
center that is NOT its representative, malformed supplied points, and
the same normal-floor recipe with an arbitrarily narrow unrelated mu.
The rational chart check uses s=t^(2*denominator(h)) to avoid numerical
root errors; it checks supplied points and never chooses actual ones.
Finite fixtures are not continuum or actual-profile proofs.
