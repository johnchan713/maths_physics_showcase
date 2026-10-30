# Q27: close the four complete-wave inputs of the first signed source

Parent: `1fcc2669852829f108e411c5648017a84a6c9f4e` (Q26).
The first signed-prefix source recipe takes four independent wave norm
callbacks T, R, Pi and U. This checkpoint replaces those callbacks by one
explicit construction from the SAME original canonical root, its actual
Haar signed target, homogeneous pulse jets, full curl operator and phase.
It does not certify the remaining lower premises or the complete source.
The analytical result is a conditional estimate; tests use manufactured
jets. Historical evidence, the field, selector, exact h and width stay fixed.

## 1. Domain, norms and two cutoff norm slots

Work on Q26's local domain Dv and every required full radial/auxiliary fiber.
Use Q25's root floors and Q23's upper propagator only with their stated
same-profile hypotheses. The canonical entrance and Haar-column jet
providers supply the original root budget's H and wave slots. The signed
input is the current Haar target Xi=epsilon^alpha zeta xi,
alpha=2-3kappa, 0<kappa<1/4. The target-edge and initial-Haar-row providers
supply its conditional jet recipe. Choose one finite d covering its
J_j(xi)<=X_j delta^(-3j-d) at ALL requests in this calculation.

J_j is the sum of factorial-normalized jets of exact total order j; the
through-order norm is their sum through n. Retain component maximum and
matrix induced infinity norms in the root recipe. Full phase component/
entry sums dominate the contractions below. All labels are frozen.

There are two wave norm slots from the SAME homogeneous vector t_sigma:
Wreg_j bounds J_j(chi*t_sigma)/P in the retained regular linear grammar;
Wcut_j bounds J_j(chi*psi*t_sigma)/P for the actual temporal cutoff.
Here chi includes the other original cutoffs and is constant along the
auxiliary pulse flow, as required for principal homogeneous cancellation.
Obtain both using the
canonical entrance provider's cutoff convolution. In particular a root
budget already built with chi*psi is not automatically the regular budget.
Rebind only the wave NORM slot to form these two templates. H, original
leading target, inverse/root gaps, mass, loghalf and complete curl operator
are identical. This changes no physical root or entrance.

The temporal cutoff is applied inside the potential before taking its full
curl. Thus the actual cut amplitude is not psi times the complete regular
velocity. This distinction is enforced below by differentiating Wcut through
one extra order. It preserves the curl-of-cutoff velocity in every nonlinear
source, while the regular linear commutator remains an additive tail.

## 2. Principal and curl polynomials, including the signed addition

For either template let A_j be `RootBudget.amplitude(j)` and Z_j be its
`RootBudget.signed(j,xi)` from the original-root formula

    delta a=H^-1(Xi/epsilon)/(2a), a=sqrt(zeta/s0)*sqrt(b).

Set B_j=A_j+Z_j and, with the original full curl operator polynomials D_i,

    C_j=sum_(i=0)^j D_i [B_(j-i)+(j-i+1)B_(j-i+1)].    (1)

This is the exact linearity of the potential and curl at the norm level;
it includes derivatives of the signed coefficient, reciprocal original
root, cutoffs and D_i. It never replaces a by sqrt(a^2+delta y).
No positivity of the updated signed amplitude is required. Xi is computed
from the initial pre-signed mean state and then held fixed during this
addition. Feeding a post-addition target back into (1) would be a different,
circular prescription; the subsequent fast/slow means are recomputed only
after the wave has been formed.

The original primary has principal epsilon order 1/2 and complete curl
remainder order mu=1-kappa. The signed addition has principal order
alpha-1/2 and remainder order alpha-kappa. Since alpha>=1 and
alpha-kappa>=mu, (1) bounds the combined principal and remainder at the
coarser orders 1/2 and mu for epsilon<=1. Their common radial/Gaussian
weight is sqrt(zeta)*P. Through order n a sufficient pole is

    p_n=max(6n+9,6n+d+6).                              (2)

Indeed the primary principal/curl poles are 6j+3 and 6j+9; the signed
principal/curl poles are 6j+d and 6j+d+6. Since delta<=1, (2) pads all
j<=n. Cutoff jet constants are included in the wave slots, not assigned
an additional epsilon gain.

For two canonical signs, each real cosine/sine Fourier pair has total
coefficient mass one. Sum both signs, keeping any coincidences by triangle
inequality. The four source callbacks will use

    T_n=2 sum_(j<=n) Breg_j,
    R_n=2 sum_(j<=n) Creg_j,
    U_n=2 sum_(j<=n) (Bcut_j+Ccut_j).                   (3)

U bounds the ACTUAL cut complete wave at epsilon order 1/2, using
mu>=1/2. No conjugate, signed square, old/new cross orientation or cutoff
curl is deleted. Further label/common-atlas conversion factors, when
needed, must be applied to these norms with the existing separated-support
atlas provider; the factor two in (3) counts signs, not arbitrary labels.

## 3. Derive the pressure slot from the full phase

Let N_j, K_j, Nprime_j be the raw component/entry-sum phase bounds, and
q0 a certified upper bound for 1/|n|^2. For raw convolution
(f*g)_j=sum_i binom(j,i) f_i g_(j-i), set

    Q_0=q0,
    Q_j=q0 sum_(i=1)^j binom(j,i)(N*N)_i Q_(j-i),
    L_j=(Q*(N*K+Nprime))_j,
    Lthrough_n=sum_(j<=n) (5^j/j!) L_j.                (4)

Differentiating |n|^-2 |n|^2=1 proves the reciprocal recurrence.
The exact homogeneous pressure is

    k*pi=i (n^T K-nprime^T)t/|n|^2.                   (5)

Multiplying by the original root or its signed coefficient preserves (5),
because both are independent of auxiliary pulse time in this FIRST signed
step. Thus a sufficient source pressure polynomial is

    Pi_n=Lthrough_n*T_n.                              (6)

The factor 5^j/j! converts raw mixed jets to normalized total jets in the
five-variable chart. Through-order product submultiplicativity proves (6),
including every derivative of the pressure row and principal amplitude.
Since k^-1<=sqrt(epsilon), Pi has epsilon order 1 and the same half edge
weight and Gaussian. It has the pole (2), assuming the phase/floor jets
are unweighted polynomial bounds as in their provider.

The moving-normal term in (5) is essential even with K=0. A purely
homogeneous pulse has no source-normal term because its source is zero;
this fact cannot be used for the next forced pulse, whose pressure retains
+n.f. Nor does (6) replace the compact mean pressure, Pi1/Pi2 background
memory, axial A1 velocity, A0 additive pressure error or moving c_rho P
moment. Those remain in the unchanged complete mean/source recipe.

## 4. Substitute into the first source without a new wave hypothesis

`SignedWaveInputs.source_wave` supplies exactly T,R,Pi,U to
`core_signed_prefix_source_budget.PrefixBudget`, preserving exact positive
polynomial coefficients across historical helper modules. Equations
(1)--(6) remove those four INDEPENDENT callbacks, relative to the listed
root/target/phase/ODE/cutoff/curl premises. The unmodified recipe recomputes
fresh fast means, then full pressure and all scalars, then the slow response;
its complete source polynomial is still Flin+Fww+Fmw.

For source order j the conservative fresh-mean route can request wave
inputs through j+15. Equation (1) consequently requests root and xi jets
through j+16. This is BEFORE expanding the target's own initial-Haar-row,
scalar, root and background dependencies; it is not the total primitive
ceiling. Use Q14's fully expanded requests and all distinct pole constants.
Do not silently reuse a target pole d chosen only for order j.

At a finite maximum wave request N, choose p>=p_N also covering the
initial mean input poles. The unchanged source recipe then uses its
2p/4p operator archive and gives source pole 3p. Subject to the retained
initial mean and fixed operator estimates it yields

    |D^I f_m| <= epsilon^mu sqrt(zeta) delta^(-3p) P
                       * j! F_j(S), |I|=j.            (7)

Q20--Q21 can convert (7) to whole-path f/P jets. The result here supplies
the wave part of that implication; it is not a claim that its remaining
fixed-map, full-background, support or actual target premises are proved.
In particular original selector input validity is not inferred from Q25
and Q26's conditional implications or from this executable composition.

Next close the remaining first-source initial-mean and fixed-map input
composition on Dv, with the actual weighted-stress and finalized-background
schedule and all domain hypotheses. If a premise cannot be supplied, state
its specific inequality instead of treating finite symbolic polynomials as
actual norm certificates. Uniform physical-scale representative bounds,
historical localization, full residual improvement and blow-up remain open.

## 5. Verification

Tests exercise the unmodified source factory with the new wave adapter,
trace its j+15/j+16 requests, preserve original root/value objects, and
check independent exact pressure, signed-square and cutoff-curl controls.
They are finite manufactured controls, not selected-field evaluations.
The provenance audit replays Q26 and its frozen dependency graph. Run
`python test_audit.py -v` and `python audit.py --verify-record evidence.json
--output FRESH.json` here. Passing tests is not a continuum proof.
