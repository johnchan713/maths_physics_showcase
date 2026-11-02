# Q98: current shifted pressure from the compact covariance flux

Parent: `e4609c5b345dfe6b940af4cdc7f7735e044361bd` (Q97).
Both parent workflows were still running without a reported conclusion
when inspected on 2 November 2026 (runs 41418521679 and 41418521802).
This checkpoint supplies a conditional sufficient coefficient for Q97's
pointwise pressure obligation on its NEW covariance pole. The flux split
gives the stronger increment order C=3/2-kappa; Q97's requested C-kappa
bound was sufficient, not false. The exact left-inverse/trace calculus is
inherited from core_compact_source_quotient and core_compact_trace_gate.
The contribution here is its current Q97 binding and finite coefficient
program, not discovery of those identities.

All mean velocities remain fixed at the incoming Q96/Q97 state. Nothing
here certifies actual native norms, whole-hull traces, fresh means, a new
complete residual rate, uniform physical pulse bounds or blow-up. Original
root, geometry, width C_setup^-20 and historical evidence are unchanged.

## 1. Preserve the local pressure memory before estimating

Let F=Delta Wrr be the FULL auxiliary-dependent changed covariance, and

    H=(Delta Wtt-Delta Wrr)/R-Dz Delta Wrz,
    delta gr=-Dr F+H,  delta P=int Haar(H) dR.             (1)

Here D1=Dr+1/R and Dz=epsilon partial_Z, with labels held fixed.
The -F/R connection is included in H. The last equality requires the
whole compact endpoint/mixed-trace contract; before cancellation retain
-[Haar F] at the exterior endpoints and any piecewise classical seam terms.
Haar cancellation does not establish a pointwise pressure estimate.

With U=R^d, V=partial_U+M vr.partial_y and f0=R^(1-d)/d,
one has f0 Dr F=V F. Integrate along the WHOLE shifted characteristic,
using the ORIGINAL moving cutoff chi and the same compact operators:

    I_minus(Dr F)=F, I_plus(Dr F)=-F,
    T0 Dr F=F, A0 Dr F=0.                               (2)

Thus the exact changed compact pressure and its radial residual are

    delta pm=-F+T0(H-rho delta P),
    Dr delta pm-delta gr=-rho delta P-A0(H-rho delta P).  (3)

The scalar rho is the same moving unit-mass compensator. It is neither
reselected nor frozen under slow derivatives. Equations (2)--(3) hold
pointwise in slow variables, so their mixed derivatives include cutoff
and rho motion. No commutation of Dz with T0 or A0 is used.
Pressure is not zero when the radial derivative is annihilated by A0:
its local contribution is exactly -F. Q97's old covariance, old compact
pressure, Pi1/Pi2, nonzero sourced harmonic pressure and previous moments
remain additive. Only the zero-angular pressure INCREMENT is bounded here;
the TOTAL pressure is not relabelled C.

For piecewise classical derivatives, on one output characteristic let
j_s=[F]_s include all internal and exterior jumps, J=sum j_s and
J_less=sum_(s<R) j_s. The actual replacement for the first line of (3) is

    delta pm=-F+(J_less-chi J)+T0(H-rho delta P).          (4)

Here delta P must also be recomputed: with exterior zero pieces included,
delta P=int Haar(H) dR+Haar(J). The Haar of a transported jump is the
same as its original fiber Haar. Also A0 Dr F=-chi_R J; hence the residual
acquires -chi_R J.
Two opposite jumps can have J=0 while J_less is nonzero. A zero Haar
moment or zero total jump therefore cannot justify deleting (4)'s defect.
For a C-order pressure bound without exact matching, the precise extra
requirement includes a SAME-weight through-n bound of order C for
J_less-chi J, the resulting T0(rho Haar(J)) pressure contribution, and
for the radial residual additionally chi_R J and rho Haar(J), with all
moving-seam mixed derivatives. Equivalently first enlarge P_n by a bound
for Haar(J), then propagate it through B_n and (6), and add the explicit
cumulative and projected jump bounds. A smooth global compact extension
makes these zero. This checkpoint does not certify that actual contract.
Discontinuous comparison fixtures are not admissible smooth velocities;
their distributional derivatives would also include interface measures.

## 2. Explicit conditional norm coefficient at the current pole

Use the five native coordinates (R,y1,y2,Z,T), factorial-through-n norm,
and Q97 raw mixed-jet coefficients D_j for EVERY covariance entry at the
common pole p=2q. Define

    F_n=sum_(j<=n) 5^j D_j/j!,
    r_n=sum_(j<=n) a^(-j-1),
    H_n=2 r_n F_n+(n+1) F_(n+1),
    P_n=(b-a) Eedge(p,width) H_n,
    B_n=H_n+rho_n P_n.                                  (5)

The 5^j/j! is the sum of 1/alpha! over length-j multiindices.
The normalized product norm is submultiplicative, and partial_Z costs
at most (n+1) times the through-(n+1) norm. Only its beneficial epsilon
was dropped. The reciprocal radius contributes r_n on the fixed positive
hull [a,b]. Eedge is the original full-hull weight absorption. rho_n is
the original weighted_rho(n) coefficient, including its weight floor.
Haar/integration gives P_n; it does not bound the shifted primitive.

Rebuild the SAME e=0 Q47/Q46/Q43 compact recipe at the EXACT new pole p,
with original exponent d, hull, shell, moving cutoff shape and width.
At output n its normalized constants KT_n(p), KA_n(p) give

    ||delta pm||_n <= F_n+KT_n(p) B_n,
    ||A0(delta gr-rho delta P)||_n <= KA_n(p) B_n,
    ||Dr delta pm-delta gr||_n <= rho_n P_n+KA_n(p) B_n.  (6)

All coefficients multiply epsilon^C zeta delta^-p under the retained
native weighted derivative/support and compact operator contracts. The
two wave Gaussians have already been treated with Q97's P_i<=1 contract;
no derivative of a weight inequality is inferred. Coefficients are exact
positive exponential-polynomials in the original S program, not evaluated
numerical enclosures. The operator chart conversions preserve the pole
and introduce no epsilon loss. This improves the sufficient pressure
increment bound over directly estimating Dr F at order C-kappa. It does
not improve that pointwise SOURCE estimate or assert the compact residual
is flat. Retain the A0 and scalar terms in (6).

ShiftedPressure accepts a Q97 SourcedCovariance factory, uses its SAME
profile and original d, and rejects a different source state/stage,
updated means, insufficient jets or an out-of-domain pole. The supplied
c0,c1,shape and inverse recipes must be those of the original construction;
their analytical validity is still a lower contract. An older compact
archive cannot be enlarged by presumed pole monotonicity. The exact new
recipe is requested explicitly. No Haar-only Q38 map is used for T0.

## 3. Requests, evidence and next obligation

Output n requests covariance through n+1, hence Q96 source through n+2,
Q94 mean-factory reserve n+4, and its conservative target/profile reserves
n+18,n+31 before primitive DAG expansion. The combined KT/KA interface
requests cutoff shape through n+1. Q97's WHOLE axial moment derivative
still needs covariance n+2; this pressure bound does not lower that separate
request. No huge physical request or actual norm has been numerically
enclosed.

Exact tests integrate a nonzero shifted Fourier flux, retain a moving chi,
verify (3) and its slow derivative and keep a nonzero rho P term. Controls
detect deleting the local pressure, the 1/R connection, or an internal
cumulative jump. Manufactured norm leaves check binding and rejection;
they are not canonical continuum data. Frozen provenance is replayed
without editing historical files. Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json`.

Next recompute changed scalars, incoming rows and complete fresh fast/axial
and slow mean updates from this changed covariance AND retained pressure,
at their new poles. Recheck global moments/reserved supports after those
updates before using a subsequent signed target. Actual compact trace,
native input and all-scale pulse contracts, complete viscosity/cutoff/
interaction residual bounds and infinite realization remain open. Passing
finite checks is not a continuum proof.
