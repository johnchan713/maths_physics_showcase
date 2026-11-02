# Q105: current projected final rows and mean-wave feedback

Parent: `f8c139d56ef53c64e1f04851cbb9718ee9b35afa` (Q104).
Its research and field-comparison CI runs 41489933754 and 41489933720
were still running when inspected on 2 November 2026. Q102's
field-comparison run 41465339945 had passed. This checkpoint binds the
frozen compact-source quotient and complete mean-difference identities
to the CURRENT Q99--Q104 state. The identities themselves are inherited;
the source-stage binding, coefficient composition, and remaining gates
are the new result. All estimates below are conditional on the inherited
native norm, whole support, mixed trace, and physical-cover contracts.

## 1. The axial remainder belongs to the state BEFORE the new fast update

Let O denote the first four mean contributions in Q104: initial fast/slow
and first-real fast/slow. Let F and S denote the new Q101 fast and Q102
slow increments. Let W denote the TOTAL old-plus-sourced complete-wave
covariance from Q97, with both cross orientations and the sourced square.
The wave remains fixed while O becomes O+F+S. Write Epre for the full rows
on (O,W), and Qcirc=I-Haar for the auxiliary-torus projection. The inherited
exact fast cancellation gives

    Efinal_theta = Haar Epre_theta + RF_theta + RS_theta,
    Efinal_z = Haar Epre_z + A1 Qcirc Epre_z + RF_z + RS_z, (1)

where RF and RS are the full mean differences EXCEPT c Nt of their
tangential increments. RF uses O,F; RS uses O+F,S. Both contain the actual
axial fast velocity gd-A1 gd, all nonlinear mean products, time, vector
viscosity, connections, and the full shifted pressure difference. The
slow increment is auxiliary-independent, so its c Nt term vanishes.

The A1 source in (1) is CURRENT Epre, after the sourced wave but BEFORE
F and S. It is neither the original primary prefast tail nor the final
postslow radial source. `pre_fields` assembles exactly O without calling
the new mean factories. The pressure compensator is Q99.ChangedScalars
(old completed postslow scalars PLUS sourced-wave changes), not Q103's
post-new-slow P. This keeps complete pressure history at the right stage.

## 2. Compact potential identity, with both viscous terms

Use the same Dr, D1=Dr+1/R, Dz=epsilon partial_Z and tstar as Q104.
Label parameters are fixed, so [tstar,D1]=[Dz,D1]=0, as required in the
frozen quotient. A variable radial c would add -Dr(c) Nt Psi to (3)
and would require a separate remainder; it is not covered by this identity.
The old compact mean potentials sum to Psi, with beta=-Dz Psi and
gamma=D1 Psi. On the ENTIRE original potential hull put

    S0=2b beta+beta^2+Wrr,
    H=b gamma+G beta+beta gamma+Wrz,
    r0=-tstar beta-S0/R-Dz H+(2Vv+v^2+Wtt)/R
                                                +epsilon Dz^2 beta,
    pm=-S0+epsilon D1 beta+T0(r0-rho Ppre),
    H1=2G gamma+gamma^2+Wzz-S0+T0(r0-rho Ppre).          (2)

The full axial row, including the prescribed axial stress Sigma_z, obeys

    Epre_z-Dz H1 = D1 Omega,
    Omega=tstar Psi+H-epsilon Dr gamma
                              -2epsilon Dz^2 Psi-Sigma_z. (3)

Indeed Dz(epsilon D1 beta)=-epsilon D1 Dz^2 Psi, while
-epsilon Dz^2 gamma=-epsilon D1 Dz^2 Psi is a SECOND equal term.
The coefficient two in (3) is essential. Local viscous pressure has been
retained before factoring; -S0 and the entire shifted primitive in H1
also remain. No complete-wave pressure or axial term is deleted.

Under the retained whole-hull compact trace/extension contracts,
A1 D1 Omega=0, hence

    A1 Qcirc Epre_z = A1 Dz Qcirc H1.                    (4)

This uses the same original moving operator, evaluated at its EXACT pole.
It does not commute Dz through A1. For a moving chi,
[partial_Z,A1]f=R^-1 chi_RZ J1 f generally survives. Separate compact
bump supports do not prove the whole Omega endpoint identity; a seam
defect must be kept if that contract fails, as specified in Q98.

The implementation bounds H1 through n+1 from O and total W through
n+3. It uses the original T0(n+1,p), rho through n+3, Q99 Ppre(n+1),
and the Q101 archive A1(n,p), with upward coefficient-pole embeddings
only. If L=pre.reduced_source and Z bounds gamma, its positive coefficient
is

    H1_bound = 2 G Z+Z^2+W+S0_bound+KT0(L+rho Ppre),
    tail_bound = 2(n+1) KA1 H1_bound.                    (5)

The factor two covers Qcirc. The explicit epsilon in Dz is retained in
the analytical grade; positive coefficient arithmetic may drop it since
0<epsilon<=1. Q104's coarse grade B=1-kappa for (2) implies the improved
conditional tail grade B+1=2-kappa. This is the already established
compact-source quotient gain, now bound to the current pre-F state.

## 3. Projected final row coefficients and exact limits of the join

For the same common pole p, multiply weighted mean increments by
1+Eedge(p,width) before products. Embed each at its own pole; unweight
F at its OWN returned pole before forming the old full O+F coefficient
for RS. Q86 supplies complete unweighted O. Use the frozen
core_fast_row_gain_join.difference at through-(n+3) inputs, including
the same KT0(n+1,p), rho and whole-hull moment constant. Then

    ||Qcirc Efinal_theta|| <= 2(RF_theta+RS_theta),
    ||Qcirc Efinal_z|| <= 2(tail_bound+RF_z+RS_z).         (6)

The grades follow from D=3/2-2kappa for F, C=3/2-kappa for S,
their radial grades D+1,C+1, and the inherited mean-difference rule
increment_grade+1-2kappa:

    RF: 5/2-4kappa; RS: 5/2-3kappa;
    Qtheta: 5/2-4kappa;
    Qz: min(2-kappa,5/2-4kappa).                         (7)

The min matters for the full retained range 0<kappa<1/4: 2-kappa
dominates only for kappa<=1/6. At the original kappa=1/100000 the axial
grade is 2-kappa. No unconditional estimate is inferred from (5)--(7).

Equation (1) still contains the POINTWISE radial functions Haar Epre.
The small global integrals P,Jtheta,Jz in Q103 do not bound these
functions. For example, on [1,2], let phi=(R-1)^3(2-R)^3 and take
e_theta=R^-2 phi', e_z=R^-1 phi', auxiliary-independent. Both relevant
weighted integrals vanish and Qcirc e=0, while both functions are nonzero.
This is a logical comparison control, not a canonical-source counterexample.
Thus (6) supplies projected rows only; it must not be advertised as a fine
bound for the COMPLETE pointwise rows. Q104's surviving radial A0 and
scalar-compensator gates are also unchanged.

## 4. Mean-update feedback against the COMPLETE current wave

Unweight the new fast F at its own pole and add unweighted S, with all
three components. Let M bound that sum and U bound the complete old-plus-
sourced wave through n+3 in the ORIGINAL half-edge weight with each label's
Gaussian. Reuse the SAME phase normal and graph from Q93's current wave
object; no replacement normal or outgoing residual callback is introduced.
Q91's mean-wave polynomial, with harmonic=2 for both +/-1 and +/-2,
gives the positive coefficient

    3[2(n+1)(g+1)+4N+2r] M U.                           (8)

Both orientations (F+S).grad(w) and w.grad(F+S), axial velocity,
connections and all cut-curl wave components are covered. Unweighted
means times the half-weighted wave retain that half weight and Gaussian;
no same-pole weighted mean squaring is silently used. The wave's epsilon
order 1/2 offsets its worst fast phase derivative epsilon^-1/2. Thus this
mean-update feedback has grade min(D,C)=D=3/2-2kappa.

This is ONLY the new mean-wave feedback group. It does not supply a bound
for the entire regular residual after the sourced solve: wave-wave
interactions (including harmonics +/-3 and +/-4), linear and viscous
response, cutoffs, new pressure errors, and their cancellations remain
separate obligations. In particular (8) cannot upgrade the whole regular
source merely by changing its stage label.

For the illustrative target t=3-6kappa, the available projected angular
gap is t-(5/2-4kappa)=1/2-2kappa, the axial gap at the original kappa is
t-(2-kappa)=1-5kappa, and the feedback gap is t-D=3/2-4kappa. These are
gaps in AVAILABLE UPPER BOUNDS, not lower bounds proving target failure.

## 5. Reserve, checks and next obligation

Output n needs O,F,S and total wave through N=n+3, final radial pressure
archive through n+1, current Q101 row reserve n+10, covariance n+12,
Q94 mean-factory n+15, conservative target/profile n+29,n+42, and the new
cutoff shape through n+12. Older operator shifts, primary tails, whole
physical covers, all high-jet request DAGs and exact-pole contracts remain;
this checkpoint does not evaluate enormous actual jets or certify their
native provider inequalities.

The tests independently differentiate (3), detect a missing factor two or
deleted local viscous pressure, check the moving-A1 commutator and whole
characteristic boundary cancellation, distinguish projected and global
moment controls from full Haar rows, and exercise current source, order,
pole and phase bindings with manufactured norm leaves. The feedback
coefficient is checked against direct exact rational arithmetic. Replay
recursively verifies the untouched frozen evidence and source hashes.
None of these checks constitutes a continuum proof. Run
`python test_audit.py -v` and `python audit.py --verify-record evidence.json
--output FRESH.json` in this directory.

Next derive the COMPLETE pointwise Haar rows in (1) from the current full
covariance, pressure and all mean histories, or exhibit the precise
inequality preventing the required bound. Join them with (6), the inner/
outer estimates and Q104's radial gates before claiming full cone closure.
The global pressure/moment identities, coefficient and reserved-mean
supports, uniform actual-pulse estimates on every physical scale, and
complete viscous/cutoff/interaction residual still require their full
actual contracts. There is no unconditional pulse bound, full residual
improvement, infinite realization or blow-up theorem. Original root,
accepted domain, flat width C_setup^-20 and historical evidence are intact.
