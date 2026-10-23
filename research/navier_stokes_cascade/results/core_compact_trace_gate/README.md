# Compact trace defects and the physical A1 support gate

Parent: `c9269c3412d4284b4c0b72ad5dd6663427e285ad`.
Its axial-potential/pressure-complement CI step passed when inspected on
23 October 2026; the full workflow was still running without a failed step.

The parent gives the conditional SAME-POSTWAVE order-two A1 error from
`Q Phi=-2G_R Q Psi+Q(gamma^2+Wzz)+Q pm`. Its next obligation includes
the actual whole compact support/extension contracts. Here we quantify
what those contracts must supply if a source is represented in several
closed pieces. We derive the exact seam defects in compact annihilation,
pressure reconstruction, and repeated integration by parts, including
moving seams. This replaces the vague support gate by explicit trace
identities or a blocking inequality. It does NOT establish them for the
canonical fields or find a defect in their retained smooth-extension premise.

The analytical conclusions below apply to classical derivatives on pieces.
A genuinely smooth global compact field has zero jumps and satisfies the
parent's identities unchanged. Discontinuous comparisons below are not
admissible continuum velocities; their distributional derivatives would
contain additional interface measures. Finite-regularity compact fixtures
are not the actual C-infinity sources. No field, pressure, root, width
`C_setup^-20`, physical scale, accepted band or historical evidence changes.

## 1. Exact primitive defects on one physical characteristic

Fix the original label M,c,epsilon, positive radial exponent d and
`U=R^d`. Use the fixed characteristic coordinate `v=y-M U vr`.
For a piecewise smooth potential h define

    H_e(r,v,Z,T)=r^e h(r,v+M r^d vr,Z,T).

Then, exactly, `r^e D_e h=partial_r H_e` on each piece. Let the
potential vanish at the two exterior endpoints and let
`j_s=[H_e]_s=H_e(s+)-H_e(s-)` at each internal seam. All occurrences
of H and j in the next formulas use the SAME output characteristic v.
Classical integration gives

    int_all r^e D_e h dr = -sum_s j_s,
    A_e D_e h = -R^-e chi_R sum_s j_s,                         (1)
    T_e D_e h = h-R^-e [sum_(s<R) j_s-chi sum_all j_s].         (2)

The SAME original chi appears in both terms of (2). The sums include
every seam inside the WHOLE hull, even in a gap between generating bumps.
Nonzero exterior values add their endpoint terms; equivalently include
the exterior zero pieces as seams. Multiplication by r^e and characteristic
transport are essential: raw chart values at different labels are not
the quantities j_s. No M,c or epsilon derivative was taken.

Thus zero total jump may suffice for (1) but not for (2): a cumulative
internal jump remains before a later compensating jump. Pressure cannot
be reconstructed from only a global zero moment and local derivative
structure. Exact matching at each seam restores both identities. For
e=1 this is also the actual `T1 gamma=Psi` identification; it is stronger
than only checking divergence and the two Haar velocity moments.

For the parent's radial decomposition `gr=r1+D_r K1`, the complement of
the complete reconstructed pressure is

    Q pm = Q K1+T0 Q r1
           -[sum_(s<R) [Q K1]_s-chi sum_all [Q K1]_s],        (3)

in the e=0 characteristic coordinates, if one uses only piecewise
classical derivatives. Here the current scalar compensator drops ONLY
because `Q(rho P)=0`. The full pressure still contains the same moving
rho P, local K0/K1, variance, six covariance entries and all earlier
pressure memory; the base Pi1/Pi2 has not been set to zero. A jump with
nonzero auxiliary Fourier mode has zero Haar and can evade both Haar
velocity moment checks while changing (3). Such a discontinuous example
does not satisfy the smooth-velocity premises of the construction.

Similarly the parent's complete axial normal form is

    Ez=D1 Omega+Dz Phi,
    Omega=Lambda+2G Psi,
    Lambda=tstar Psi+Hrz-epsilon D_r gamma
                              -epsilon Dz^2 Psi-Sigma_z.

Consequently the projected axial identity has the additional term

    A1 Q Ez=A1 Dz Q Phi -R^-1 chi_R sum_s [r Q Omega]_s.       (4)

It vanishes under the parent's actual whole-hull trace contract. Q removes
the auxiliary-independent assigned stress; no assertion of its separate
compactness is required. Both axial viscous contributions, including the
one in full pm, remain as in the parent. No moving inverse is commuted
through Dz. Equations (1)--(4) are a trace audit of the conditional identity,
not an asserted correction to the actual smooth field or its frozen record.

## 2. Higher traces in the nonstationary compact tail

For nonzero auxiliary Fourier mode k put `alpha=vr.k`,
`omega=2pi M alpha`. The characteristic moment has the radial integral

    I_k=sum_j int_(a_j)^(b_j) F_j(U) exp(i omega U) dU,
    F(U)=R^(e-d+1) f_k(R)/d.                              (5)

The original output factor is `R^-e chi_R exp(-i omega U)`;
it must not be differentiated as a fixed constant in higher output jets.
Fixed multipliers belong INSIDE F before taking traces, with their
original profile bounds. For the new A1 route `f=Dz Q Phi`, so F includes
the epsilon and the full partial_Z flux, not just Phi.

At every seam, including the exterior zero pieces, define
`Delta_r(s)=F_left^(r)(s)-F_right^(r)(s)`. Repeated integration by parts
on the individual pieces, with these signs, gives

    I_k = sum_s exp(i omega s)
                  sum_(r=0)^(p-1) (-1)^r Delta_r(s)/(i omega)^(r+1)
          +(-1)^p (i omega)^(-p)
                  sum_j int F_j^(p)(U) exp(i omega U) dU.    (6)

Artificial seams cancel when their transported weighted jets match; no
flatness of each separate piece at an internal seam is required. Checking
only that the source is zero outside a hull does not prove this matching.
For a given finite p, zero traces through p-2 suffice for an M^-p bound
at output zero, PROVIDED the remaining order-(p-1) traces and the bulk
p-jet have their required norm bounds. That last trace can be bounded
at the bulk order, so demanding zero traces through p-1 is unnecessary
for this single estimate. For p=1 no lower trace is forced to vanish.
The parent's smooth global contract of course implies all are zero.

In the original fixed-mode scaling, write a for the epsilon order of f,
and suppose a surviving Delta_r has additional order b_r. The two original
flooring bounds for M give, up to the retained fixed S-polynomials,

    bulk: epsilon^(a+p kappa),
    r-th trace: epsilon^(a+b_r+(r+1)kappa).                 (7)

A sufficient boundary contribution condition for target t is

    a+b_r+(r+1)kappa > t                                  (8)

for each trace term, or an EXACT cancellation identity for the combined
boundary sum. Retain the original bounded-overlap count, or bound the SUM
of all trace coefficients by a band-uniform S-polynomial; individual bounds
with an uncontrolled number of seams would not suffice. If any exponent is smaller, the separate positive norm
argument cannot certify that target. This alone is NOT a lower bound on
the sum: distinct seam phases may cancel. Equality is left open to extra
coefficient information; a strict margin is the existing sufficient rule
that absorbs arbitrary fixed S-polynomials. Summing modes and taking
higher output derivatives still requires the ORIGINAL Diophantine,
cutoff, raw-profile, weight and operator programs. A bound on one mode or
on finitely many seam samples is not their substitute.

An exact comparison shows that lost traces can actually matter. Convolve
n unit boxes and shift to a positive radial interval. The resulting compact
B-spline is C^(n-2), with its first derivative jump at n-1. Its Fourier
integral is

    exp(i omega) [(exp(i omega)-1)/(i omega)]^n.

At the fixed torus mode alpha=1 and half-integer M, its magnitude is
`2^n/|omega|^n`, nonzero. Multiplying by the reciprocal positive kernel
in (5) gives a compact source with the same finite regularity, and a
nonzero auxiliary cosine gives zero Haar. Hence these piecewise inputs
do not admit arbitrary flatness from support and lower jets alone.
This is a comparison of hypotheses: it is not smooth canonical data,
and its exterior edges do not satisfy the original flat-weight class.
It does not impose the full current scalar/mean/covariance history or
show failure of any actual source bound. The half-integer M sequence is
a diagnostic kernel comparison, not the original rounded-carrier sequence.
A fixed smooth
mollification restores smooth extension; its higher constants must then
be bounded, rather than inferred from the finite-regularity comparison.

## 3. Moving seams cost a further radial trace for Dz Phi

At `U=s(Z,T)`, differentiation of a trace identity gives

    partial_Z Delta_r(F)
       =Delta_r(partial_Z F)+s_Z Delta_(r+1)(F).            (9)

For a moving integral, the classical piecewise derivative likewise has
the extra contribution `s_Z Delta_0(F) exp(i omega s)`.
A continuous moving hinge has Delta_0(F)=0 but nonzero Delta_1(F);
then its partial_Z derivative has a nonzero value jump. Thus continuity
of Phi does not by itself supply value matching for Dz Phi.

Sufficiently, require the weighted Phi traces to match IDENTICALLY in
slow variables through order p-1. Equation (9) then makes the weighted
partial_Z Phi traces through p-2 vanish. This assumes the needed one-sided
mixed traces exist and their original norm bounds hold; it is not inferred
by differentiating a finite list of sampled zero values. An equivalent
route is to verify those partial_Z Phi traces directly. Cutoff motion
still has the parent's separate `[partial_Z,A1]` term, so it cannot be
discarded even after all seam jumps vanish.

Use the parent's unchanged zeroth-output physical example only:
`h=1/200`, `N=0`, `kappa=1/100000`, target t=302 and a=2.
The original strict count and deferred requests are

| Requirement | Order |
|---|---:|
| Nonstationary integrations p | 30,000,001 |
| Sufficient zero traces of weighted partial_Z Q Phi | through 29,999,999 |
| Sufficient identical weighted Q Phi traces for moving seams | through 30,000,000 |
| Parent raw Q Phi request | 60,000,006 |

The last nonzero trace at p-1 can be kept at the bulk order as in (6).
These trace indices are not a count of total original profile derivatives
and do not lower the parent's earlier source/root/A0/A1 request graph.
The normalized direct order-two estimate has p=0 and introduces NO new
nonstationary trace-jet request; its compact quotient/value trace and
original raw norm contracts still apply. No implication for higher
physical output orders is newly claimed by this zeroth-output ledger.

## 4. Actual gate and next obligation

For the actual source, establish on every original closed physical cover:

* matching characteristic value traces of Q Omega and Q K1, including
  every full-hull exterior trace;
* matching weighted partial_Z Q Phi traces at the requested finite jets,
  or the explicit bounds (8) for every surviving boundary contribution;
* all mixed chart, weight, fixed-profile and integer-cover transports,
  with the original band-uniform S-polynomial constants.

Compare traces AFTER the original physical transport. With the original
`J=((3,1),(1,5))`, determinant 14, `Pi_i f(Y)=f(J^i Y)` sends mode k
to `(J^i)^T k`. Factor the common covering first; the retained overlap
differences are at most three. Traces match when their coefficients at
each resulting PHYSICAL frequency match, not when their raw k lists
match. Haar preservation controls only the zero mode. The finite exact
transport controls retain this matrix and detect the raw-index shortcut;
they do not verify the actual infinite Fourier traces. If the
canonical assembled source has a proved smooth global compact extension,
these qualitative matching identities follow immediately; its quantitative
high-jet estimates and all original norm/domain conversions remain separate.
The present evidence does not materialize that assembled source or certify
these identities. `certify_actual_traces` refuses that unsupported upgrade.

Next verify this gate from the defining actual potentials, covariance and
lifted source on the closed covers, then validate their reduced high-jet
norms. No further generic signed cycle is used to assume those inputs.
The full regular order remains E=2-3kappa, with negative physical exponent;
complete absolute residual decay, infinite summation, forcing regularity
and complete velocity growth remain unproved.

Reproduction: `python test_audit.py -v`, then `python audit.py
--verify-record evidence.json --output FRESH_PATH.json`. Exact symbolic
checks and historical replays support these formulas, not a continuum
enclosure. Frozen historical files are bound without modification.
