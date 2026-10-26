# Conditional finite-stage compact request recurrence from the corrected seed

Parent `f27853fe050706ac1834d3d5892cd66526cfcce8` synchronizes the
corrected post-signed induction seed but leaves the later stage-indexed
compact graph open.  This checkpoint supplies that finite recurrence.  It
uses the actual integer request traces of the existing generic-signed and
sourced mean modules and the older half-line order recurrence.

The conclusion remains conditional on every named all-fixed-order and
all-fixed-pole background, root, inverse, support, pressure and compact-map
contract.  No large derivative, factorial, operator constant or actual
selected pulse is evaluated.  Tests check exact arithmetic and frozen
provenance; they do not prove those continuum contracts.

## 1. Corrected stage rates

Let

    mu=1-kappa,
    B_j=3/2-2 kappa+j delta,
    C_j=B_j+1/2,
    delta=1/2-3 kappa>0.                                 (1)

The corrected post-signed seed is `j=0`.  The inherited half-line
certificate applies because its old minimum is `3/2-6 kappa`, while
`B_0` is larger by `4 kappa`.

One repeated cycle consists, forward in time, of the sourced wave and
fresh means followed by the generic signed wave and fresh means.  The
radial compact source keeps the coarse order `mu`.  Define the normalized
integration counts

    p0s(j)=strict(mu,C_j),       p1s(j)=strict(B_(j+1),C_j),
    p0g(j)=strict(mu,C_(j+1)),   p1g(j)=strict(B_(j+1),C_(j+1)),       (2)

where `strict(a,b)=max(1,floor((b-a)/kappa)+1)`.

The actual sourced `request_trace` reaches the preceding state after

    d_s(j)=2 p0s(j)+18,                                       (3)

and the actual generic `dependency_trace` reaches the sourced state after

    d_g(j)=2 p0g(j)+17.                                       (4)

The `A0` branches dominate the corresponding `A1` requests; the code
still computes all four counts and both complete traces.  At `j=0`,

    (p0s,p1s)=(100000,4),       d_s=200018,
    (p0g,p1g)=(149997,50001),   d_g=300011.                  (5)

These offsets include complete covariance, pressure/scalars, actual axial
`A1` velocity, fresh fast/slow means and inverse-normal input.  Constants
and harmonic sets may grow with stage, but the derivative graph remains
finite.

## 2. Backward derivative recurrence

Let `R_j(n)` be a safe first-real primitive order sufficient for every
stage-`j` state quantity through order `n`.  The corrected seed gives

    R_0(n)=n+500038.                                         (6)

Because each backward cycle decreases the stage index exactly once,

    R_(j+1)(n)=R_j(n+d_g(j)+d_s(j)).                         (7)

Thus

    R_j(n)=n+500038+sum_(i<j)(d_g(i)+d_s(i)),                (8)

which is finite for every fixed `j,n`.  For example,

    R_0(0)=500038,       R_1(0)=1000067.                     (9)

There is no use of the desired current `A0` output on an earlier edge.
The recursion terminates at the already closed first-real primitive graph.

For physical jet `m` and desired `q^N`, put

    T=(N+3/2+2h+2m)/h,
    p0=max(1,floor((T-mu)/kappa)+1),
    p1(j)=max(1,floor((T-B_j)/kappa)+1).                    (10)

The current generic trace first reaches its sourced predecessor; (3),
then (7), expands the remaining graph.  At `h=1/200,m=N=0`,

| post-signed stage | raw `A0` | raw `A1` | first-real ceiling |
|---:|---:|---:|---:|
| 0 | 60,200,007 | 60,100,009 | 60,500,042 |
| 1 | 60,200,007 | 60,000,015 | 60,900,077 |

Both branches have strict physical exponent margin `1/20000000` under
the retained integer rule.  The ceiling is sufficient, not a necessary
lower bound for every possible proof.

## 3. Pole, operator, support and harmonic archive

For every derivative request `L` and original pole `p`, retain the
existing sufficient common-pole rule

    q=max(q_previous,6p,6L+9,6L+12p+6).                    (11)

At every distinct `q`, archive the separate `6q/12q` pressure/fast-map
shift and the oriented `6q` primitive constant.  The common estimate
argument contains the SUM of all archived constants.  No operator constant
at a larger pole is assumed to dominate one at a smaller pole.

The old finite-stage support induction applies to each node: derivatives
preserve the fixed closed supports, cut potentials and pressures have
smooth endpoint extensions, distinct labels retain closed separation, and
mean fields may fill the auxiliary torus while every nonzero interaction
retains a wave factor.  These are inherited analytical contracts.  The
angular set remains finite with the conservative bound

    M_j<=4*2^j.                                             (12)

Equations (8), (11) and (12) give a finite derivative, pole, operator,
support and harmonic ledger for each fixed stage.  They do not give
uniform constants as `j` tends to infinity.

## 4. Conditional physical consequence and remaining boundary

Under the primitive contracts, induction on (7) gives the actual-source
inequalities required by the parent checkpoint at every fixed `j,m,N`.
The same-domain polynomial absorption then makes the same-state `A0/A1`
tails `O(q^N)`.  Combining them with the inherited regular rate certificate
gives, for each fixed finite stage and physical jet,

    |R(U^[j])|_m <= C_(j,m) q^[h B_j-(3/2+2h+2m+eta)]       (13)

after a chosen positive polynomial reserve `eta`.  Constants may depend
arbitrarily on the fixed stage and derivative order.

For `h=1/200`, `m=N=0`, `eta=1`, the corrected recurrence first gives a
nonnegative exponent at `j=1002`.  This is an exact existence threshold,
not construction or numerical evaluation of 1002 fields.

The quantifiers remain

    for every fixed target there is a finite conditional state,

which does not construct one field meeting all targets.  A fixed-`q`
contraction, stage-uniform constants, convergence of the uncut series and
the substitution `j=j(q)` are not proved.  The next obligation is to join
this corrected current-seed hierarchy to a complete-potential diagonal,
then audit all cutoff curls, pressure gradients, viscosity, nonlinear
interactions, endpoint/localization forcing, energy and the core lower
bound.  No unconditional residual improvement or blow-up is asserted.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
