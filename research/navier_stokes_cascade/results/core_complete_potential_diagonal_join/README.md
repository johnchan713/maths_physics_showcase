# Corrected complete-potential diagonal from the synchronized finite stages

Parent `85272d0234a1b4da33b315b19fbd974b339f5f96` supplies the
stage-indexed compact request, pole and support recurrence for the corrected
post-signed seed.  This checkpoint joins those finite states to the existing
complete-potential cutoff mechanism.  It also corrects the last paragraph of
the parent: stage-uniform constants, fixed-`q` contraction and convergence of
the raw correction series are not necessary for this implication.

The result remains conditional on every named all-fixed-order background,
root, inverse, support, pressure, compact-map and physical-transfer contract.
No large derivative, actual pulse coefficient or cutoff scale is evaluated.
Passing exact controls and frozen provenance do not prove those continuum
contracts.

## 1. Corrected finite states and complete representatives

Keep

    kappa=1/100000,       delta=1/2-3 kappa,
    B_j=3/2-2 kappa+j delta,       C_j=B_j+1/2.             (1)

The parent synchronizes every fixed-stage compact tail with the same state
used by the older complete finite-stage induction.  Relative to the retained
analytic contracts, that induction gives real divergence-free states
`U^[j]=(u^[j],p^[j])` on one fixed `q` domain and

    |R(U^[j])|_m <= C_(j,m) q^[h B_j-K_m],
    K_m=2A+3/2+2m,       A=1/2+h.                            (2)

Equation (2) is the complete fixed-stage residual: pressure, viscosity,
cutoff-curl terms, both old/new crosses, new squares, Gaussian tails and all
mean components are included.  Constants may depend arbitrarily on the
finite stage and jet.

For the actual cycle `j-1 -> j`, retain the operator-selected gauge and write

    Delta u_j=curl(A_j)+b_j e_theta,       partial_theta b_j=0,
    pi_j=p^[j]-p^[j-1].                                  (3)

`A_j` contains the complete sourced and signed wave potentials and the fast
and slow azimuthal mean potentials.  The direct mean swirl stays outside the
curl.  `pi_j` telescopes every intermediate pressure change; complete
covariance memory, the axial `A1` velocity, moving axial pressure moment and
both exceptional same-base pressure tails are not discarded.

A conservative normalized ledger, with `B=B_(j-1)` and `C=B+1/2`, is

| Actual representative | Lower order |
|---|---:|
| Complete wave velocity | `B-kappa` |
| Wave potential | `C-kappa` |
| Mean potential, direct swirl, fixed-gauge pressure difference | `C-2 kappa` |
| Mean radial velocity | `C+1-2 kappa` |

Every row has exact slope `delta` in `j`.  The weakest row exceeds
`j delta` by exactly one.  The fixed closed annular supports, smooth zero
extensions, common chart identifications and uniformity over every physical
band at a fixed stage are inherited from the synchronized support and
finite-stage contracts.  The harmonic set is finite at every stage, although
its constants may grow.

Physical normalization and phase differentiation therefore give, for the
tuple `Z_j=(A_j,b_j e_theta,pi_j)`,

    |Z_j|_m <= C_(j,m)(1+|log q|)^P_(j,m)
               q^[h delta j-ell_m],       ell_m=2A+2m.     (4)

The constants in (4) need not be uniform in `j`.

## 2. Cut complete potentials, direct swirl and pressure

Fix a smooth `chi` equal to one on `s<=1/2` and zero on `s>=1`.  For
`chi_a=chi(aq)`, define

    V_(j,a)=curl(chi_a A_j)+chi_a b_j e_theta,
    P_(j,a)=chi_a pi_j.                                   (5)

This is the required complete cut.  Cutting only `Delta u_j` would generally
destroy divergence.  In particular,

    V_(j,a)=chi_a Delta u_j+grad(chi_a) cross A_j,          (6)

so every wave commutator and the radial mean term
`-a chi'(aq) q_z Psi_j` remains.  Direct swirl is divergence-safe because it
and `q` are axisymmetric.  Pressure cutoff terms `pi_j grad(chi_a)` also
remain.

For the implicit coordinate `q=1-t+z^2 q^(2h)`, exact weighted
differentiation gives

    |partial_z^r partial_t^s chi(aq)| <= C_(r,s) q^(-Dr-s),
    D=1/2-h,                                               (7)

with a constant independent of `a` on the transition interval.  Combining
(4)--(7), reserving one curl, yields

    |(V_(j,a),P_(j,a))|_m
       <= C~_(j,m)(1+|log q|)^P~_(j,m)
          q^[h delta j-L_m],       L_m=ell_(m+1).          (8)

The fixed potential and pressure gauges are essential.  An arbitrary
stage-dependent potential gradient produces
`grad(chi_a) cross grad(phi_j)`, and an independent time pressure gauge
produces `c_j(t) grad(chi_a)`.  Equal uncut states alone therefore do not
imply (8).

## 3. One locally finite diagonal absorbs arbitrary constants

At stage `j`, collect the finitely many constants in (8) for jets through
`j`.  Put `gamma_j=h delta j/2`.  If `c>=log C` and `P>=0`, choose
`s_j=log a_j` strictly above

    max{s_(j-1)+1, log(qstar^-1)+1, 1,
        2(c+j)/gamma_j+P^2/gamma_j^2}                       (9)

for every current constraint.  The inequalities
`log(1+s)<=sqrt(s)` and
`P sqrt(s)<=gamma_j s/2+P^2/(2 gamma_j)` show on the whole interval
`0<q<=1/a_j` that

    C(1+|log q|)^P q^(gamma_j) <= 2^(-j).                  (10)

Such an `a_j` exists for every finite list, regardless of how fast its
constants grew.  Also `a_j>2a_(j-1)`.  Define

    u*=u^[0]+sum_(j>=1) V_(j,a_j),
    p*=p^[0]+sum_(j>=1) P_(j,a_j).                          (11)

On every compact subset of `q>0`, only finitely many summands are active.
Thus (11) is one smooth real divergence-free local field.  No convergence
of the raw series and no fixed-`q` contraction are invoked.  Multiplication
by `chi(a_jq)` and the local curl do not enlarge the inherited closed radial
supports.  Hence every correction still vanishes on the protected inner
region and in the heat exterior: the inherited inner velocity and exterior
field are unchanged, conditionally on those support contracts.

For fixed `J>=max(1,m)` and `q<=1/(2a_J)`, the first `J` cutoffs equal one
and telescope in the fixed gauges.  The remaining locally finite tail obeys

    |(u*,p*)-U^[J]|_m
       <= 2^(-J) q^[h delta (J+1)/2-L_m].                  (12)

The comparison stage is fixed after choosing the requested jet and power;
there is no uncontrolled substitution `j=j(q)`.

## 4. Complete residual comparison

Set `(e,pi)=(u*,p*)-U^[J]`.  The exact difference is

    e_t-Delta e+grad pi+(u^[J].grad)e
      +(e.grad)u^[J]+(e.grad)e.                            (13)

This lists time, viscosity, pressure, both nonlinear crosses and the new
square.  The complete cut curls are already inside `e`.  With
`H_m=L_m+1`, choose one finite `J>=m+2` satisfying

    h delta (J+1)/2 >= N+1+L_(m+2)+H_(m+2),
    h B_J-K_m >= N+1.                                     (14)

Both left sides have positive slope.  Equations (2), (12)--(14) imply for
the same field (11)

    for every fixed m,N,       |R(u*,p*)|_m=O(q^N).        (15)

At `h=1/200,m=N=0`, the deliberately coarse sufficient comparison stage is
`J=12816`; the finite-state residual condition alone needs only `J=1402`.
These are exact exponent budgets, not evaluations of thousands of fields.

This proves a **conditional local diagonal implication**.  It corrects the
parent's suggestion that stage-uniform constants were the obstruction.
The actual obstruction to an unconditional theorem remains the independent
proof of the primitive background, root, inverse, support, pressure,
compact-map and physical-transfer contracts.  Closed-time one-sided jets,
localization and initialization, an admissible smooth compact force,
whole-space energy and the global comparison remain open.  The local inner
velocity is protected as above, but turning that conditional local fact into
the final global lower-bound/blow-up argument is not done.  No unconditional pulse bound,
global complete residual theorem or blow-up is asserted.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
