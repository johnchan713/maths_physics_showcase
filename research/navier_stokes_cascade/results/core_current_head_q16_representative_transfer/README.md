# Q16: termwise physical transfer of the complete representative ledger

Parent: `0efd95500b4291224821036370c63acd78c4a3a7` (Q15).
Status: `conditional-complete-representative-physical-transfer`.

This checkpoint connects the normalized representative orders to the physical
bound required by Q13. It keeps the corrected delta=1/2-3*kappa,
kappa=1/100000, and root h=exp[-8(exp(64)+10)]. The older canonical
realization uses a different slope; its numerical stage rates are not
imported. Its decomposition and physical-transfer mechanism are used only
with the explicitly corrected orders below.

## 1. Actual components and normalized inputs

At the cycle j-1 to j, write B=B0+(j-1)delta, B0=3/2-2*kappa.
The complete diagonal provider's conservative normalized ledger is:

| Component | Normalized order alpha | alpha-delta*j | Physical normalization |
|---|---:|---:|---:|
| Complete wave velocity (including curl error) | B-kappa | 1 | Q^-A |
| Sourced and signed wave potentials | B+1/2-kappa | 3/2 | Q^-h |
| Fast/slow azimuthal mean potentials | B+1/2-2*kappa | 3/2-kappa | Q^-h |
| Direct axisymmetric swirl | B+1/2-2*kappa | 3/2-kappa | Q^-A |
| Complete fixed-gauge pressure difference | B+1/2-2*kappa | 3/2-kappa | Q^-2A |
| Mean radial velocity | B+3/2-2*kappa | 5/2-kappa | Q^-A |

Here A=1/2+h. The actual representative is
A_j=A_j^wave+Psi_j e_theta, B_j^sw=b_j e_theta,
pi_j=p^[j]-p^[j-1], with Delta u_j=curl A_j+B_j^sw.
The last and first rows audit velocity consistency; they cannot substitute
for the potential rows. The pressure row is the telescoping difference of
all intermediate reconstructions in a single fixed gauge. Complete
covariance, gd-A1 gd, the derivative of J_z+c_rho P, and both same-base
pressure memories remain. An unweighted old pressure term cannot be
assigned a flat-edge weight merely because a new wave has one.

## 2. Physical graph and phase transfer

The physical-pulse provider uses Q=2^-ell, epsilon=Q^h, S=ell^2, and

    partial_r^phys=Q^(-1/2)(partial_R+M_i d_r R^(d_r-1)L_i),
    partial_z=Q^-D partial_Z,
    partial_t^phys=Q^(-1-h)(-epsilon partial_T+c_i N_i).

The graph term remains for amplitudes. It cancels for the phase through
L_i Phi=0. Its displayed primary/forced phase costs are 1/2+h/2 in space
and 1+3h/2 in time. The older full-representative transfer uses the more
conservative worst cost

    nu=1+5h/2 < 41/40 < 2.

Using nu for every derivative safely covers both dictionaries, including
cylindrical frame and graph differentiation, under their coefficient/phase
hypotheses. Harmonic ranges are finite at a fixed stage; their size affects
constants and polynomial degrees, not the q exponents. The discrete carriers,
cover and label parameters are frozen before differentiation.

For a normalized row of order alpha and physical normalization Q^-a,
the conditional total-m-jet bound therefore has power

    E = h*alpha-a-m*nu.

Subtract the target E_target=h*delta*j-2A-2m. Exactly,

    E-E_target=h*(alpha-delta*j)+(2A-a)+m*(1-5h/2).       (1)

All six rows have positive margin. The smallest is the complete pressure row:

    E-E_target >= h*(3/2-kappa)+m*(1-5h/2)>0.             (2)

For comparison, the wave-velocity margin minus the pressure margin is
1/2+h/2+h*kappa>0; all potential and remaining velocity rows are larger.
This establishes that the requested coarse physical exponent is compatible
with every component, including pressure, at the exact root parameter.
There is no invocation of h=1/200 and no uniform positive margin as h->0.

## 3. From each band to one physical seminorm

On 1/4<=q/Q<=4, every real E satisfies

    Q^E <= 4^ceil(|E|) q^E.

The factor is finite after fixing j,m. For ell>=16 the physical-tail
provider gives S<=16|log q|^2. If P(S+K)=sum c_i(S+K)^i with c_i>=0 and
fixed K>=0, then for d=degree P,

    P(S+K) <= [sum c_i(16+K)^i] (1+|log q|)^(2d).         (3)

Thus no S>=K assumption or derivative-dependent starting band is needed.
Polynomial degrees remain finite but may depend on stage and order. Keep
these log powers; do not spend (2) or a compact/residual reserve twice.

Wave amplitude estimates use sqrt(zeta), whereas compact mean estimates
may use zeta. Q15 applies directly to the latter. For the former its same
proof with c=w^2/2 gives

    sqrt(zeta) delta^(-p)
       <= n! 2^n w^(-2n) delta^s,
       n=ceil((p+s)/2), 0<w<=1.                          (4)

The interior delta=1 is again handled separately. The fixed width is
C_setup^-20; no new width is chosen. Equations (3)--(4) preserve scale
uniformity if the input coefficient bounds were scale-uniform. The Gaussian
amplitude factor is at most one and may be discarded for this size bound;
its decay remains available separately for temporal cutoff tails.

The wave support provider states at most one nonzero actual wave label at
a point and common-chart agreement on closed supports. A fixed stage has
finitely many sourced/signed substeps. Mean torus fields need their own
fixed-stage overlap/support bounds; wave-label separation must not be used
as if means had disjoint labels. With those overlap bounds, summing the
finite component estimates proves

    |D^m(A_j,b_j e_theta,pi_j)|
       <= C_(j,m)(1+|log q|)^P_(j,m) q^[h*delta*j-(2A+2m)]. (5)

Smooth extensions and the one-sided tau=0 domain must be the same as in the
input contracts. Nothing here extends the field through q=0.

## 4. Exact remaining frontier

Equation (5) is a conditional transfer, not an independently established
actual-representative theorem. Each row still requires a band-uniform
normalized coefficient/phase estimate at Q14's expanded primitive ceilings
and archived poles, for the actual same-choice sourced/signed and mean
formulas. Q15 removes flat-edge poles only when the appropriate weighted
estimate exists. Pressure terms lacking that weight require their own
unweighted bound and same-base cancellation/support identity.

Next: verify these row-specific input estimates and the mean/pressure support
and overlap statements from the actual current formulas, rather than infer
them from finite graph size or rate tables. This identifies the remaining
boundary precisely: the physical exponents and log/edge conversion now close,
but the actual normalized coefficient and support hypotheses are not proved
here. Q13 supplies cutoff uniformity once (5) is available; endpoint closure,
unconditional residual improvement and blow-up remain unproved.

## 5. Validation

Exact checks compare all six rows with the live corrected diagonal provider,
verify the smallest margin and polynomial-log conversion, and keep the root
symbolic. Q15 is regenerated and the full inherited evidence graph rehashed.
Historical frozen evidence and unrelated files remain untouched. On
30 October 2026, Q15's focused CI had succeeded and its broad run was still
in progress at inspection. Passing tests is not a continuum proof.

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q16_representative_transfer/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q16_representative_transfer/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q16_representative_transfer/evidence.json --output q16_fresh.json
```
