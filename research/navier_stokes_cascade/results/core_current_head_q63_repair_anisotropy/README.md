# Q63: low angular repair order and explicit high radial bump bounds

Q62 bounds the axial target by leading velocity derivatives. The sufficient
isotropic order n+99998 can be replaced, in this principal-moment chain,
by the downward-closed profile index set

    i+j ≤ n+99998,       j ≤ n+2,                         (1)

where i counts X derivatives and j counts η derivatives of E₀,U₀.
This does not bound all physical wave derivatives at that restricted order.
Here it removes a high implicit-repair derivative input: the five repair
coefficients and their normalizing amplitude need only n+2 η derivatives,
while the high radial derivatives are explicitly bounded fixed bumps.
The pre-modulation profile and moving loop still require the corresponding
mixed inputs. No parameter or support is changed.

## Trace the restriction through the same leading comparison

Use Q57's two choices d=n+1+s, J=d+q; at most e=n+1 slots are Z,T
combined. In the frozen band X=R²/(2s_scale), η=Z s_scale^(-D).
A derivative landing on η must consume a Z or T slot; radial chart
changes never create one. Thus slow field compositions with d slots
request total profile order ≤d and η order ≤e.

The phase of Q60 uses first slow derivatives H_R,H_Z, and the matrix K
uses F_R,G_R. Consequently the normal and generator require total order
≤d+1 and η order ≤e+1. The path coordinate is independent of the slow
parameters in this fixed coefficient chart. Positive path derivatives
of the affine normal add no field derivative. The prescribed reference
rate is a function of path time and frozen label constants. Taking
the matrix products and reciprocals in the generator distributes
derivatives among subindices; it does not increase their Z,T count.
The fixed-path differentiated ODE has precisely the same property:
a forcing term uses a subindex on the generator and its complement
on the vector. Its inductive dependency set is downward closed.
Do not apply a total-order maximum over *all* slow directions here;
keep the multivariate version of the Q59 and Q60 recurrences.

The same reasoning applies to the fixed-v Haar column, differentiated
matrix solve and partition products. The original entrance uses the
fixed-label normal and its tangential inverse, not a slow derivative
of a newly chosen carrier. The low-order ambient propagator, floors
and comparison branch of Q61 remain its existing hypotheses.
For the target and four weights use Q62 and Q58: their explicit first
field/primitive derivatives cost at most one further η derivative.
Since d≤n+99997, (1) suffices. This includes the finalized pressure
memory in the target, via the terminal moments of Q62. The low-order
full-background correction of Q61 and the separate curl term retain
their own previously stated bounds.

## The exact repair has no radial coefficient derivatives

The original first reserved patch has y=log(X/ρp), with ρp frozen.
As specified in [stress realization, §5](../stress_realization_audit/README.md),

    ΔE=e(η) Σ_(three E bumps) c_l(η) β_l(y),
    ΔU=e(η) Σ_(two U bumps) c_l(η) β_l(y),
    β_l(y)=10 σ'(10(y−center_l)+1/2).                    (2)

The five closed supports are disjoint, have width 1/10 and lie inside
[0.45,2.55]⊂(0,5). They are not differentiated as η-dependent intervals.
Their coefficients solve the SAME accepted quadratic map
Bc+Q(c,c)=z_N(η), with constant B,Q after normalization. The normalization
contains e and e² and the transformed target contains λ⁻¹; all their
η derivatives remain in z_N. This is not a zero-pressure target.

Write βinv≥||B⁻¹||, qquad≥||Q|| and
||z_N^(k)||/k!≤Z_k/N, with the original accepted low-order root and
pointwise Jacobian bound ||J_c⁻¹||≤2βinv. Define

    D₀=2βinv Z₀,
    D_k=2βinv [Z_k+qquad Σ_(j=1)^(k−1) D_j D_(k−j)].      (3)

The [all-order repair theorem](../modulation_repair_all_order/README.md)
then gives ||c^(k)||/k!≤D_k/N, using N≥1.
Only k≤n+2 is needed in (3) for (1). If
e_k≥sup|e^(k)|/k!, put C_k=Σ_(j=0)^k e_j D_(k−j).
For each fixed bump, radial and η derivatives separate exactly:

    |∂y^r ∂η^k (e c_l β_l)| ≤ (k! C_k/N) ||β_l^(r)||∞.  (4)

Since the supports are disjoint the same maximum bound works for each
of the sums in (2), with flat endpoints included. No increasing count
of labels or factor of the modulation frequency is introduced here.
The modulation itself still has the N^(r−1) behavior proved earlier.

## A closed factorial bound for the actual bump

Use the actual squared-flat σ=g(x)/(g(x)+g(1−x)),
g(x)=exp(−1/x²), extended by flat endpoints.
The [cutoff certificate](../core_borel_cutoff_certificate/README.md)
gives g^(j)=P_j(1/x)g, P_(j+1)=−u²P_j'+2u³P_j.
Its degree is ≤3j and its coefficient absolute sum is ≤5^j j!,
because the next step multiplies this sum by at most 3j+2≤5(j+1).
The monomial bound u^a exp(−u²)≤ceil(a/2)! for u≥1 and
(2j)!≤4^j(j!)² give the raw derivative bound

    sup |g^(j)| ≤ 20^j(j!)³.                             (5)

It also holds at j=0. For H=g(x)+g(1−x), H≥1/81 and
|H^(j)|≤2·20^j(j!)³. Normalize reciprocal derivatives by (j!)³.
In Leibniz, binom(j,k)(k!)³((j−k)!)³/(j!)³=binom(j,k)^(-2)≤1.
The reciprocal recurrence is therefore dominated by
b₀=81, b_j=162 Σ_(k=1)^j 20^k b_(j−k).
Its exact positive solution for j≥1 is 81·162·20^j·163^(j−1),
so b_j≤81·3260^j. Another product convolution, with the same
binomial inequality and Σ_(k≥0)(20/3260)^k<2, yields

    ||σ^(j)||∞ ≤ 162·3260^j (j!)³,
    ||β_l^(r)||∞ ≤ B_r:=162·32600^(r+1)((r+1)!)³.        (6)

These conservative integer bounds use no numerical maximization,
Taylor convergence at an edge or high-order polynomial expansion.
They are uniform over the five translations of the original bump.

For X derivatives use the complete signed Stirling identity
∂X^r=X^(-r)Σ_(i=0)^r s(r,i)∂y^i. On the repair supports X≥a>0,
B_i≤B_r and Σ_i |s(r,i)|=r! (also at r=0). Thus

    |∂X^r ∂η^k Δ(E,U)| ≤ a^(-r) r! B_r k! C_k/N.         (7)

All constants are fixed-profile constants independent of the physical
band. In particular the radial part of the repair at r=n+99998 is
now explicit even without forming the huge integer B_r.
Equation (7) does not claim smallness uniform in r or change N with r.

The remaining leading-velocity task is to enclose the pre-modulation and
moving-loop mixed derivatives on (1), and the low η-order defect and
normalizer data Z_k,e_k for the same repair. Low physical value/inverse
gates must still be checked together. The radial pressure tail, uniform
actual-pulse and complete residual estimates and blow-up remain unproved.

Reproduce test_audit.py and audit.py --verify-record evidence.json
--output /path/to/fresh.json. Exact manufactured controls test the
jet counts, factorial majorant and separated repair derivatives.
