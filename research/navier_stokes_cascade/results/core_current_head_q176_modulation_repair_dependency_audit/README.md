# Q176: audit compact modulation and exact repair at the source level

Parent: `341e2842c6fa455fa3912797fc64bc47e393051f` (Q175).
Status: `Q176-compact-modulation-exact-repair-audited-with-refined-dependency-order`.

Q175 restored compact modulation and exact repair as the first independent
validation frontier. This checkpoint audits that route in the frozen source
records. It establishes a conditional same-map implication, verifies the
selected Jacobian and cone inequalities, traces all five moment rows and the
complete pressure memory, and corrects an overly coarse dependency edge in
`foundation_bundle_route`. Historical evidence is unchanged.

No field, root, pressure gauge, width, support or schedule is changed. The
retained selected-profile continuation and low-order regional suprema are not
independently recomputed here, so this is not an unconditional construction,
pulse estimate, residual improvement or blow-up result.

## 1. The frozen records have a finer dependency order

The coarse foundation graph places its `matching_roots` node before
`modulation_repair`. That is a useful route summary, but its named provider
`matching_root_foundation` is not an upstream input to the frozen all-order
map theorem. The actual records and their source hashes give the following
order:

1. `prescribed_continuation_foundation`, `loop_modulation_bounds` and
   `repair_state_bounds` supply the retained smooth input and low-order maps.
2. `regional_cone_join` accepts the finite-frequency leading-value map from
   those low-order inputs.
3. `modulation_repair_all_order` freezes the regional record and proves that
   this same accepted map has every fixed derivative.
4. `annular_exterior_foundation` freezes both records, supplies the whole
   compact smooth input, and composes the conditional all-order map.
5. `matching_root_foundation` is a later quantitative reaudit: it freezes the
   modulation and annular records and sharpens their inverse estimates.

Thus the low-order repair acceptance and the later matching-root reaudit must
be represented by distinct nodes. With that split, all twelve audited edges
point forward and there is no theorem cycle. The dates in the live history
are consistent with the frozen source bindings: regional (2026-10-16),
all-order modulation (2026-10-17 04:09), annular/exterior (05:11), and the
matching-root reaudit (19:28).

## 2. Same finite frequency and same accepted repair root

The signed variance equation uses its smooth removable zero rather than an
unsigned square root. Its derivative in the signed root is positive. The
phase map has positive angular derivative, so the moving inverse is smooth.
For every fixed radial order `r`, the frozen calculation distinguishes

\[
  \|D_X^r(\hbox{field error})\|=O_r(N^{r-1}),\qquad
  \|D_X^r(\hbox{moment error})\|=O_r(N^{\max(r-2,-1)}).
\]

The estimates are finite at every fixed order; high radial field derivatives
are not uniformly small in `N`.

For the exact repair, using the same original normalizer makes the pointwise
map

\[
  Bc+Q(c,c)=z_N(\eta),\qquad J_c=B+2Q(c,\cdot).
\]

The all-order provider retains `beta<=1000`, `q<=10^5`,
`D=H^16`, `N>H^32` and `||c||<=2 beta D/N`. Therefore

\[
 \|B^{-1}(J_c-B)\|\le 2\beta q\|c\|
 \le {4\beta^2qD\over N}
 < {4\cdot1000^2\cdot10^5\over H^{16}}<\tfrac12
\]

already for `H>=8`, and `||J_c^-1||<=2000`. Differentiating this same equation
then gives, for normalized fixed-order jets,

\[
 d_k\le2\beta\left(z_k+q\sum_{j=1}^{k-1}d_jd_{k-j}\right).
\]

The highest linear terms are inside `J_c`. No contraction in a new `C^k`
ball and no frequency reselection is used. The later matching-root record
sharpens the same calculation to `beta<=800`, `q<=25000` and inverse at most
`1600`; Q176 keeps the conservative theorem and the later refinement
separate.

## 3. Complete pressure memory and the compact cone join

The compact transformed defect rows are

\[
 T_c(M,R,J,S,C_p)=\left(M,{M-J\over\lambda},R,-S,C_p\right),
 \qquad \lambda>0,
\]

with exact inverse

\[
 (r_1,r_2,r_3,r_4,r_5)\mapsto
 (r_1,r_3,r_1-\lambda r_2,-r_4,r_5).
\]

Zero transformed rows therefore restore all five original defects as
functions of `eta`. The fifth row preserves the axis datum and accumulated
pressure; it cannot be inferred from the other four. The annular/exterior,
compact/global and core-pressure records retain the same root, gauge,
positive-order `Pi1/Pi2` memory and backward pressure identity.

For `R=4G^3`, the state perturbation loses at most
`256 R^10 G^-64`. Absorption into half of the `G^-4` gap reduces to

\[
 512\,4^{10}<G^{30},
\]

which already holds at `G=4` and hence at the selected much larger value.
This joins the inherited inner collar, compact modulation/repair region and
outer leading-value estimates conditionally. It does not supply uniform
physical-pulse derivatives or the full viscous/cutoff/interaction residual.

## 4. Remaining obstruction and next obligation

The source records still retain the actual selected-profile finite-order
continuation, loop/repair inverse and regional gap suprema. Exact recurrences,
manufactured jets and passing tests verify algebra and evidence integrity;
they do not independently evaluate those continuum inputs. The precise next
obligation is to validate those selected-profile suprema, then audit the
weighted first-stress edge quotient and the reserved coefficient/mean
supports on the same roots. Only after those inputs close can the all-scale
actual-pulse, complete covariance/mean and full residual estimates be
advanced.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q176_modulation_repair_dependency_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q176_modulation_repair_dependency_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q176_modulation_repair_dependency_audit/evidence.json \
  --output fresh.json
```
