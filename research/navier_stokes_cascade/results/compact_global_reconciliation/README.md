# Explicit-cutoff compact repair and global pressure reconciliation

Parent: `b295e7c7179e638f4f2122b5477a8edbab7d7d03`.
Status: `explicit-cutoff-to-same-root-compact-cone-and-global-pressure-reconciliation`.

The parent fixes an explicit source-admissible variance cutoff and proves the
two raw bounds that were missing from the compact construction. This checkpoint
follows those bounds through the already derived conditional interfaces,
retaining their quantifiers and pressure memory. Relative to the retained
core/native/regional/root premises, the **same** frequency and five-row repair
now satisfy the complete compact fourth pressure/flux budget, preserve and
join the inner/compact/outer leading value cones, and restore the five global
moment defects so the forward and backward exterior stresses agree.

This is an analytical implication between frozen checkpoints, not a new
numerical profile or an independent verification of every inherited premise.
It is not a complete all-scale pulse/residual theorem or an unconditional
blow-up proof. No parameter, width `C^-20`, support, gauge, root, `K=H^64` or
`N=1+floor(H^32)` is changed. Historical evidence is untouched.

## 1. Both loop inputs now hold for the selected cutoff

Put `A=C^4096`, `G=exp(A^32)`, `H=exp(G^4)`. The parent proves for its
explicit flat cutoff

`Z4<2,624,400`, `Z5<213,101,280`,

and verifies the fourth target gate and the fifth gate after combining the
latter with the native estimate `R5<A^32`. The fourth gate implies the mixed
third gate, so the three-variable moving-loop bridge gives full joint `J3`
primitives below `H`. Independently, the joint-fourth checkpoint gives
`Gamma<2H` and the compact field bound `J4<K^2=H^128`. These implications
use the explicitly selected member of the cutoff class; they do not revive
the corrected claim that equation (C.9) had already fixed a shape.

The fifth gate gives

`log Lambda5 < G^2`.                                      (1)

The complete compact-map interface asks for

`Lambda5 <= N H^260/A^1044`.                              (2)

Since `N>H^32`, the logarithm of the right side of (2) is greater than
`292G^4-1044 log A`. Using `log A<G`, it is enough that

`G^2+1044G < 292G^4`,

which already holds at `G=2` and strengthens thereafter. Thus (1) implies
(2). The interface's full moment bound is below `H^261`, and its complete
ordinary fourth jets of `F,U,Pi,V,p1,p2` are below `H^925<H^1024`. Both
angular pressure terms, axial velocity, axis pressure, accumulated `Cp`, and
incoming compact memory are present. This closes the compact fourth-map
budget; it does not manufacture a global all-order leading-profile estimate.

## 2. Same frequency, same root, and the regional cone join

The mixed bridge differentiates the modulation before evaluating the fast
phase, hence retains

`a_N=a_L-2 Acal_y/N`,

`b_N=exp(-Acal/N)(b_L+2 Bcal_y/(N E))`.                  (3)

At the unchanged `N`, modulation plus the same exact repair has state error
less than `G^-64/2`. With the original loop radius `R=4G^3`, the cone loss is
bounded by `256R^10 G^-64`. It is below `G^-4/2`, since after clearing
positive powers it is enough that `512*4^10<G^30`; the actual `G` is much
larger than four. The repaired compact gaps are therefore at least `G^-4/2`,
and `a>=G^-3/2`.

The repair is the same previously accepted `C1` root, identified by pointwise
uniqueness; no high-norm contraction or new root is introduced. Its endpoint
remains exactly `Xrep=rhop exp(5)`. The regional ledger uses this same
frequency, root and endpoint. Its ten-region value estimates therefore join
the untouched inner collar, compact modulation/repair region, first outer
ramp, axial transition and compensated exterior. The corrected signed
outer-ramp comparison is the one used here.

This is a leading value/cone join in the similarity profile. It is not a
claim that every physical pulse derivative or every positive-order
coefficient has acquired a uniform all-scale bound.

## 3. Exact five-row memory and exterior stresses

Write the original defect vector as `d=(M,R,J,S,Cp)`, where `R` is the
renormalized heat defect. The annular and compact row maps are

`Ta d=(M,J-4 eta R,R,-(S-8 eta M),Cp)`,                 (4)

`Tc d=(M,(M-J)/lambda,R,-S,Cp)`.                        (5)

Their exact inverses are

`(r1,r2,r3,r4,r5) -> (r1,r3,r2+4 eta r3,8 eta r1-r4,r5)`,

`(r1,r2,r3,r4,r5) -> (r1,r3,r1-lambda r2,-r4,r5)`.      (6)

The second identity uses the unchanged `lambda>0`. Thus zero transformed
targets are precisely the five original zero defects. Because the future
heat tail is identical on both sides of the compact join, finite endpoint
heat equality is the global renormalized equality. All identities hold as
functions of `eta`, not only at sampled angles.

The original axis datum is `Pi_axis=-Cp_id`. Every edit either preserves the
five totals or is paired with the same exact five-row repair, so

`M_infty=R=J_infty=S_infty=0`, `Pi_axis+Cp_infty=0`.      (7)

Hence the completed leading pressure is the backward tail

`Pi(X)=-integral_X^infty E^2/(2x) dx`.                   (8)

Substitution of (7) in the complete exterior difference formulas makes
`Delta V0`, `Delta(r^2 T_theta)` and `Delta(r T_z)` zero. This joins the
signed forward stress to the backward stress without resetting pressure.
The fifth row is essential: with all other defects zero and `p(eta)=eta`,
the axial pressure contribution at `eta=0` is `-X`, not zero. Smallness or
vanishing at one angle cannot replace exact `Cp` restoration.

## 4. Consequence and remaining limits

The compact cutoff premise in the current global pressure/moment foundation
is now discharged. Its finite positive-order coefficient consequence can be
reached with the same leading state: every finite coefficient uses the exact
linear five-row solve, has the documented common velocity support, and leaves
the exceptional first two pressure tails on the reserved mean patch. This
does not discharge those checkpoints' retained exact-core, native/regional
existence or root premises, and gives no norm uniform in coefficient order.

The next obligation is to audit the remaining higher-coefficient and
all-scale pulse-family defining/existence/inverse/gap/value/derivative/support
estimates against this reconciled leading state. Complete-wave covariance
means must retain radial and axial velocity and both exceptional pressure
tails. The complete residual must still include viscosity, every cutoff and
all interactions, followed by the distinct realization/comparison steps.

The diagnostic suite checks (1)--(6) with exact integer/rational arithmetic,
binds the eight frozen dependency records, and includes the nonzero pressure-row
control. These checks verify bookkeeping and source immutability; they are
not a continuum proof of the inherited analytical premises.

Reproduce with:

```sh
python test_audit.py -v
python audit.py --verify-record evidence.json --output fresh.json
```

The output path must be fresh; frozen evidence is never overwritten.
