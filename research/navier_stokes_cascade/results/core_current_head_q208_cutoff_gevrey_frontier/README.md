# Q208: existing Borel-cutoff provider join and narrowed primitive frontier

Parent: Q207, `929c5f5390e183500f0e0e030e2b0d615f11db68`.
The live branch, Q207, its active CI, Q204--Q205, Q120, Q192, the
weighted-stress foundation and the regional gate were read first. Historical
evidence, the selected profile, pressure gauge, roots, supports, cutoff
sequence and unchanged `C^-20` activation width remain fixed.

Q207 leaves a direct stage-602 primitive-seminorm, gap and inverse bundle at
derivative ceiling

\[
 R=36\,479\,350\,300.                              \tag{1}
\]

Its discussion treats the fixed cutoff shape as an unevaluated finite
constant. The live dependency closure contains Q120's all-order Euler bound
for the exact implemented squared-flat Borel cutoff, and Q204--Q205 choose
each background or stage scale only after the relevant finite constants are
known. This checkpoint computes the corresponding high-order cutoff bit
budget and identifies the exact missing premise for applying it: Q120 leaves
actual-field identity uncertified, while Q207 does not explicitly identify
its qualitative `chi` with Q120's function.

This is a dependency correction, not a new cutoff theorem. It also records a
discarded draft error: replacing the squared-flat seed by `exp(-1/x)` would
change the construction and repeat an obligation already completed in Q120.

## 1. The exact existing cutoff

Q120's fixed function is

\[
 g(x)=\begin{cases}e^{-1/x^2},&x>0,\\0,&x\le0,\end{cases}
 \qquad
 \sigma(x)={g(x)\over g(x)+g(1-x)},
 \qquad
 \chi(s)=1-\sigma(2s-1).                            \tag{2}
\]

It proves, analytically and without finite sampling, the Euler-seminorm
bound

\[
 K_m:=\max_{0\le j\le m}\|(s\partial_s)^j\chi\|_\infty
 \le 82\,4097^m(m!)^3.                             \tag{3}
\]

The proof is already frozen in Q120: the Gaussian-polynomial recurrence
gives `12^j(j!)^3`, the reciprocal recurrence gives
`81*2048^j(j!)^3`, Leibniz gives the quotient bound, and the Stirling
expansion converts ordinary derivatives to Euler derivatives. Q208 does not
reprove or replace that result.

For `m>=3`, Q120's integer-only estimate is

\[
 \lceil\log_2 K_m\rceil
 \le b_m:=7+13m+3m\lceil\log_2m\rceil.             \tag{4}
\]

At (1), `2^35<R<2^36`, so

\[
 b_R=7+(13+3\cdot36)R
     =4\,414\,001\,386\,307.                       \tag{5}
\]

If Q207's cutoff is explicitly bound to (2), its complete stage-602 cutoff
shape can be carried symbolically as `K_R<=2^4414001386307`. No factorial,
derivative array or integer of that bit length is materialized.

## 2. Same-object and schedule join

The frozen dependency graph contains Q120 inside Q207's closure. A dependency
path proves provenance, but it does not by itself turn Q120's explicit false
`actual_field_identity_certified` flag into a same-function theorem. The
closure also contains the current diagonal records:

1. Q204 selects each background Borel scale after all finite coefficient
   seminorms at that order are known.
2. Q205 selects each stage cutoff after the finite stage constants and log
   powers are known.
3. Q207 retains the potential, direct swirl and pressure under the same cut,
   with all six residual difference groups.

After the missing same-cutoff identity is supplied, no pre-existing numerical
scale list is needed to absorb (5). For each fixed row, (5) is one of the
finite constants supplied before the rule-defined scale is selected. This
does not produce a uniform constant in the stage or an uncut convergent
hierarchy.

The correction is specific to the Borel cutoff (2). The compact-X and
variance cutoffs have different identities and keep their own inherited
bounds. A qualitative statement such as “some smooth cutoff” would not be
enough. The precise first cutoff obligation is therefore

\[
 \chi_{Q207}=\chi_{Q120}\quad\hbox{as the actual field cutoff}. \tag{6}
\]

No frozen boolean or exact identity in the reviewed Q207-to-Q120 chain
certifies (6). Q204's prose calls the schedule “implemented,” but does not
repair Q120's deliberately false field-identity flag. Q208 therefore does
not promote (3) to the actual Q207 field.

## 3. What remains unresolved

Even after (6), Q207's direct continuum bundle remains open. On the actual
common closed domain one still needs a displayed finite `A_R`
such that, schematically,

\[
 \max_{|I|\le R+1}\sup_{\mathcal D}
 w_I S^{-d_I}
 \left|D^I(F_0,U_0,\Pi_0,M_0,V_0,
                 E_1,U_1,W_1,T_1/\zeta)\right|
 \le A_R.                                           \tag{7}
\]

The polynomial factor `S=ell^2` must remain; Q207 already explains why it
cannot be discarded at the same pure-power exponent. Full pressure memory,
`Pi1/Pi2`, the whole moving pressure product, axial velocity and every
selected root/frame coefficient remain in the defining maps.

Equation (7) is only the seminorm part. The direct reconstruction must also
give one common positive lower bound `g_R` for all selected Jacobian, normal,
root, determinant, variance and support gaps, and every compact/slow inverse
constant in the archived pole family through Q207's common pole
`218,876,101,818` and its pressure/moment descendants.

Q192 and the downstream finite maps propagate jets once primitive arrays and
floors are supplied. The weighted-stress and regional records prove
orderwise finiteness relative to their retained exact-core and finite inputs.
They do not display (7), the common high-order floor bundle and every inverse
constant directly from the selected Volterra, ODE, pressure and moment maps.
Evidence hashes and executable integer recurrences cannot prove those
continuum suprema or minima.

Therefore the exact obstruction is narrower and more precise than Q207
stated: Q120 supplies the cutoff estimate and Q204--Q205 supply the
scale-choice quantifier, but one actual-field cutoff identity still has to be
bound. The selected-profile/root/frame/source/pressure/mean seminorms, common
gaps and operator archive also remain conditional.

No unconditional pulse hierarchy, complete residual, admissible force,
competitor exclusion or Navier--Stokes blow-up is proved. The next obligation
is to certify (6), then derive (7), the common gap floors and the archived
inverse bundle from the explicit selected maps before replaying Q207's finite
acyclic recurrence.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q208_cutoff_gevrey_frontier/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q208_cutoff_gevrey_frontier/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q208_cutoff_gevrey_frontier/evidence.json \
  --output fresh.json
```
