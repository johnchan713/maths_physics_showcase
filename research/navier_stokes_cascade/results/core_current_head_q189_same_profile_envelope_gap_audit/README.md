# Q189: same-profile compact envelope and regional gap reconstruction

Parent: Q188, `7457c5108177a9767b422fc9cc92efa8dcde5812`.
The live branch, Q188, its named analytical providers, frozen evidence and
successful current-head CI run were read first. Every selected parameter,
root, support, pressure gauge, frequency and the `C^-20` activation width is
unchanged.

Q188 reduced compact modulation and exact repair to strict hierarchy
inequalities and isolated two remaining displayed inputs on one selected
profile: the compact derivative/denominator envelope at most `A=C^4096` and
the applicable regional cone gaps at least `A^-1=C^-4096`. This checkpoint
independently reconstructs the exponent arithmetic and gap comparisons from
the named analytical regional records. Both inequalities close. The closure
is relative to those continuum theorems; this checkpoint does not rederive
every primitive regional theorem from the initial ODEs and roots.

## 1. Independent compact exponent ledger

Use the factorial logarithmic-radial jet algebra of the actual compact-input
envelope. Starting with field, moment and pressure exponents

\[
 64,\qquad 256,\qquad 257,
\]

and retaining differentiation, sums, products, reciprocal floors and the
full pressure coordinates gives the following largest entries:

| quantity | exponent |
|---|---:|
| `E^-1` and `H^-1` in `J1` | 73, 75 |
| normalized moments in `C2` | 546 |
| pressure coordinates `p1,p2` in `J1` | 351, 414 |
| shear coordinates `a,bs` in `J1` | 140, 139 |
| `a^-1` in `J1` | 301 |
| loop coordinates `ts,vs` in `J1` | 557, 1371 |
| `Pc,Jc` in `J1` | 1070 |

The maximum is therefore

\[
 1371<4096.                                             \tag{1}
\]

The executable recomputes these costs without importing the earlier envelope
implementation. In particular, reciprocal differentiation retains the
positive floors, and `p1,p2,Pc,Jc` retain the cumulative pressure row. Thus
Q188's displayed bound (11) holds for the actual compact input theorem.

## 2. Exact correction of the outer-ramp comparison

The regional theorem corrected an older kernel sign mistake and needs only
`p1>.5 XR`. The present audit replaces its decimal check by exact rational
bounds. The alternating series and the elementary series estimate `e<3`
give

\[
 e^{-3/5}>1-\frac35+\frac{(3/5)^2}{2}-\frac{(3/5)^3}{6}
 =\frac{68}{125},\qquad e-1<2.
\]

Consequently

\[
 \frac{p_1}{X_R}>
 \frac{11}{10}\frac{68}{125}-\frac1{98}\,2
 =\frac{17701}{30625}>\frac12.                          \tag{2}
\]

This proves the weaker floor actually consumed downstream. It does not revive
the superseded stronger expression that multiplied the negative source term
by the positive-data kernel lower bound.

## 3. Regional gaps on their correct subsets

On the activation collar the analytic provider gives

\[
 a>C^{-22},\quad P_c>2.2,
 \quad P_c-v>C^{-401},
 \quad 2(P_c-v)^2-(v-2)J_c^2>C^{-802}
\]

where the last two gates apply, namely where `v>=2`. The joined common bound
is `C^-1024`. Since

\[
 22,401,802,1024<4096,                                  \tag{3}
\]

each applicable lower bound exceeds `A^-1=C^-4096`. On the continuation,
shape and joining rows `v<1`, so imposing the last two gates there would be a
category error rather than a stronger audit.

The later regions have order-one slack. For the axial stage the retained
finite-radius factors give

\[
 2(35)^2\frac{47}{50}\frac{417}{1000}>1,
\]

and on the intermediate stage

\[
 2(12)^2\frac{627}{2000}>1.                              \tag{4}
\]

The protected right collar has `v-2=2 lambda>C^-1`, also larger than
`C^-4096`. Equations (2)--(4), together with the activation floors, close
Q188's displayed bound (12) on precisely the subsets where it is required.

## 4. One full-pressure profile across the joins

This is not a splice of unrelated estimates. The core contraction uses the
same exact pressure family and zero leading stress; the growth checkpoint
keeps the same core and `C^-20` activation width; continuation retains the
complete forward pressure terms; and the shape/join checkpoint uses the same
exact `C2` root and backward pressure ledger. At the compact repair endpoint

\[
 X_{rep}=\rho_p e^5=X_*                                  \tag{5}
\]

and all five cumulative rows `(M,I,J,S,Cp)` are restored. The global moment
ledger then retains `Pi=Pi_axis+Cp`, including nonzero intermediate memory,
and propagates the same gauge into the compensated exterior. Field equality
alone would not establish this conclusion.

The resulting statement is therefore: relative to the named analytical
regional theorems, the selected profile simultaneously satisfies Q188's
`C^4096` compact envelope and `C^-4096` regional gap requirements, and the
already-audited single finite frequency and exact repair apply to it.

## 5. Correction, limits and next obligation

Q188 correctly said that it had not reconstructed (11)--(12), but its phrase
"surviving independent-validation boundary" is now discharged at that
level. This checkpoint independently redoes their exponent and gap arithmetic
and binds their outputs to the same full-pressure profile. It still relies on
the named continuum provider theorems; frozen flags and passing tests do not
constitute a new first-principles proof of every regional ODE estimate.

No prescribed historical first band, unconditional uniform physical pulse
bound, complete residual improvement, admissible force, competitor exclusion
or Navier--Stokes blow-up is proved here. The next justified obligation is to
bind this same-profile envelope/gap object through Q177--Q185 and identify the
first downstream premise not already supplied by an analytical provider.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q189_same_profile_envelope_gap_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q189_same_profile_envelope_gap_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q189_same_profile_envelope_gap_audit/evidence.json \
  --output /tmp/q189_same_profile_envelope_gap.json
```
