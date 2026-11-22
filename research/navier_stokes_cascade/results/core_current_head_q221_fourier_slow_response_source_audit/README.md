# Q221 zero-Haar Fourier and fixed slow-response source audit

Parent: `376836e797c2dac1899fef7c83a89786dc953a5e` (Q220).

Q220 source-audited the shifted two-edge compact kernel.  This checkpoint
audits the two remaining inverse maps consumed by Q216: the zero-Haar
Fourier inverse and the fixed five-profile slow response.  Both preserve
their declared weighted domains at every fixed order and every fixed finite
pole.  The result retains the actual axial repair and complete pressure
inputs.  It also records the exact boundaries of the theorem: zero Haar is
necessary, and the slow axial constant is not uniform as the selected
positive `lambda` tends to zero.

## 1. Exact quadratic divisor and the four-derivative loss

On the auxiliary two-torus use the fixed directions

```text
v_r=(1,-b),  v_t=(b,1),  b=sqrt(2)-1.
```

For a nonzero integer mode `k=(k1,k2)`, multiplication by the quadratic
conjugate gives the nonzero integers

```text
(v_r.k)(v_r^*.k)=(k1+k2)^2-2 k2^2,
(v_t.k)(v_t^*.k)=(k2-k1)^2-2 k1^2.
```

Neither integer can vanish away from `k=0`: an integer square equal to twice
an integer square would make `sqrt(2)` rational (equivalently, infinite
descent removes a common factor two).  Since the conjugate factor is bounded
by `(1+sqrt(2))|k|_1<3|k|_1`,

```text
|v_r.k|^-1, |v_t.k|^-1 < 3 |k|_1.                    (1)
```

For an output torus derivative of order at most `m`, integrate the Fourier
coefficient by parts `m+4` times in a coordinate with
`|k_j|>=|k|_1/2`.  The output derivative and the divisor in (1) leave a
per-mode factor `|k|_1^-3`.  The `l1` shell of radius `r` has `4r` points,
and

```text
sum_(r>=1) r^-2 <= 2                                  (2)
```

by comparison with `1/[r(r-1)]` for `r>=2`.  Therefore

```text
N_m(c^-1 N_t^-1 f)
  <= c^-1 Gamma_m N_(m+4)(f),
Gamma_m=24*2^(m+4).                                   (3)
```

This proves the stated fixed-order loss without evaluating the enormous
constant at the largest archived order.  Slow and radial derivatives commute
with the Fourier multiplier because the directions are fixed and the
multiplier acts only on the torus.  Jets of the varying reciprocal `c^-1`
remain separate finite coefficient factors; they are not silently assigned
to the source.

The zero-Haar hypothesis is essential.  At `k=0` the multiplier
`(v.k)^-1` is undefined.  Moment cancellation in radius or compact support
does not remove this obstruction.

## 2. Fast field and complete axial row

For the current zero-Haar tangential and axial rows, the exact update remains

```text
u=-c^-1 N_t^-1 F_theta,
g_d=-c^-1 N_t^-1 F_z,
psi=T_1 g_d,
beta=-epsilon partial_Z psi,
gamma=D_1 psi=g_d-A_1 g_d.                            (4)
```

Thus the axial velocity includes the generally nonzero `-A_1 g_d` term.
Combining (3) with the Q220 compact estimate requests source order `m+4`
for `u,g_d`, `m+5` for the radial potential, `m+8` after the fresh scalar
and slow steps, and `m+9` for the whole axial pressure moment.  These are
fixed finite requests, not order-uniform constants.

## 3. Fixed five-profile slow response

Let the selected smooth base bump generate the three angular profiles
`F_P,F_z`, the axial profile `F_g`, and its compact potential `Psi`.  With
`A=1/2+h`, the distinct physical powers are

```text
u=a^-1[s^A P F_P+s^(A-1) J_z F_z],
gamma=a^-1 s^(A-3/2) J_theta F_g,
psi=a^-1 s^(A-1) J_theta Psi,
beta=-epsilon partial_Z psi.                          (5)
```

The fixed matrix identities give the five moment rows and
`D_1 beta+epsilon partial_Z gamma=0`.  The powers in (5) cannot be replaced
by one common scale.  The profile potential may fill gaps between its bumps,
so the support is the whole convex potential hull.  Q30 places that hull at
a fixed positive distance from both shell edges.  On it, the flat weight has
a positive lower bound.  Finite profile and coordinate chain rules therefore
map every fixed finite input pole to the same output pole.

The coefficient `a^-1` uses the retained positive background floor and its
fixed-order reciprocal recurrence.  The operator theorem consumes those
coefficient jets; it does not prove the background floor or the scalar
inputs.

## 4. Precise conditioning obstruction

For `g=exp(d)>1`, `v=exp(-2d lambda)` and the selected `lambda>0`, the axial
two-profile determinant contains

```text
g(1-v)=g(1-exp(-2d lambda))>0.                         (6)
```

Hence the fixed axial response is finite.  Its singular factor satisfies

```text
(1-exp(-2d lambda))^-1
  <= exp(2d lambda)/(2d lambda),                       (7)
```

because `exp(x)-1>=x`.  Conversely it is asymptotic to
`(2d lambda)^-1`.  Therefore no constant uniform as `lambda` decreases to
zero follows.  This does not obstruct the construction, which selected one
positive `lambda` before requesting the inverse constant.  The three-profile
angular matrix remains nonsingular at `lambda=0`; the obstruction is in the
axial block.

## 5. Pressure memory and scope

The current scalars entering (5) are recomputed from the complete state.
The axial row retains

```text
epsilon partial_Z(J_z+c_rho P)
 =epsilon[J_(z,Z)+c_(rho,Z)P+c_rho P_Z],              (8)
```

including `c_(rho,Z)P`, the actual axial velocity, all six covariance
entries, old pressure, and both `Pi1/Pi2` tails.  The Fourier and slow
operators preserve these inputs and their fixed support; they do not create
the weighted source estimates that feed them.

Equations (1)--(8) support Q216's remaining inverse quantifiers on the typed
domain: all fixed derivative orders and all fixed finite poles.  Passing
tests check exact norm identities, request arithmetic, source typing and the
frozen dependency closure; they are not a continuum proof.  This checkpoint
does not establish nonzero-Haar inversion, `lambda`-uniform response,
unconditional generated-source bounds, a complete residual improvement,
force exclusion or blow-up.

The next obligation is to audit the global pressure and moment identities on
the joined compact/fast/slow object, retaining both exact tangential moments,
the whole moving-pressure derivative, axial velocity, all six covariances,
and `Pi1/Pi2`.
