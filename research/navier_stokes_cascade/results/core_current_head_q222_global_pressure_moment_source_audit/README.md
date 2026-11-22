# Q222 joined global pressure and moment source audit

Parent: `69153a8653e0c07794f5b4cce95ec943b9d1a649` (Q221).

Q221 completed the source audit of the compact, fast and slow inverse maps.
This checkpoint joins those maps to the pressure and moment identities they
consume.  It distinguishes two ledgers that use similar notation but have
different meanings:

1. the leading profile's five global defects `(M,R,J,S,Cp)`; and
2. the evolving mean scalars `(P,J_theta,J_z)` and their two integrated
   tangential residual rows.

The identities close on the frozen same-profile chain.  The audit corrects
the unsafe shorthand that identifies the finite endpoint cumulative row `I`
with the renormalized global heat moment `R`.  It also records that integrated
moments do not determine the full pointwise Haar rows.

## 1. Finite matching rows and the renormalized global row

The annular and compact matching maps act on the finite row

```text
(M,I,J,S,Cp),
```

with exact transforms

```text
T_a(M,I,J,S,Cp)=(M,J-4 eta I,I,-(S-8 eta M),Cp),
T_c(M,I,J,S,Cp)=(M,(M-J)/lambda,I,-S,Cp).             (1)
```

Their displayed inverses round-trip all five rows for `lambda>0`.  In the
terminal heat region the global angular row is instead

```text
R=int_0^infinity(H-H_pow)dX
 =lim_(X->infinity)[I(X)-X H_pow(X)/(1-h)].            (2)
```

The raw `I(infinity)` diverges.  If two repaired profiles are compared at
`X_rep`, then

```text
Delta R=Delta I(X_rep)
       +int_(X_rep)^infinity Delta(H-H_pow)dX.          (3)
```

Thus `Delta I(X_rep)=0` does not imply `Delta R=0` by itself.  The implication
used here is valid because the frozen repair provider proves that the two
profiles have the same future heat tail and the same power subtraction.  A
nonzero tail integral is a precise obstruction to dropping that hypothesis.

With the common tail retained, the same repaired profile satisfies the five
functional identities

```text
M_infinity=R=J_infinity=S_infinity=0,
Pi_axis+Cp_infinity=0.                                 (4)
```

In the unchanged gauge,

```text
Pi(X)=-int_X^infinity E(x)^2/(2x) dx.                  (5)
```

Intermediate heat, shape and compact edits can have nonzero pressure
memory; (5) follows from the completed five-row cancellation and is not a
pressure reset.  All five rows and their angular derivatives are required
for the exterior stress join.  A four-row repair can leave an independent
pressure derivative defect in the axial stress.

## 2. Compact mean pressure identity

For the actual current radial Haar source set

```text
P=int gbar dR,  Qg=int R^2 gbar dR,
pbar=I_0(gbar-rho P),  c_rho=(1/2)int R^2 rho dR.      (6)
```

The unit-mass compensator makes the primitive input have zero radial total.
Smooth compact extension removes its boundary terms.  Since
`partial_R pbar=gbar-rho P` on this auxiliary-independent mean input,
integration by parts gives

```text
int R pbar dR
 =-(1/2)int R^2 partial_R pbar dR
 =-Qg/2+c_rho P.                                      (7)
```

For general pointwise inputs the compact map retains its `A0` remainder;
(7) is not permission to delete that remainder from the pointwise radial
equation.  The original `Pi1/Pi2` pressure tails also remain in the complete
base.  They cancel only when subtracting the same complete base and gauge.

## 3. Both integrated tangential identities

Let `E_theta,E_z` be the current complete zero-angular residual rows after
recomputing every wave, fast and slow update.  Direct integration of their
conservative radial parts gives

```text
int R^2 Haar(E_theta)dR=epsilon partial_Z J_theta,      (8)
int R   Haar(E_z)dR
 =epsilon partial_Z(J_z+c_rho P).                      (9)
```

The radial identities are exact because

```text
int R^2 (partial_R+2/R)F dR=0,
int R   (partial_R+1/R)G dR=0                         (10)
```

for smooth compact extensions.  Here `J_theta` contains the complete
angular-axial flux.  The definition of `J_z` includes axial velocity
squares, `W_zz`, and `-Qg/2`; combining it with (7) yields (9).

The derivative in (9) acts on the whole moving product:

```text
epsilon[J_(z,Z)+c_(rho,Z)P+c_rho P_Z].                (11)
```

The middle term generally does not vanish.  Freezing `c_rho`, using only a
new pressure increment, or replacing the complete pressure by the leading
profile changes the identity.

## 4. Preservation through the joined maps

Nonzero-harmonic waves have zero angular mean and do not change the two Haar
velocity moments.  The Q221 fast temporal-inverse curl has zero auxiliary
Haar in all three components.  The fixed slow map imposes

```text
int R^2 u dR=0,  int R gamma dR=0,
D_1 beta+epsilon partial_Z gamma=0.                    (12)
```

Therefore (7)--(11) persist after sourced, signed, fast and slow updates.
All six covariance entries and the actual axial velocity are recomputed
after each update.  Products of zero-Haar fields can have nonzero Haar
correlations, so only the linear response averages vanish; their covariance
and pressure contributions remain.

For the two projected target rows, with `e=2` angular and `e=1` axial,

```text
M_2=epsilon partial_Z J_theta,
M_1=epsilon partial_Z(J_z+c_rho P),
int R^e sigma_e dR=1,
h_e=E_e-sigma_e M_e.                                  (13)
```

Thus `h_e` has exact zero weighted moment and the correctly oriented Q220
primitive preserves both flat edges.  The projection removes the scalar
moment from the primitive input; it does not remove it from the pressure or
mean system.  Exact zero moment is required at the opposite flat edge.

## 5. Scope and next obligation

This source audit joins the leading five-row global identities and the
evolving compact/fast/slow mean identities on the same frozen pressure gauge
and profile.  It retains the common-tail premise, the whole moving pressure
derivative, axial velocity, all six covariances, and `Pi1/Pi2`.  Exact tests
verify transforms, signs, product differentiation and frozen provenance;
they do not prove the continuum providers.

The result does not infer pointwise Haar-row estimates from moments, erase
the `A0/A1` tails, prove uniform generated-pulse bounds, improve the complete
residual unconditionally, or establish force exclusion or blow-up.

The next obligation is to audit the finite coefficient solves and
reserved-mean supports on this joined pressure/moment object, including the
nonremovable `Pi1/Pi2` pressure floor.
