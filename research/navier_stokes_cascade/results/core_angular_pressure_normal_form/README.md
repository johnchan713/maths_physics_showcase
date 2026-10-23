# Current angular-row normal form and the complete mean pressure

Parent: `29195306c4bccab9dc19dd5528ae03c2a6b84d51`.
Its mixed compact/time check was pending when inspected on 23 October
2026; the workflow was running without a failed step.

The parent bounds the centrifugal source through Gtheta=2(V+v)(-c Nt v).
We now use the already DERIVED final angular row instead of requesting
that enormous nonlinear source norm. The graph derivative of its radial
flux becomes an exact radial derivative under the linear centrifugal
multiplier. Its compact A0 image vanishes, while its pressure primitive
is retained explicitly. The remaining flux has order one rather than B.

This yields a conditional M_1 estimate for the COMPLETE CURRENT MEAN
pressure increment pm, improving its former coarse M_B estimate. It is
not an estimate for total pressure (which includes the base and nonzero
wave pressures), nor a new regular residual order. Actual weighted input,
Haar, high-jet, fixed profile and whole-support inequalities remain open.
The regular source stays E=2-3kappa and the physical exponent is still
negative. No absolute decay, unconditional pulse, infinite realization
or blow-up is proved. Historical frozen evidence and physical fields,
roots, width C_setup^-20, h, carriers and accepted band are unchanged.

## 1. Use the FINAL angular equation, which has no current pressure

Use one common current chart, with all coefficient derivatives holding
the label M,c,epsilon fixed. Set

    B=1-kappa, C0=2-2kappa, E=2-3kappa, 0<kappa<1/4,
    J=b v+V beta+beta v+W_rtheta,
    L=G v+V gamma+gamma v+W_theta,z.

W is the SAME complete six-entry current covariance, including every
old/new product. In the retained regular mean system the exact row is

    Etheta=tstar v+D2 J+Dz L
                   -epsilon(D_r D1 v+Dz^2 v)-D2 Sigma_theta,
    tstar=-epsilon partial_T+c Nt, Dz=epsilon partial_Z.       (1)

Both viscous parts remain: D_r D1 v=D_r^2 v+D_r v/R-v/R^2.
Sigma_theta is auxiliary-independent, as are b,V,G and the moving cutoff.
Put

    F=-Etheta-epsilon partial_T v+Dz L
                              -epsilon(D_r D1 v+Dz^2 v),
    H=c^-1 Nt^-1 Qcirc J, K=c^-1 Nt^-1 Qcirc F.               (2)

Nt commutes with D2 and Dz on this fixed chart. Equation (1) gives

    -c Nt v=D2 J+F-D2 Sigma_theta,
    Qcirc v=-D2 H-K.                                       (3)

No auxiliary-dependent multiplier has been commuted through Nt^-1.
In particular this is a linear identity for the FULL accumulated v;
the square v^2 will remain a separate term. No raw sum of differently
normalized old prefast rows is taken, so no cover/label transport is
silently assumed. Original physical maps/chart intertwiners still apply.
Only Qcirc F enters an inverse or compact estimate. Its whole-hull
support/extension is the required contract; no compactness assertion for
an unprojected auxiliary-independent stress part of F is needed.

The input Etheta in (2) is the parent's already derived FINAL angular
row at C0, from the exact fast cancellation and full slow difference:

    Etheta_final=Haar(Etheta_postwave)+Rf_theta+Rs_theta.       (4)

The A1 term is in the AXIAL component, not (4). The angular equation
and its full row differences contain no current pm and no current A0
tail. Thus (4) is available before the new radial estimate, so this is
not a circular inference from that radial estimate. The older incoming
source DAG still contains earlier compact/A1/pressure requests and must
not be truncated. No fine output callback replaces (4)'s lower recipes.
Additive Gaussian/base/exterior residual tails are kept in their assigned
ledger; they are not inserted into this regular-row temporal inverse.

## 2. A second exact radial quotient, with its pressure primitive

Keep the parent's FINAL S,Hrz,Rrest,K0 (Hrz is not the H in (2)):

    S=2b beta+beta^2+Wrr,
    Hrz=b gamma+G beta+beta gamma+Wrz,
    Rrest=-tstar beta-S/R-Dz Hrz+Wtt/R+epsilon Dz^2 beta,
    K0=-S+epsilon D1 beta,
    gr=a v+v^2/R+Rrest+D_r K0, a=2V/R.

The FULL radial multiplier is variable. Define

    a1=2a/R-partial_R a=6V/R^2-2V_R/R.

Since a is auxiliary-independent, D_r a=partial_R a, and the product rule
gives a D2 H=D_r(a H)+a1 H. Therefore

    K1=K0-a H,
    r1=a Haar(v)+v^2/R+Rrest-a1 H-a K,
    gr=r1+D_r K1.                                         (5)

Under the SAME whole-hull smooth compact extension, A0 D_r K1=0 and
T0 D_r K1=K1. The exact compact radial source is now

    A0 Qcirc gr=-A0(a1 H)-A0(a K)
                              +A0 Qcirc(v^2/R+Rrest).     (6)

Dropping V_R/R is incorrect even though V has no auxiliary dependence.
Dropping -a H from pressure is also incorrect. With the SAME moving
unit-mass compensator rho the complete current mean pressure is

    P=int Haar r1=int Haar gr,
    pm=K1+T0(r1-rho P).                                  (7)

H and K have zero Haar. Hence, EXACTLY pointwise in slow variables,

    Haar r1=a Haar(v)+Haar(v^2)/R+Haar Rrest.               (8)

All fast variance and old/fast correlations are retained. (8) equals
the previous Haar r0; the old derivative K0 has the same vanishing
global endpoint trace. P is not reset or weakened: its previously
derived finer C0+1 scalar bound remains available. The entire moving
rho P, local K0, complete axial pressure force and old memory stay in
(7). In particular Dz(-a H) generally is nonzero. The balanced base
still contains Pi1/Pi2, and nonzero wave pressure, assigned exterior
stress and additive tails remain as before. A0's annihilation does not
delete a pressure primitive from the full PDE.

## 3. Conditional order-one mean pressure from the retained inputs

Use the SAME original weighted current inputs, with fixed unweighted
background/profile multipliers retained separately:

    beta : M_(B+1), v,gamma : M_B, W_ab : M_1,
    b=epsilon b0, V,G bounded fixed multipliers,
    Haar(v) : M_1 from the preserved finite Haar-slow sum,
    Etheta_final : M_C0 from (4).

The sharper Haar bound is the EXISTING derived initial-slow plus later
slow sum, not a relabelled full M_B norm. Its actual same-hull weighted
identification is still a required inequality. Direct orders in (2) are

| Group | Order |
|---|---:|
| b v, V beta, beta v, W_rtheta in J | B+1, B+1, 2B+1, 1 |
| -Etheta | C0 |
| -epsilon partial_T v, Dz L | at least B+1 |
| -epsilon D_r D1 v | B+1-2kappa=E |
| -epsilon Dz^2 v | B+3 |

The retained graph floor M<=epsilon^-kappa justifies BOTH radial
viscous derivative losses. It follows that J in M_1, Qcirc F in M_E, H in
M_1 and K in M_E. Temporal inversion costs four raw derivatives and a
polynomial in S, not an epsilon power. Since E>1 and 2B>1,

    Rrest, K0, K1, r1, pm in M_1.                        (9)

Here the K0 viscous part has order B+2-kappa>1. Original T0 and moment
bounds prove (9) from (7). This is a NEW CONDITIONAL mean-pressure
bound using the extra angular identity, not evidence of actual all-scale
input validation. It does not improve total pressure, pointwise gr
(which still has D_r K1), the complete regular source E, or the parent's
direct A1 order B+1. The old coarse claims remain valid.

`pressure_envelope` supplies a positive factorial-through-order program
for (9) at output n, using current fields/backgrounds through n+6,
Etheta through n+4 and the original four-jet temporal constant with raw
factorial conversion, projection and weight costs included. For a
through-(n+6) bound g on j=d R^(d-1), the sufficient graph factor is
D=3(n+7)(1+g), retaining both j partial_U and M j vr.partial_y.
The a1 envelope also includes g in V_R=j V_U; replacing it by an
unscaled U derivative would lose a fixed-profile factor. Its actual
input coefficients are not evaluated or certified by this program.
The temporal argument is the original recipe, not an arbitrary inverse
norm: a sufficient raw/through factor is
48*6^n*2^(n+4)*(n+4)! times c^-1 and the retained weighted/profile
conversions. Here Gamma_n=24*2^(n+4), Qcirc costs two and 6^n bounds the
five-variable through-n conversion. `temporal_recipe` keeps these as
deferred integer instructions at huge n.

## 4. Explicit compact requests: still huge, on different sources

Apply the parent's SAME mixed compact/time rule to J with multiplier
a1 and to F with multiplier a. Apply the original scalar compact rule
to v^2/R and Rrest. For target t=C0, or the unchanged physical target
t=(N+3/2+2h+2m)/h, choose strictly positive integer counts

    pJ=floor((t-1)/kappa)+1, pF=floor((t-E)/kappa)+1,
    pQ=floor((t-2B)/kappa)+1, pR=pJ,
    JJ=m+2pJ+2, JF=m+2pF+2,
    JQ=m+2pQ+3, JR=m+2pR+3.                              (10)

For the displayed targets all counts are positive. At kappa=1/100000:

| Request, m=0 | J raw jet | F raw jet | v^2/R raw jet | Rrest raw jet |
|---|---:|---:|---:|---:|
| Normalized C0 | 200000 | 6 | 5 | 200001 |
| Physical h=1/200,N=0 | 60200004 | 60000010 | 60000009 | 60200005 |

The parent instead requested Gtheta through 200002/60200006; that is a
DIFFERENT function. No source high-jet gate has simply disappeared.
The new J kernel needs V through JJ+1 and the full original profiles:

    f0 a=2V/(d U),
    f0 a1=6V/(d U R)-2 V_U/R, U=R^d.                    (11)

All those reciprocal/profile derivatives remain explicit in the mixed
kernel norm N_J(f0 a1) or N_J(f0 a). F through JF needs current fields
through JF+2 and Etheta through JF. Rrest through JR needs beta through
JR+2. In these examples the largest LOCAL field index is still
200003/60200007, exactly the parent's local maximum. The entire older
source/root/normal/A0/A1 request graph is retained conservatively; no
total minimal original-profile derivative reduction is proved.

Retain the parent's FINAL 12q common domains, original weighted-to-raw
and output edge conversions, hull length, width/gap, cutoff derivatives,
carrier floors, source factorials and all archived shifts. All local
requests fit the older conservative graph. Any fixed multiplier pole
must be added before choosing q; operator monotonicity at different
poles is not assumed. No finite huge constant, factorial or actual
profile is materialized, and no derivative-dependent band is selected.

## 5. Remaining frontier and exact diagnostics

The required actual inequalities are now those for the current J,F,
v^2/R,Rrest at (10), including (4)'s same-field angular-row inputs,
the sharper preserved Haar norm, the V_R kernel and WHOLE compact
potential/flux traces. Low-jet manufactured fields do not prove them.
The parent's finite-jet obstruction remains valid for arbitrary mixed
row inputs; (5) uses additional angular-equation structure and does not
contradict it or establish that canonical source inequalities fail.

Next validate these actual reduced flux/profile/pressure inputs and
support identities before advancing the generic signed cycle. The full
regular order E, negative physical exponent, infinite diagonal,
localization/forcing regularity and velocity growth obligations remain.

Exact controls independently reconstruct the full angular and radial
rows, retain both viscosities, detect lost slow time/V_R/local pressure,
compare characteristic Fourier moments and pressure primitives, and
check independent epsilon norm budgets. Finite polynomial/Laurent
fixtures have finite endpoint regularity and are not canonical waves,
continuum norm enclosures or a PDE solution. Analytical derivations
(1)--(11), conditional on the stated retained contracts, carry the
mathematical content. Historical frozen records are bound unchanged.

Reproduce: `python test_audit.py -v`, then `python audit.py
--verify-record evidence.json --output FRESH_PATH.json`.
