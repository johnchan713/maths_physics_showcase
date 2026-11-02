# Q107: current sourced Haar flux and compensated pressure

Parent: `7a67d57f3bf7434323bef38b5c4a86a828844622` (Q106).
On 2 November 2026 its research CI 41504982812 and field CI 41504983092
were still running. Q105 field CI 41494912987 had passed. The inspected
remote head was Q106. This checkpoint binds the wave-induced part of
Q106's stress-subtracted Haar fluxes to the ACTUAL Q97/Q98 coefficient
factories. It does not assert that the older flux or the full final rows
have now been controlled. Historical evidence, root, accepted cover and
flat width C_setup^-20 are unchanged.

## 1. Subtract at the current fixed-mean stage

Use Q97's `post_first_sourced_wave_fixed_means` bundle. Both old and new
waves contain +/-1 and +/-2 with complete cut curls; delta W contains
BOTH complex cross orientations and the new square. The four incoming
mean histories, assigned stresses, base, rho and physical domain are
identical on both sides. Q106's averaged products therefore cancel
exactly in this subtraction. Put Dij=Haar(delta Wij),
De=partial_R+e/R and Dz=epsilon partial_Z. Then

    delta Htheta = D2 Drt + Dz Dzt,
    delta Hz = D1 Drz + Dz(Dzz+delta pbar).                (1)

These are the exact changes of the stress-subtracted fluxes, not an
estimate obtained by resetting the stress or discarding old products.
The full old Haar rows remain additive. Likewise both mean viscosities
and time terms cancel ONLY because the mean fields are fixed; they are
still required for Q105's later fast/slow differences.

For a smooth periodic coefficient f, with the original graph multiplier
independent of the auxiliary variables,

    Haar(Dr f)=partial_R Haar f.                          (2)

The graph term is a torus derivative with zero Haar integral. Thus (1)
has no graph multiplier epsilon^-kappa. This improves Q100's FULL row
increment grade D=3/2-2kappa only for its Haar part; the full nonzero-mode
rows and their existing bounds are unchanged. Haar contraction costs one,
not the factor two required for Qcirc.

## 2. Assemble pressure before using the zero moment

Q98's full covariance pressure identity, averaged on the SAME whole
characteristic hull [a,b], gives

    H=(Dtt-Drr)/R-Dz Drz,
    delta P=int_a^b H dr,
    delta pbar=-Drr+L0(H-rho delta P),
    L0 f=int_a^R f dr-chi(R,Z,T) int_a^b f dr.             (3)

The formula for delta P requires the original whole covariance trace
identity: the boundary of Drr vanishes. Seams or unmatched endpoints
require Q98's explicit defects and invalidate this simplified equation
until those defects are accounted for. All mixed trace identities remain
analytical requirements, not a consequence of a test passing.

On that contract, the SAME moving rho has int rho=1, so

    int (H-rho delta P) dr=0,
    delta pbar=-Drr+int_a^R(H-rho delta P) dr.             (4)

This is an exact simplification of the HAAR pressure only. Q106 retained
the whole moving L0 before imposing this identity; after imposing it the
chi term and its slow derivatives vanish together. They must not be
removed from an arbitrary unbalanced primitive. Nor does (4) remove
rho_Z delta P+rho partial_Z delta P from Dz(delta pbar). The local
-Drr memory survives and contributes -Dz Drr to the axial row. Old
pressure including the original base/Pi1/Pi2 subtraction remains in
the incoming row; nothing is reset to zero by this increment calculation.

In particular partial_R delta pbar-delta gbar=-rho delta P. There is
no surviving Haar A0 term once these exact contracts hold, consistent
with Q104. The nonzero-Haar A0 source can still be nonzero: (4) gives no
improvement for that compact modulation error. The implementation uses
Q98's already proved full weighted pressure envelope as a safe upper
bound for the averaged pressure; it does NOT substitute an unproved
weighted norm for the ordinary primitive in (4).

## 3. Bind coefficients to the current Q97/Q98 factory

Let F_n=Q98.through(bundle.change,n) bound every changed covariance entry
through n, with its raw binomial convolutions, original five-variable
normalization and covariance pole. Let Pm_(n+1) be the SAME Q98 shifted
pressure bound through n+1. On the common positive hull, let
r_n=sum_(j=0)^n a^(-j-1). Q106's derivative/Leibniz rule applied to (1)
gives sufficient positive coefficients

    delta Htheta_n <= (n+1) F_(n+1)+2 r_n F_n
                                        +epsilon(n+1) F_(n+1),
    delta Hz_n <= (n+1) F_(n+1)+r_n F_n
                          +epsilon(n+1)(F_(n+1)+Pm_(n+1)). (5)

All norms are on the SAME proved weighted chart. Haar averaging preserves
the weight because that original edge weight is auxiliary-independent.
No unweighted-to-weighted promotion or monotonicity of operator norms
is inferred. Coefficient arithmetic in `envelope` drops epsilon<=1 only
after retaining its favorable analytical power in the grade calculation.

`HaarChanges` requires Q100.SourcedRowChanges and its identical Q98
ShiftedPressure object. It calls pressure.bound(n+1,bundle) directly,
thereby retaining the original state, stage, pole archive, shape, time
inverse and mixed-trace gates. It requests covariance through n+2 and
uses bundle.change, not total W or an independently supplied fine Haar
callback. No new mean, target, root, pressure or support is introduced.
Output n needs the Q98 pressure archive through n+1, covariance jets
through n+2, and all the inherited Q97/Q96 response and native-profile
request DAGs. The larger final slow-map reserves in Q105 are not erased
by this smaller increment request. Native constants and actual provider
inequalities remain uncertified.

Under those original all-fixed-order input contracts, F and Pm have
order C=3/2-kappa. Radial derivatives in (1) keep C; both slow-gradient
terms have order C+1. Therefore

    delta Haar rows have sufficient order C=3/2-kappa,
    improvement over Q100 full-row bound = kappa.          (6)

This does NOT imply the previous generic conditional mean-row target
C0=2-2kappa. The remaining exponent gap is C0-C=1/2-kappa>0 throughout
0<kappa<1/4. This is a gap in available upper bounds, not a lower bound
on the actual row. To recover C0 one needs special cancellation or a
stronger bound for the ACTUAL combinations D2 Drt+Dz Dzt and
D1 Drz+Dz(Dzz+delta pbar), and the older Haar row must also be controlled.
A zero auxiliary mean of the sourced velocity does not force these
covariance means to vanish: opposite-sign equal harmonics survive in
both cross orientations, and the square survives as well.

## 4. Join and remaining obligation

Writing Hold for Q106's Haar row on the completed incoming first-real
state BEFORE the new sourced wave, Q105's exact final identity now reads

    final_theta=Hold_theta+delta Htheta+RFtheta+RStheta,
    final_z=Hold_z+delta Hz+A1 Qcirc Epre_z+RFz+RSz.       (7)

No stress term or old pressure appears twice: they are in Hold and cancel
in the fixed-mean increment. The A1 source is still Q105's current
pre-new-fast row; it is not replaced by Hold or by (1). Equation (7)
keeps all mean histories and full differences, not just their projections.
It is not full cone closure. The radial A0/scalar gates, original inner/
outer estimates, and physical scale transfers retain their obligations.

Tests compare subtraction of the full Q106 rows with (1), exercise (2)
with a nonconstant graph, check the unit-mass pressure cancellation and
surviving moving-rho/local-pressure derivatives, and verify cross/square
Haar contributions for both harmonics. Manufactured operator leaves run
the REAL Q98 bound and guards; they are not actual continuum data. The
positive envelope is checked against independent rational arithmetic.
Frozen replay binds all new source files to Q106 and its unchanged history.
Passing these checks does not prove the native inequalities.

Next bound Hold's ACTUAL assigned-stress cancellation from the original
current target/root and averaged full mean products, and determine whether
the sourced combinations in (1) have extra structure beyond (6). Q106's
positive-covariance comparison remains an information-gap control, not a
canonical-wave counterexample. Uniform actual-pulse/support estimates,
complete viscous/cutoff/interaction residual improvement, infinite
realization and blow-up remain unproved.
