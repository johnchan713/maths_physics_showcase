# Complete Haar rows, the radial mean bound, and the actual stress identity

Parent: `eb39e7efe518ada0319fe8c26041f2a58d9ac286`.
The parent bounds the nonzero angular regular source and complementary
zero-Haar tangential rows at D=3/2-2kappa. This checkpoint supplies the
remaining Haar-row estimate **in terms of the actual selected stress
defect**, and derives a stronger fresh postslow scalar/radial Haar bound
from the existing coarse fields and sharper slow increment.

The important distinction is between the canonical principal stress identity,
already proved for a separately rebuilt canonical family, and an estimate
for any other inherited or localized field. Global scalar moments, normal
transversality, positive covariance and nearby frame/root values do not
establish that identity for a different field. We give its exact necessary
inequality and an independent null-moment tensor control. We do NOT assert
that the prescribed actual pulses fail it.

All historical evidence and actual fields, roots, supports, width C^-20 and
accepted band remain unchanged. No new acceptance gate is imposed. Original
coefficient/profile/support premises remain conditional; no complete residual
improvement, unconditional pulse family, infinite field or blow-up is certified.

## 1. The full current Haar equations, without initial-state shortcuts

Keep B=1-kappa, C=B+1/2, D=C-kappa, 0<kappa<1/4. The existing actual
mean bounds are tangential B and radial B+1, with b=epsilon*b_N
and V,G bounded. Let W be the CURRENT complete angular covariance,
including the original real curl square, first signed update and sourced
update. It has order one. Put bar for normalized auxiliary Haar average
and define the two physical selected flux defects

    Ftheta=bar W_rtheta-Sigma_theta,
    Fz=bar W_rz-Sigma_z.                                (1)

These contain the actual full stress, complete cutoff velocities and all
old/new crosses and squares. The fixed summed base and Pi1/Pi2 are identical
in the subtraction. No leading pressure is substituted for the full one.
The original exact velocity moments int R^2 bar(v)=int R bar(gamma)=0
are retained; they do not assert that these Haar velocities vanish pointwise.

With De0=partial_R+e/R and L0=partial_R^2+R^-1 partial_R
+epsilon^2 partial_Z^2, the exact generic Haar equations are

    gbar=-bar(t*beta)-D1^0[2b bar(beta)+bar(beta^2)+bar Wrr]
         -epsilon partial_Z[b bar(gamma)+G bar(beta)
                            +bar(beta gamma)+bar Wrz]
         +[2V bar(v)+bar(v^2)+bar Wtt]/R
         +epsilon(L0-R^-2)bar(beta),
    P=int gbar, pbar=T0(gbar-rho P),
    Ebar_theta=D2^0 Ftheta+Ntheta,
    Ntheta=-epsilon partial_T bar(v)
       +D2^0[b bar(v)+V bar(beta)+bar(beta v)]
       +epsilon partial_Z[G bar(v)+V bar(gamma)
                             +bar(gamma v)+bar Wtheta,z]
       -epsilon(L0-R^-2)bar(v),
    Ebar_z=D1^0 Fz+Nz,
    Nz=-epsilon partial_T bar(gamma)
       +D1^0[b bar(gamma)+G bar(beta)+bar(beta gamma)]
       +epsilon partial_Z[2G bar(gamma)+bar(gamma^2)+bar Wzz+pbar]
       -epsilon L0 bar(gamma).                          (2)

The sign convention for slow time is the parent's t*=-epsilon partial_T+c Nt;
only its absolute value enters the estimates. Haar kills the torus derivatives,
including graph derivatives of products. It does not turn bar(beta v) into
bar(beta)bar(v), set old Haar velocities to zero, or erase A1 correlations.
Haar pressure is the SAME full compensated pressure because T0 commutes with
Haar. Its A0 remainder has zero Haar average, not zero pointwise value.

## 2. Quantitative reduction to one spatial stress estimate

At output order n use current state/fixed/covariance norms through n+3.
Put d=n+4, r>=||R^-1||, vis0=2d^2+rd+r^2. The routine `haar_rows`
forms positive versions of (2), with

    S=2bx+x^2+W, H=bz+Gx+xz+W, T=2Vu+u^2+W,
    Gr=dx+(d+r)S+dH+rT+vis0*x,
    Pd=Kmom Gr, Pm=Kp(Gr+rho Pd),
    Rt=bu+Vx+xu, Rz=bz+Gx+xz,
    At=Gu+Vz+zu+W, Az=2Gz+z^2+W+Pm,
    Ntheta_poly=du+(d+2r)Rt+dAt+vis0*u,
    Nz_poly=dz+(d+r)Rz+dAz+vis0*z.                     (3)

Here symbols are nonnegative norm polynomials, not products of Haar means.
All mean correlations are bounded by full-field products AFTER multiplication.
The complete covariance W bounds every entry, not just the selected pair.
Gr and pbar have the coarse order B. Time, viscosity and the whole pressure
gradient gain one epsilon after Haar; the nongap radial mean fluxes have a
radial factor or b. Therefore

    Ntheta,Nz in M_(B+1), B+1>C>D.                      (4)

The possible mean-square terms have order 2B+1, and covariance axial terms
order two, both stronger than (4). There is no graph-kappa viscosity loss in
(3); it was necessary in the parent's UNAVERAGED estimate. The code's Haar
grammar is independent of the graph and fast-time multipliers.

Consequently a sufficient remaining spatial estimate is

    ||(Ftheta,Fz)||_n <= epsilon^D Q_n(S) zeta delta^-p.  (5)

For compact selected fluxes it is also necessary up to the existing finite
jet, pole and fixed-operator costs. Indeed (2) gives
De0 Fe=Ebar_e-Ne. The right side has exactly zero R^e moment because Fe
is compact. Apply the SAME two-edge primitive to this WHOLE difference:
Fe=Ie(Ebar_e-Ne). Its original weighted bound preserves the epsilon order.
One must not separately assume a nonzero-moment primitive has compact support.
This necessity argument uses no inversion of an edge weight and no new
first-band selection. It concerns all fixed jets, not a C^n contraction.

This identifies the missing spatial quantity rather than replacing it by
the three integrated scalar defects or a new full-row hypothesis.

## 3. Canonical completion and the same-field principal identity

Write W0 for the original real primary covariance, W0p for its principal
quadrature tensor, dWs for the COMPLETE first signed change, dWf for the
COMPLETE sourced change, and Sigma0 for the original leading stress.
Exactly, at the final mean step (mean updates do not change W),

    Fe=Le+Re,
    Le=bar(W0p)_re-Sigma0_e,
    Re=bar(W0-W0p+dWs+dWf)_re-(Sigma-Sigma0)_e.          (6)

The retained original-root and sourced bounds give respectively

| Remainder in (6) | Order |
|---|---:|
| REAL primary curl square | 2-2kappa |
| First signed complete covariance change | 2-3kappa |
| Sourced complete covariance change | C |
| Higher background stress | 3 |

Thus Re has order C for 0<kappa<1/4. This uses the exact real sine/cosine
primary identity for the first row ONLY; it is not used on the complex sourced
wave. The positive coefficient sum is implemented by `flux_remainder`.

If the canonical identity Le=0 holds for the SAME complete localized
primary and stress on the domain in question, (2)--(6) give Haar rows at C.
That is stronger than the parent's Q1 row D and nonzero-angular source D.
The exact canonical H*y=T construction and squared physical partition
identity establish this principal match for the separately rebuilt canonical
family in `canonical_initial_mean_join` and `core_primary_root_budget`,
conditional on their original valid core/column/root/cutoff premises.

There is NO need to identify two different homogeneous frames if the entire
subsequent construction uses that canonical family and its freshly recomputed
means. Conversely those theorems cannot be transplanted to an inherited
actual-normal primary merely by closeness, or from a local leading core to
a different localized tensor on the full support. All cutoffs, angular halves,
Haar normalizations, positive roots, squared partitions and full stress must
belong to the same fields. If an unmatched cutoff term is put into a tail,
provide its exact subtraction and tail bound; do not set it to zero in (6).

For a family without this identified exact match, the precise outstanding
inequality is

    ||epsilon^-1 Le||_n
       <= epsilon^(D-1) Q_n(S) zeta delta^-p,
    D-1=1/2-2kappa>0.                                  (7)

This is the actual principal tensor/stress discrepancy, not a request for
an independently assumed outgoing Haar row. The new two-rate executable
envelope retains Le at order one and Re at C. It cannot certify (7) by a
Boolean flag or subtract unrelated positive norm bounds. Original actual
field/support identification remains unverified here.

## 4. Exact null-moment stress control

On a positive compact annulus choose a nonzero smooth flat g, f=R^4 g^2,
er=(1,0,0), t=(0,R,-1), a=er+t, b=er-t and 0<eta<1. Let

    W0=epsilon*f[(3/2) a tensor a+(1/2) b tensor b],
    Weta=epsilon*f[(3/2+eta/2) a tensor a
                             +(1/2-eta/2) b tensor b]. (8)

Both tensors are positive semidefinite, real, smoothly compact and transverse
to the SAME n=(0,1/R,1). Their traces agree exactly because er.t=0.
Their difference has only

    D_rtheta=epsilon*eta*R*f,
    D_rz=-epsilon*eta*f,

and symmetric counterparts. It is independent of slow Z. Therefore its
full radial source and compensated pressure changes vanish exactly, as do
ALL THREE changes of P,Jtheta,Jz. Nevertheless

    dEbar_theta=D2^0(epsilon*eta*R*f),
    dEbar_z=D1^0(-epsilon*eta*f)                         (9)

are nonzero. Both integrated residual moments are zero by compact boundary
integration. The global moments miss a spatial stress error. The reference
selected stress is matched in this control; the perturbed one is not.
The moving pressure is kept in both states and is identical; its axial
gradient need not vanish even in the reference state.

Taking eta=1/L gives arbitrarily close positive covariance weights. With
L<=S=ell^2 and epsilon=exp(-h ell), an error of size epsilon/L cannot
in general be bounded by epsilon^D times a fixed polynomial in S: the
ratio contains exp(h(D-1)ell) divided by a polynomial. A mere O(1/L)
comparison therefore cannot prove (7). It may still be that the prescribed
family has an exact zero discrepancy; this control does not decide that.

These are covariance-map inputs, NOT asserted to be the selected pulse,
its invertible canonical cone/column matrix or its actual full compact curl.
Real amplitudes for (8) exist by its positive rank-one representation; compact
wave potentials add the separate real curl square at order 2-2kappa. That
higher-order correction cannot generically erase the order-one selected flux
defect, but no full-field lower bound for the prescribed construction is claimed.

## 5. The fresh postslow radial Haar obstruction does improve

The parent scalar join gives actual slow tangential increment C, radial C+1.
Old postfast total means still have coarse B/B+1 bounds. Apply the SAME exact
three-scalar subtraction after the fresh linear five-map cancellation.
Its tangential old/new correlations have order B+C, new squares 2C, and
time/base-radial terms at least C+2. Hence

    Pnew,Jtheta,new,Jz,new in S_(B+C), B+C=5/2-2kappa.   (10)

This improves the earlier coarse 2B reset by using the derived C increment;
no sharper old Haar-velocity estimate or incoming scalar callback is needed.
`postslow_scalars` computes the old complete postfast fields and this SAME
slow increment from existing recipes and uses the full nonlinear pressure
and axial moment formula. The original G=0/power-law V contract on the
ENTIRE slow potential hull is still required for its signed linear rows.

In the same base-subtracted mean system, for D0 T0 h=h-A0 h the exact radial residual is
-rho Pnew-A0(gr-rho Pnew).
The compensated source has zero radial Haar moment, so A0 has exactly zero
Haar average. Thus its radial HAAR residual is bounded at B+C by rho Pnew.
No pointwise A0 estimate or pointwise pressure B+C estimate is inferred.
The axial integrated residual keeps epsilon partial_Z(Jz+c_rho P), including
(partial_Z c_rho)P; it has one additional epsilon and one higher input jet.

Field/pole inputs are kept at the parent's 6p, products at 12p and the same
moment domains, with the explicit lower-pole embedding already derived.
The selected flux row joins use the proper two-edge pressure and primitive
operators; all constants retain the actual profiles, width and potential hull.
All finite jets required by the parent's A1 recipe still need original valid
inputs; neither this argument nor finite diagnostics independently encloses them.

## 6. Scope and next obligation

Established: the explicit full Haar-row two-rate envelope, the reduction of
Haar completion to (7), canonical conditional completion where its exact
same-field identity applies, and the actual postslow radial Haar bound (10)
relative to the retained mean-map/input contracts. No new analytical assumption
or accepted-band gate is marked proved. The null-moment tensor diagnoses a
failed implication; it does not refute the actual canonical construction.

Next resolve the same-field principal identity and localization bookkeeping
for the family actually used, or consistently carry out all joins in the
already rebuilt canonical family. Then join all exact tail bounds, including
pointwise A0/A1, Gaussian and spatial cutoffs, to the full physical residual.
Original actual input/support validity, forcing, infinite realization and
complete-velocity growth remain unproved.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output /fresh/path.json`.
The controls are exact polynomial/matrix, mixed-rate and provenance checks;
they are not continuum suprema or evaluations of selected enormous fields.
