# Q28: derive the initial mean input by a second full mean-map composition

Parent: `011261b106abba9fbb23a0443218ff555d2f751c` (Q27).
Q27 removes four independent wave inputs relative to the original canonical
root/target/phase providers. Here the independent INITIAL MEAN callback of
the first signed source is removed as well. The same complete mean recipe
is applied first to the primary with zero old means, then to the signed
wave with those derived means retained. This is a conditional analytical
size estimate, not an improvement of the full residual or a certification
of all underlying fixed-map norms. No field or historical record changes.

## 1. Identify the actual two passes, without changing the target

The initial scalar seed provider (2)--(4) starts from zero mean velocities,
forms all six entries of the actual complete primary covariance, applies
the zero-Haar fast map, recomputes pressure and the three scalars, and then
applies the fixed slow five-map. This is exactly the operation order of
`core_signed_prefix_source_budget.PrefixBudget.final` with initial=0.
Its general complete mean equations reduce to the seed provider's equations
when beta=v=gamma=0. There is no old mean time or viscosity term to drop in
that first evaluation; later evaluations retain them in full.

Use the Q27 wave construction with xi identically zero for this primary
pass. In particular use the ACTUAL cut complete amplitude U in covariance,
not just its principal part. The original canonical roots and full pressure
remain. All common-atlas/label conversion factors must already be included,
as specified by the original aggregate provider; zero initial means do not
supply a missing label-overlap bound.

After this initial pass compute the original signed target from the sharper
initial Haar/scalar providers, retaining the moving c_rho P moment. Q28
uses their existing target estimate; it does not replace it by the coarse
size estimate derived here. Form the signed wave once using that PRE-signed
target and the original roots. The SECOND mean pass takes the first pass's
complete means as its old state. Its new fast/slow increments are ADDED.
There is no fixed-point equation and no post-signed target feedback.

## 2. A coarse bound suffices for the initial means

Put mu=1-kappa, 0<kappa<1/4. The complete primary is bounded at order
1/2 with weight sqrt(zeta)P, so each covariance entry is bounded at order
1 with full weight zeta and a finite pole. With zero means, the full radial
row and tangential rows are bounded at order mu by (3) of the signed-prefix
provider. This retains the whole reconstructed pressure

    Pdef=int Haar(gr), pm=T0(gr-rho Pdef),

before its axial derivative is used. Both background stress components
remain in the full rows. The actual fast inverse uses their zero-Haar parts.

There is no missing projection constant in applying the retained Fourier
inverse estimate to a full row E: for every nonzero torus frequency,
(E-Haar E)_k=E_k. The proof of Gamma_m=24*2^(m+4) in the fast-inverse
provider integrates this coefficient by parts at least four times in a
torus variable. The torus derivative of Haar E is zero, also after slow
parameter differentiation. Thus exactly the same Fourier sum is bounded
by raw full-row jets through m+4. This observation applies to Nt^-1(I-Haar),
not to an arbitrary operator on a projected function. In particular it is
not a claim that the sup norm of I-Haar is one.

The weighted shifted compact operators and temporal inverse consequently
give tangential fast means at order mu and radial fast mean at mu+1.
Recompute gr, pm, Pdef, Jtheta and Jz from the updated actual state and the
unchanged complete primary. The same full-row polynomial bounds their
scalars at order mu. The fixed slow response then has tangential order mu
and radial order mu+1. Adding both contributions proves

    |beta_initial|_n <= epsilon^(mu+1) zeta delta^-2p X_n(S),
    |v_initial|_n    <= epsilon^mu zeta delta^-2p U_n(S),
    |gamma_initial|_n<= epsilon^mu zeta delta^-2p Z_n(S). (1)

Here (X,U,Z)=first_pass.final(n) is an explicit finite positive-polynomial
recipe from the original primary and fixed operators. This conservative
bound is weaker than the seed's sharper slow rates and initial Haar defect
rates. Those sharper results remain necessary for the signed target.
No later Haar cancellation is inferred from (1).

The initial weighted-map proof, rather than unweighted norm smallness,
supplies the operator premises. The real axial velocity is gd-A1 gd.
Slow profiles restore the full weight only on their entire reserved
potential hull, not just on individual bump supports. The constants include
its weight floor, actual cutoff gaps, fractional physical scales and lambda.
They need not be small, or uniform as these fixed parameters vary.

## 3. Two pole archives and the full request depth

Choose primary input pole p covering all requested primary jets. The first
pass uses covariance/fast-row pole 2p, postfast means 2p, fresh pressure/row
pole 4p, and slow response on a fixed interior hull. It yields (1) at pole
2p. If p_signed covers all Q27 signed-wave requests, the second pass must
start with

    r=max(2p,p_signed),
    final mean pole=2r, source pole=3r.                 (2)

Retain fast constants at the distinct input poles {2p,2r}; retain pressure
and moment constants at {2p,4p,2r,4r}, with the correct derivative requests.
A constant for the largest pole is not assumed to dominate all the others.
Separate provider dictionaries in `compose` preserve these two archives;
if common constants are used, pad by their sum at each requested order.
The helper's pole ledger does not manufacture the actual constants.

For a mean output n, the first recipe's slow radial branch calls fresh
rows at n+1, which read postfast fields at n+4, whose radial fast branch
reads rows at n+9 and primary wave inputs at n+12. The signed source at
order j requests its initial means through j+15. Therefore this route uses

    primary complete-wave inputs through j+27,
    primary root/amplitude inputs through j+28,
    signed complete-wave inputs through j+15,
    signed root/target inputs through j+16.             (3)

The additional one is the complete curl. Counts (3) are BEFORE expanding
the signed target, homogeneous ODE, phase and background dependencies.
They do not replace Q14's much larger fully expanded primitive ceiling.
They do show why checking primary jets merely through j+16 is insufficient
for this particular two-pass route.

## 4. What the composition establishes, and the remaining gate

With (1)--(3), the unmodified signed source polynomial now takes neither
an independent initial mean bound nor independent T,R,Pi,U wave bounds.
It yields the conditional estimate

    |D^I f_m| <= epsilon^mu sqrt(zeta) delta^(-3r) P
                        * j! F_j(S), |I|=j.            (4)

All old/new wave products, fresh pressure, Pi1/Pi2 base subtraction,
A0 additive error, axial A1 velocity, viscosity and cutoff-curl terms keep
their original routing. Q26 supplies the fixed-q domain statement; Q20--21
supply the subsequent Gaussian quotient and lifted-path transfers.

The remaining operator gate is substantive: at every distinct archived
pole and requested order, the ACTUAL shifted T_e/A_e, pressure and moment
maps and slow five-map must satisfy their weighted jet inequalities on
the full radial hull with band-independent polynomial constants. The
initial-Haar-row provider (3)--(4) gives a recipe in fixed profile/gap/chart
seminorms; those hypotheses must be bound to the same exact-h construction.
Also retained are the actual weighted full stress and finalized background
schedule, canonical value gaps, common atlas and sharp initial-target
premises. Neither two-pass composition nor finite test fixtures proves them.

Next bind these fixed operator/profile constants and support hypotheses to
the actual full-hull construction at both archives, or document the exact
failed inequality. Complete all-scale source realization, historical
localization, full residual improvement and blow-up remain unproved.

Tests check both stage traces, nonzero retained old means, fresh pressure,
raw-source composition, exact Fourier projection and separate pole archives.
They use manufactured coefficients, not selected-field evaluations. Run
`python test_audit.py -v` and `python audit.py --verify-record evidence.json
--output FRESH.json`. Frozen Q27 evidence and its dependencies are replayed.
