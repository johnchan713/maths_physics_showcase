# Frozen-phase coefficient inputs and the precise normal-floor gate

Parent: `12d02e2b988c538a66b0f9d93bcd5e0444c6d743`.
The preceding source and inverse budgets require normal, moving-normal,
matrix, logarithmic-profile and eikonal coefficient jets. This checkpoint
replaces those independent jet inputs by positive polynomial recipes from
the FULL background's slow jets and the prescribed frozen phase. It also
derives the normal floor from a concrete whole-path value inequality.

The lower full-background seminorms, leading uniform growth/shear gaps,
phase-comparison constant, tangent propagator and support premises remain
conditional. No actual giant parameter or selected field norm is sampled.
No accepted band, original profile, tilt, carrier, width C^-20, radius,
cutoff, support or historical evidence is changed. U below is merely a
rational UPPER BOUND on the already fixed tilt u, not a replacement tilt.

## 1. Freeze the actual carrier before differentiating

Use the existing positive shell Rmin<=R<=Rmax. The label freezes p,pz,x0,
lambda0,L and sign sigma, with 0<lambda_min<=lambda0<=lambda_max and
L>=Lmin>0. The path coordinate v lies in [0,L]. In common coefficient
coordinates (R,Z,T,H1,H2), v is affine in the auxiliary coordinates.
Let V_j(S) bound its RAW mixed derivatives, including V_0>=L and the
actual common-lift first-derivative cost. Higher affine jets vanish.
Their constants retain the existing covering and c_i^-1 costs. All fields
F,G,b=epsilon*btilde below are the FULL summed background in these
coordinates; fixed slow jet bounds are the remaining lower inputs.

Let gmin be a positive lower bound for the frozen leading |g0|. Set
Bplus=max(1,lambda_max), which bounds sqrt(lambda_max), and

    qz=Bplus[1+U/(Lmin*gmin)],
    qtheta=Rmax*qz+1, qrad=Bplus*U/2.                    (1)

These bound |pz|,|p|,|x0| in the actual phase construction. Indeed Bs^2=
lambda0/[epsilon k^2(1+u^2)^(3/2)]<=lambda_max, epsilon k^2>=1, and
the carrier's leading shear correction has length at most U/(Lmin*gmin).
The final +1 retains |p-ptilde|<=1/k<=1 for nearest NONZERO rounding.
The invalid universal half-unit rounding error is not used. There is no
positive lower bound on |p| or on a tangential COMPONENT of the normal.
Frozen carriers have zero continuous derivatives, including at tie labels.

## 2. Actual normal, moving normal, K and Eik from full F,G,b

The actual phase is

    Phi=p theta+pz Z/epsilon+x0 R-v H, H=pF+pzG,
    n=(x0-v H_R, p/R, pz-epsilon v H_Z),
    n'=(-H_R,0,-epsilon H_Z),                            (2)
    K=[[0,-2F,0],[2F+R F_R,0,0],[G_R,0,0]].

The graph/lift identities Dr v=Dz v=0, t* v=1 are the inherited actual
geometry premises. They give the exact COMPLETE identity

    Eik/epsilon=v(H_T-G H_Z)+btilde*n_r.                 (3)

The b term is not zero. Both G transport and the axial normal term stay.
No leading background is substituted before taking derivatives. Formula
(3) keeps the epsilon factor before any S estimate.

For raw mixed jet maxima define convolution

    (A*B)_j=sum_(i=0)^j binom(j,i) A_i B_(j-i).

Let F_j,G_j,B_j bound the full slow jets of F,G,btilde, with no assumed
small higher derivatives. Define R_0=Rmax, R_1=1, R_j=0 for j>=2 and
r_j=j! Rmin^(-j-1). Sufficient explicit recipes are

    h_j=qtheta F_(j+1)+qz G_(j+1),
    nr_j=1_(j=0)qrad+(V*h)_j,
    nt_j=qtheta r_j,
    nz_j=1_(j=0)qz+(V*h)_j,
    N_j=nr_j+nt_j+nz_j, Nprime_j=2h_j,
    K_j=4F_j+(R*Fnext)_j+G_(j+1),                       (4)
    E_j=(V*[h+G*h])_j+(B*nr)_j,
    bR_j=B_(j+1), Bz_j=B_(j+1)+(R*Fnext)_j+G_(j+1).

Here E_j bounds epsilon^-1 Eik, bR bounds epsilon^-1 b_R,
and Bz bounds epsilon^-1 Dz(b,RF,G); epsilon<=1 has only been dropped
in stronger terms. Euclidean vector/operator norms are bounded by these
component/entry sums. Repeated raw mixed derivatives obey the displayed
binomial convolution by grouping multibinomial coefficients by order.
No extra spatial graph derivative is charged to an ordinary coefficient
or to n'. These upper bounds need not capture n'=O(1/S) sharply enough
for the separate moving-frame perturbation proof; that proof remains a
lower dependency, not a consequence of these coarse jet norms.

Each order j of (4) uses the background only through j+1. All S losses
come from supplied finite jet polynomials and the actual common-lift
V_j. There is no new negative epsilon power. To feed the parent's
FACTORIAL through-order source norm from these RAW maxima use

    C_kind,n=sum_(j=0)^n (5^j/j!) kind_j.                (5)

There are five coefficient variables; sum_(|alpha|=j)1/alpha!=5^j/j!.
Angular basis-connection derivatives belong to the already retained
cylindrical operator, rather than an extra independent coefficient variable.
`source_coefficients` supplies b,G,V=RF,r,N,K,Eik,bR,Bz using (4)--(5).
The remaining graph/time/cutoff, stress and weighted-map constants remain
their existing inputs. This is a substitution into the earlier source
recipe, not a different source or a hidden new current residual norm.

## 3. The logarithmic pulse rate has varying growth at EVERY order

The inherited reference profile satisfies

    P'/P=lambda(v)-dref(v),
    lambda(v)=lambda0/sqrt(1+s_sigma(v)^2),
    dref=lambda0(1+s_sigma(v)^2)/(1+u^2)^(3/2),
    s_sigma(v)=sigma*u*(1/2+v/L).                        (6)

The growth in (6) is NOT the frozen value lambda0. In particular the
damping is quadratic in v but the growth generally has nonzero third
and higher derivatives. Those derivatives cannot be set to zero in the
Gaussian-adapted extension. This is explicit bookkeeping of the actual
inherited formula, not a change to its physical profile.

Let s_0=3U/2, s_j=U V_j/Lmin for j>=1, and A_j=1_(j=0)+(s*s)_j.
Since the ACTUAL A=1+s_sigma^2 is at least one, finite Taylor composition
at each real point gives raw inverse-square-root bounds

    J_0=1,
    J_n=n! sum_(k=1)^n [binom(2k,k)/4^k]
                  [z^n](sum_(j=1)^n A_j z^j/j!)^k.      (7)

Only finitely many positive powers occur. This is a finite jet formula;
no analytic continuation or convergence of the full background Taylor
series is assumed. A sufficient logarithmic-rate bound is

    p_j=lambda_max(J_j+A_j),                            (8)

where (1+u^2)^(-3/2)<=1 was safely dropped only in the damping bound.
These p_j replace the independent `profile_rate` inputs in the prior
ambient extension. Actual lambda0,L and u are still frozen per label;
their continuous derivatives have not been invented or discarded.

## 4. The floor is a whole-path value gate, not nonzero angular frequency

Use the retained actual phase comparison

    |n-Bs(s_sigma,Kunit)|
      <= Cphi[S^-1+S(S^-3+epsilon+epsilon^2+k^-1)].       (9)

If S>=1 and S^2(epsilon+epsilon^2+k^-1)<=1, (9) is at most 3Cphi/S.
The rounded error acts on the WHOLE path; its S factor is not dropped.
Since epsilon k^2<=4 and u<=U, the rational outward bound

    beta=lambda_min/[4(1+U^2)^2] <= Bs^2

is valid. The squared denominator exceeds the power 3/2 denominator.
The precise remaining comparison is

    S^2(epsilon+epsilon^2+k^-1)<=1,
    36 Cphi^2 <= beta S^2.                              (10)

Then 3Cphi/S<=sqrt(beta)/2. The reference tangential vector has length
Bs, hence triangle inequality gives |n_tan|>=sqrt(beta)/2 and

    |n|^2>=beta/4, |n|^-2<=4/beta.                       (11)

This supplies the actual inverse's full-normal floor from the value
comparison. An arbitrarily small but nonzero angular carrier alone does
not give (11). `normal_gate` records both inequalities and the rational
floor; the forced-inverse constructor rejects a failed or inconsistent
gate. Passing manufactured gate fixtures does not verify actual Cphi,
lambda_min or their substitution on the selected domain.

The SAME accepted-domain selection must already imply (10) for actual
constants, as required by its original normal/frame value gates. No
new band, derivative-dependent threshold or profile is selected here.
Until that substitution is independently verified, (10) is the precise
remaining inequality, and (11) is conditional. The core leading gaps,
full background bounds and tangent propagator are not closed by (11).

## 5. Join and remaining obligations

`PhaseBudget.forced` supplies raw N,Nprime,K,p jets and the floor (11)
to the existing Gaussian-adapted inverse. It retains moving-normal and
UNPROJECTED source pressure contributions. The ambient lift at value
order needs normal jets through one, hence background jets through two
even for output order zero; higher fixed orders use the finite larger
requested set. `source_coefficients` supplies (5) to the signed-prefix
source norm, with exact positive-polynomial coefficient transfer between
independently loaded historical helpers. No rounding to floats occurs.

The actual F,G,btilde input constants and Cphi remain to be derived and
substituted from the full core/modulation/repair estimates on the unchanged
domain. The homogeneous tangent propagator, original-root wave amplitudes,
closed-label/full-path support and remaining fixed operators also retain
their stated lower premises. After those inputs, the first sourced
complete covariance and fresh full mean/pressure estimates must be joined
quantitatively. Uniform physical-scale realization, additive cutoff/pressure
tails, complete residual improvement and blow-up remain unproved.

Reproduce with `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output /path/to/fresh.json`. Exact rational phase identities,
varying-growth jets and independent coefficient-sum bounds over a whole
unit box test the recipes on manufactured fields. The controls do not
evaluate the actual selected giant pulse or certify the continuum inputs.
