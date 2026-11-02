# Q96: zero-entry normalized response, full pressure and complete cut curl

Parent: `fa4ace3740260083fc0f34a2a05258a8d07e9702` (Q95).
This checkpoint defines the adjacent sourced response by Duhamel and derives
explicit positive jet bounds from Q95's SAME first-real source. All conclusions
are CONDITIONAL on the original native input, normal/frame/propagator and
support contracts. The program computes norm polynomials, not numerical
selected fields. Historical evidence, widths, fields and labels are preserved.

## 1. Normalize the equation without dropping the moving constraint

Use core_full_path_ode_reconstruction (1)--(9). Set N=|n|^2, h=n/N,
H=h n^T, Pi=I-H, p=P'/P, lambda_m=m^2 epsilon k^2 N and
A0=-K+h(n^T K-(n')^T). The full unprojected source is f, and g=f/P is
precisely Q95's normalized source, at B=1-kappa. Put u=t/P. Then the
actual tangent zero-entry response satisfies

    u'=(A0-lambda_m I-p I)u-Pi g, n.u=0, u(0)=0.       (1)

For estimating mixed derivatives use the ambient extension

    C=A0-lambda_m I-p I+(lambda_m+p)H.                 (2)

It agrees with (1) on tangent vectors and (n')^T+n^T C=0. Thus n.u is
constant along every homogeneous ambient trajectory. Differentiated
solutions need not be tangent: n.Du=-Dn.u. Applying a tangent propagator
to arbitrary differentiated forcing would be unjustified.

The earlier ForcedBudget bounds RAW t derivatives divided by P; Q95 bounds
DERIVATIVES of f/P. Feeding these into an old recurrence and calling its
output derivatives of t/P would omit the -p I in (2). Q96 uses (2) and
later restores P by the full Bell product, keeping both conventions explicit.

## 2. Ambient propagator and normalized response jets

Retain the actual tangent bound ||V_T(v,w)||<=T P(v)/P(w), T>=1, and length
<=L, on EVERY subinterval of the same original whole path. For y'=Cy,
z=n.y is constant and y=v+h z with v tangent. Then

    v'=(A0-lambda_m I-p I)v+(A0 h-h')z,
    Kamb=T[1+L Q0 n0]+h0 n0, Q0=A0_0 h0+h1.            (3)

The normalized tangent propagator is bounded by T. Duhamel proves the
ambient bound Kamb, without exp(L*coefficient), a new smallness condition,
a reciprocal Gaussian minimum or a new smaller-q selector. Normal jets
through one extra order are needed for h1 and n'. Smooth coefficients and
source give existence/uniqueness by the finite-interval Picard factorial
series. The zero-entry field is the unique Duhamel solution; no entrance
amplitude or positive square root is replaced.

For raw mixed-jet norm sequences use binomial convolution *. The normal,
moving-normal and K jets come from the SAME Q93 phase provider, transported
to the whole path as in Q95. The reference rate p_j is a separate original
whole-path UNWEIGHTED bound. If only log P jets are available, p jets at
order r require log P through r+1. This costs reference jets, not new source
jets. A provider with coefficient edge poles requires additional pole
propagation; unweighted coefficient hypotheses cannot be replaced by a tag.

Use the reciprocal-square recurrence from Q93, then h=n*Ninv,
v=n*K+nprime, A0bound=K+h*v, H=h*n, damping=4 m^2(n*n). A sufficient
coefficient bound for (2) is

    Cj=A0bound_j+damping_j+[(damping+p)*H]j+p_j.        (4)

The last p_j bounds the normalized -p I term. Fixed labels keep epsilon k^2
constant, <=4; a chart varying it needs its full jets. Pi0 has Euclidean
operator norm one; higher Pi jets are bounded by H jets. The component-max
Q95 g bound is safely doubled for Euclidean vector norms; let F=Pi*G.

For each j, define Uj as the sum of

    Kamb L[Fj+sum_(i=1)^j binom(j,i) Ci U_(j-i)]

and, for j>=1,

    F_(j-1)+sum_(i=0)^(j-1) binom(j-1,i) Ci U_(j-1-i). (5)

For parameter-only jets use ambient Duhamel. For a jet containing a pulse
 derivative, differentiate (1)--(2) one fewer time. Adding the two positive
bounds covers every mixed jet. This is triangular, uses source jets only
through j, and proves |D^j u|<=epsilon^B sqrt(zeta) delta^-q Uj. All
coefficient jets in (4) must actually be unweighted; Q95's optional log
poles are separately tracked in q and the restoration below. No uniformity
in j, harmonic range or construction stage is asserted.

## 3. Full pressure and potential, then restore the Gaussian

The exact local pressure and potential are

    pi/P=i[(Ninv*vrow).u+h.g]/(k m),
    A/P=i h cross u/(k m), vrow=n^T K-(n')^T.           (6)

The positive bounds use v=n*K+nprime and yield

    Pressure_j=|m|^-1[(Ninv*v)*U+h*G]j,
    Potential_j=|m|^-1[h*U]j,                          (7)

at order B+1/2, since k^-1<=sqrt(epsilon). Full G is retained: a purely
normal source can have zero velocity and nonzero pressure. Neither the
moving-normal term nor n.f is dropped. This local wave pressure does not
replace Pi1/Pi2 background memory, cumulative mean pressure, axial gd-A1 gd
or the differentiated c_rho P moment in the subsequent mean update.

Restore RAW field jets before taking spatial curls. If Hj is Q95's Bell
bound for D^j P/P, each raw t/P, pi/P and A/P jet is bounded by H*U,
H*Pressure and H*Potential respectively. The common pole increases from
Q95 q by the extra Bell pole when log jets carry poles. With the original
unweighted log contract it stays unchanged.

## 4. Evaluate the packet and cut its complete potential and pressure

Whole-path response jets are not automatically physical packet jets.
Evaluate w=v(H) on the SAME original fixed lift. This is affine; the
nontrivial column sum is bounded by 1+376 c_i^-1, using ||J^Delta||<=188
and ||lambda_t||<=2. The original compact archive supplies I0>=c_i^-1/S,
so a sufficient bound is

    E=1+376 S I0.                                     (8)

Multiply each raw order-j response coefficient by E^j. No extra source
derivatives or Gaussian ratio appear: P is evaluated at the same point.
Moving labels, entrances or a different nonaffine chart require their own
chain terms. This return to the five native coefficient coordinates is
separate from Q95's six-coordinate redundant path norm.

Use the original cutoff jets chi on that full native R chart, cutting BOTH
potential and pressure. Raw cut jets are chi*(restored evaluated jets),
including every derivative of the cutoff. The exact complete curl is

    curl(Acut exp(i k m Phi))=chi t exp(i k m Phi)
                                  +curl_* Acut exp(i k m Phi).      (9)

Coefficient curl contains -Dz Atheta, Dz Ar-Dr Az and (Dr+1/R)Atheta.
At output order j it requires potential jets through j+1. A sufficient
raw l1 bound is two Dr products plus two Dz derivatives plus the reciprocal
radius product. The graph archive's factorial-through norms give raw
coefficient jets i! g and i! r0. The resulting curl remainder has order
B+1/2-kappa, and complete cut velocity has order B after dropping only
its positive reserve. The program returns principal, curl, full velocity,
pressure and potential coefficients separately.

Support follows from whole slow/transverse paths and uniqueness, not one
entrance value. Q90's entrance-flat cutoff grammar remains an actual
premise; no raw exit zero is inferred. The outgoing potential and pressure
cutoffs and all their residuals remain. Magnitude bounds here do not
establish physical corner extension or a complete cutoff/residual estimate.

## 5. Scope, finite request and next obligation

For complete cut output through j, use Q95 source ceiling R=j+1. Q94's
mean-factory reserve is j+3, with its conservative target/profile requests
j+17,j+30 before full primitive expansion. Reference p jets through j+1
and cutoff jets through j+1 are additional original providers. The phase
moving-normal branch requests one extra normal jet; Q94 already reserves
more through its wave pressure chain. I0, T and L are ORIGINAL same-domain
bounds, not independently certified numbers or altered lengths.

Established: conditional Duhamel existence and explicit normalized, restored,
evaluated and complete cut response estimates for Q95's first-real source.
All actual native/source/frame/propagator/rate/cutoff inequalities remain
hypotheses. Next form the changed full complex wave covariance, all six
tensor entries, both cross orientations and the new square, and recompute
pressure, moments, fast/slow means and stage-specific compact sources.
The earlier first-real tails cannot be transplanted. No unconditional pulse
family, complete residual improvement, infinite realization or blow-up is
claimed. Exact manufactured checks and frozen provenance are not continuum
proofs. Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`.
