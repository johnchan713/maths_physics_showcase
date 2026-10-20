# Fourth scalar/moving loop: explicit conditional cutoff-shape gate

Parent: `66ff8c49955a1b17e04a617d42681633cc19b51b`.
Status: `fourth-angular-signed-root-moving-phase-and-conditional-shape-gate`.

This checkpoint bounds q through order FIVE, the regular signed variance
root through order four, and the complete fourth moving-phase composition.
With the parent's actual native pressure/mixed inputs, it proves a
conditional fourth-angular primitive bound below the unchanged H. It also
shows conditional headroom in the earlier compact angular-C4 inequality,
retaining every five-row normalizer and the same exact small root.

The remaining cutoff condition is stated explicitly below. It is NOT
verified for the actual variance cutoff: its plateau locations and smooth
range do not provide the derivative budget. No cutoff is replaced by the
explicit RADIAL step. No unconditional compact acceptance, global mixed
leading C4, first C3, second H2,3, physical pulse, covariance/mean, complete
residual or blow-up claim follows. These are conditional analytical
derivations; finite diagnostics are not a continuum proof or independent
foundation review. Historical evidence and all roots, gauges, supports,
width C^-20, C, G, H and N=1+floor(H^32) remain unchanged.

Sources: the retained [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(C.9)--(C.16), and the transcribed maps in `loop_modulation_bounds`,
`modulation_repair_all_order`, `fourth_pressure_input_boundary` and
`fourth_target_normalization`. PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The manuscript's final theorem is not a premise.

## 1. Inputs and the precise unproved condition

Set A=C^4096>=2^16, Bn=A^256, G=exp(A^32), H=exp(G^4),
deltaL=exp(-A^20), d0=1/(16A). All slow derivatives in this note
are ANGULAR eta derivatives at fixed y=log X. We do NOT claim new
joint y/eta/phase order-four bounds. Norms sum suprema of normalized
derivatives; angular endpoints use one-sided limits. Uniformity in y
is on the original pre-modulation compact rectangle.

The preceding native checkpoints give p=p2 in C4 below A^1546,
a,bs,ts,vs in C4 below powers at most A^2583, and E in C4 below
A^2. The original value lemmas give positive floors, a>=1/A,
2<v<=A, and mu,|mu p|,|t|<G. Their pointwise cap proof uses only
the cutoff range, not a derivative budget. Since A^10000<G, all the
needed native C4 norms are below G. In particular the new pressure
envelope is NOT silently assigned the old Bn. Full axis and radial
pressure memory remains in p2; p1 is not another varying argument of F.

For the SAME actual target

`s=zetaL(vs)*sqrt((2+deltaL/2-vs)/a)`, extended flatly by zero,

the additional condition used here is

`(S)  |s|Ceta4 <= G`.

On the nonzero support the numerator is at least deltaL/4. The prior
shape audit defines Kshape(t)=zetaL(2+deltaL/8+(deltaL/8)t) and
Z=max(1,sup|Kshape^(j)|,1<=j<=4). Its justified bound is

`|s|4 <=1e12 Z deltaL^-8 A^20673`.

Thus the precise sufficient inequality still to verify is

`(Z)  log Z +12log10+20673log A+8A^20 <= A^32`.

For example Z<=1e100 would suffice, but is an AUXILIARY sufficient
budget, not an enclosure of the retained cutoff. The compressed smooth
monotone cutoffs in the earlier checkpoint show that plateau/range
conditions alone cannot supply (Z). A separate actual-composition bound
proving (S) would be an alternative. Every conclusion below is conditional
on (S), not an assertion that the selected profile fails it or satisfies it.

## 2. q5 and the removable signed variance root

Write q(z,theta)=exp(z sin theta)/I0(z), |z|<=G. Positivity of the
integral for I0 gives |I0^(j)|<=I0 and I0>=exp(-G).
The normalized reciprocal coefficients satisfy b0=1/I0 and
|bj|<=2^(j-1)/I0 for 1<=j<=5. Their norm sum is <=32exp(G).
For the numerator, the nonconstant joint Taylor polynomial of
z sin theta through total degree five has norm <=40G. Hence its
exponential jet is bounded by exp(G) sum_(m=0..5)(40G)^m/m!.
The universal comparison

`1e12 G^10 <exp(G)`, for every G>=2^16,

proves the uniform joint (z,theta) bound |q|5<exp(4G). It follows
coefficientwise from differentiation under the integral that

`R(z,theta)=integral_0^1 q_z(rz,theta)dr`, `|R|4<exp(5G)`.

The factor from q_z is at most FIVE in normalized order-four jets.
Using only q through four would miss the highest input. No division by
z or p is made. For R2(z)=average_theta R(z,theta)^2,

`R2>=exp(-2G)/2`, `|R2|4<exp(10G)`.

The lower bound follows from the positive Bessel numerator series,
whose constant after removal of z^2 is 1/2, and I0^2<=exp(2G).
Finite square-root Taylor expansion gives |sqrt(R2)|4<exp(60G).
For example its value costs exp(5G), inverse base R2 costs <exp(3G),
and four powers of the normalized increment cost at most exp(52G).

Use the SIGNED extension

`F(mu,p)=d0 mu sqrt(R2(mu p))`.

At any actual base point, the joint (mu,p) Taylor increment of mu p
has norm <=3G. Thus |F|C4_partial<exp(64G). The retained tilted
variance lower bound gives Fmu>=d0 exp(-4G)/32, hence
Fmu^-1<exp(10G). At mu=0 use EXACTLY Fmu=d0/sqrt(2), Fp=0.
Neither Vmu=0 nor p=0 is inverted. This proves the quantitative
fourth partial estimate that the old third-order ledger did not supply.

As an exact diagnostic of the removable formula,

`[z^5]q=sin(theta)^5/120-sin(theta)^3/24+3sin(theta)/64`,
`R2(z)=1/2-3z^2/32+11z^4/576+O(z^6)`.

Dropping q5 changes the displayed fourth R2 coefficient. The positive
series and bounds, not this finite expansion, justify the global estimate.

## 3. Fourth signed-root jets from the actual p4 input

At fixed y and a base eta, use normalized jets. Remove only the
unknown mu_n and take the nth coefficient of F(mu_<n,p)-s;
divide by Fmu at the base point. This retains Fp p_n and all lower
cross terms. The exact positive Taylor count has respectively
2,5,12,28 terms (including the target) for orders 1 through 4.
Every coefficient count is absorbed by exp(G). Bounds of the form
exp(cG) therefore close as follows:

| n | largest lower input exponent | sufficient c | chosen c for |mu_n| |
|---|---:|---:|---:|
| 1 | 1 | 76 | 80 |
| 2 | 160 | 235 | 240 |
| 3 | 320 | 395 | 400 |
| 4 | 480 | 555 | 560 |

Here the sufficient exponent is 10+max(1,64+largest)+1.
Summing all five normalized jets gives |mu|4<exp(600G).
This is a fresh conditional C4 estimate, not an invocation of the old
unbudgeted cutoff derivatives or a C3-to-C4 promotion.

## 4. Complete fourth moving inverse and primitive bounds

Before reparametrization, t=ts+d0 mu R(mu p,theta). Composition of
the preceding finite jets gives a joint (eta,theta) norm through four
below exp(4000G): the R composition costs at most
exp([5+4(600+1+1)]G), then multiplication by mu and addition of ts.
Keep the sharper pointwise |t|<G separately.

The EXACT variance equation F=s gives
v=vs+a s^2, so |v|4<=2G^3<G^4; its value remains 2<v<=A.
Finite reciprocal jets give |1/v|4<G^17. For

`w=a(1+t^2)/(2pi v)`, `Phi(eta,theta)=integral_0^theta w(eta,r)dr`,

we obtain |w|4<exp(10000G) and |Phi|4<exp(10010G) on
one fundamental interval 0<=theta<=2pi. In the latter bound pure
angular derivatives cost at most 2pi and positive theta derivatives
come from w. Value floors give G^-3<w<G^3.

The mean-one identity is exact: average R=0 and F^2=d0^2 mu^2 R2
give average a(1+t^2)=v. Therefore integral w=1, and the SAME
inverse circle lift satisfies Phi(eta,Theta(eta,varphi))=varphi.
Apply the normalized implicit recurrence, with divisor w and
partial budget exp(10010G). Safe inverse angular jet exponents are
11000,33000,55000,77000. The complete norm of its positive-order
jets is <exp(80000G). Every lower cross/phase partition is included;
Phi40=0 alone would not remove its mixed fourth remainder.

The two rational scalar factors 1/(1+t^2) and -t/(1+t^2) have
normalized derivatives of every order bounded by one on the real
line, from their partial fractions at +/-i. Thus, BEFORE moving the
angle, the sources

`gA=-(v/(1+t^2)-a)/2`, `gB=E*(-v t/(1+t^2)-bs)/2`

have joint (eta,theta) C4 norm below exp(20000G). Full finite Taylor
composition with Theta, not a frozen angle, gives

`|gA(eta,Theta)|4, |gB(eta,Theta)|4 <exp(350000G)`.

Indeed 20000+4(80000+1)<350000. This bounds the EXACT fourth
moving composition in `fourth_pressure_input_boundary` (7), including
all its R4 terms. The lift value is bounded on one period only;
positive-order slow jets and the pulled periodic sources extend to all
phases. Their phase means vanish exactly: integrating aL with w gives
a, while integrating bL gives -a ts=bs. This is a leading loop identity,
not a claim about complete physical-wave covariance corrections.

Zero-mean periodic integration on period one has operator norm <=2
and commutes with angular differentiation at fixed EXTERNAL phase.
Consequently the actual angular C4 primitive norms obey

`|Acal|4, |Bcal|4 <exp(400000G) < H`.

The last comparison is 400000<G^3. This supplies a4,b4<H under (S).
It does not bound fourth y derivatives or new full joint J4 norms.

## 5. Conditional compact angular-C4 gate at the SAME N and root

No old quantitative joint-C3 loop budget is used to obtain these angular
bounds. They themselves supply the angular C3 H premise of the tame
fourth-density argument. Use its COMPLETE five target rows, including
Cp, e(eta)=e_star/(1+eta^2), every inverse-e derivative and the
inverse-lambda transformed M-J row. The original value geometry gives
rhop^-1,e_star^-1,lambda^-1<=A; take the valid auxiliary square-root
enclosure s_rho=A^2. This changes no physical scale or patch.

The five C4 reciprocal normalizers are <=A^5. Every raw moment
affine coefficient is <=A^515, using Bn=A^256. The transformed
row coefficients therefore satisfy

`ZH,Za,Zb <= A^522`.

These are positive-coefficient bounds on the exact `target_affines`
map, not omitted normalizations. Its full target C4 norm is bounded
by 3A^522 H/N, and hence by H^16/N. The constant-Q repair map
has inverse<=1000, quadratic<=1e5; its C4 contraction is valid at
N>H^32. It gives the same prescribed small root by C1 uniqueness,
with norm <=2000H^16/N=:r. This justifies r angularly without
borrowing unverified fourth or mixed cutoff budgets. All five exact
functional matches, and therefore their leading future pressure memory,
are preserved. The separate positive-order Pi1/Pi2 tails are not erased.

For the earlier compact interface, let K=H^64 as before. Its known
part is Kknown=Bn+6400Ar+960000000000Ar^2. Since r<1,
Kknown<A^259<G. The three coefficients
alpha=3200000A Za, beta=1+3200000A Zb,
ch=1+3200000A ZH are each <A^525<G.
Substituting a4,b4<H into the exact positive majorant gives

`compact_cost <=Kknown+(ch H+alpha a4+beta b4)/N`
`                 <G+3GH/N <G+1 <H <K`.

Equivalently 2alpha L_A+2beta L_B<=R has headroom using the actual
source bound above. This proves a CONDITIONAL angular compact-input
certificate under (S), not verification of (S). It does not establish
the separate radial/state cone acceptance conditions, FINAL full mixed
leading bounds, actual first C3, second H2,3, retained physical scales or
complete pulse/covariance/residual estimates.

## 6. Meaningful checks and next obligation

Run `python test_audit.py -v`, then `python audit.py --verify-record
evidence.json --output FRESH.json`. Frozen history is byte-bound.
Universal parameter comparisons are checked at 80/110 digits at the
minimum A,G; monotonic ratios prove the all-larger statements. No actual
G,H,N or enormous flow is materialized. Auxiliary algebraic target tests
call the ORIGINAL five-row affine/gate implementation with nonzero
normalizers and inverse lambda, not a rewritten simplified substitute.

Exact checks detect omitted q5 and frozen inverse phase. The implicit
phase control gives fourth coefficients 7/81,135/4096,114/15625 at
eta=-1,0,1 for the manufactured lift theta+eta sin(theta)/2.
Numerical signed-root diagnostics compare normalized fourth jets with
nine independently solved root values at zero pressure and zero variance.
Their tolerance is 1e-18; raw version-dependent roundoff is not frozen.
These are diagnostics, not continuum enclosures of the selected profile.

The next obligation is to identify an actual derivative enclosure for
the UNCHANGED variance cutoff satisfying (Z), or independently prove
(S). Then join the remaining radial/state/FINAL leading and first-
coefficient inputs and retained-scale gates. The missing premise cannot
be supplied by passing tests or by the radial step's explicit formula.
