# Mixed compact/time inverse and the complete centrifugal source

Parent: `57aabbb712e34ce4e4a81363281a924f09f63308`.
Its compact-source quotient CI step was still pending when inspected on
23 October 2026; the workflow was running without a failed step.

The parent reduces the FINAL radial compact source to r0, but its order B
centrifugal term remains. Here we audit the additional structure of the
actual fast swirl map, rather than applying a generic scalar estimate
twice. For the ORIGINAL irrational directions the two divisors satisfy
an exact joint integer identity. It gives the sufficient mixed estimate

    N_m[A_e(a c^-1 Nt^-1 F)]
         <= c^-1 C_(m,p) M^-p N_(m+2p+2)(F), p>=1.    (1)

The five-derivative saving is relative to COMPOSING the existing separate
upper bounds, which reads F through m+2p+7. Relative to the ordinary
compact estimate on its own input, the mixed row-input count is one less.
A fast-map-compliant Pell comparison proves necessity J>=m+2p-1 for this
uniform finite-J row estimate. Thus (1) is within three derivatives of
that boundary; the enormous requests have not disappeared.

We also derive the exact nonlinear centrifugal rewrite, retain its Haar
variance in pressure, and split the remaining radial source at orders B
and 1. These are conditional analytical estimates and exact identities.
Actual all-scale field/row/high-jet/support inequalities remain unverified.
The full regular source stays E=2-3kappa. No absolute decay, unconditional
pulse, infinite field or blow-up is established. Historical evidence,
fields, roots, width C_setup^-20, carriers, cutoffs, h and accepted band
are unchanged.

## 1. The two ORIGINAL divisors are conjugates, not independent losses

Put b=sqrt(2)-1, vr=(1,-b), vt=(b,1). For nonzero integer k let

    alpha=vr.k=(k1+k2)-sqrt(2)k2,
    tau=vt.k=b[(k1+k2)+sqrt(2)k2].

Consequently

    alpha*tau=b[(k1+k2)^2-2k2^2].                    (2)

The bracket is a nonzero integer: a zero would give a rational sqrt(2),
unless k=0. Thus |alpha*tau|>=b. With r=|k|_1 and the original bound
|alpha|^-1<3r, for integer p>=1,

    |alpha|^-p |tau|^-1
        <=b^-1 |alpha|^(-(p-1)) < 3^p r^(p-1).      (3)

Both divisors MUST refer to the same Fourier mode. An auxiliary-dependent
multiplier can mix modes, so cannot be pulled through Nt^-1 for free.
The fixed multiplier a in (1) is auxiliary-independent. In the actual
nonlinear application all auxiliary-dependent factors remain inside the
new source Gtheta; (3) is applied to its OUTPUT mode, not separately to
two unrelated internal factors.

## 2. Explicit raw kernel bound, including all fixed profiles

Use the parent's U=R^d kernel on the same fixed positive input hull:

    A_e f=k_e chi_U int_all f_e(U+a')
                 f(U+a',y+Ma'vr,Z,T) da'.

Let F have zero auxiliary Haar and smooth whole-line radial extension.
All coefficient derivatives hold M,c fixed. Define J=m+2p+2,
B_m=N_m(k_e chi_U), A_J=N_J(f_e a), and retain the ORIGINAL hull length L.
A valid sufficient constant in (1) is

    C_(m,p)=2^(m+J) L B_m A_J [8*3^p*2^(m+p+2)].     (4)

In a derivative with a_y torus factors, p radial integrations give M^-p,
and the temporal multiplier adds (2pi i tau)^-1. Bound the product of
divisors by (3). Fourier integration by parts a_y+p+2 times in a coordinate
with |k_j|>=r/2 leaves r^-3. There are 4r lattice points with l1 size r,
and sum r^-2<=2, giving the factor eight. The input uses

    (m-a_y)+p+(a_y+p+2)=m+2p+2 derivatives.

The leftover 2pi power is (2pi)^(-2p-3)<=1. Leibniz on the output profile
and on (f_e a)F gives 2^(m+J). Smooth compact extension removes every
boundary term. This proves (1) for the complete coefficient derivatives,
including U,Z,T and both torus coordinates; no M derivative is taken.

For the radial application a=R^-1, so (4) needs N_J(f_0/R), not just
N_J(f_0). Here the ACTUAL new profile simplifies exactly:

    f_0/R=R^-d/d=1/(d U),
    partial_U^j(f_0/R)=(-1)^j j!/(d U^(j+1)),
    N_J(f_0/R)<=J! max(1,u0^-1)^J/(d u0), U>=u0>0.   (5)

The SAME original positive hull floor u0 and d are retained. This supplies
an explicit analytic bound for this fixed kernel factor at every requested
jet, without numerically allocating its giant factorial. It does not bound
the unknown actual source seminorm. All non-U derivatives of this profile
vanish. Qcirc=I-Haar costs a separate factor at most two
when a source is not already supplied with zero Haar. Raw-to-factorial
conversion keeps J! as a deferred finite factor. Weighted source
unweighting, output cutoff-support edge floor, original width/gap,
c^-1<Tg S, and M^-p<(Tg S)^(rho_g p)epsilon^(p kappa) are unchanged.
Every resulting S polynomial is absorbed with the SAME strict-margin
rule on the same accepted band. No giant constant is numerically enclosed.

The separate temporal estimate costs four derivatives and the separate
compact estimate costs m+2p+3. Applying both literally would therefore
request m+2p+7 derivatives of F. Equation (4) avoids five of those, but
does not claim either separate theorem was false or minimal.

## 3. Sharp comparison which actually obeys the fast swirl map

Take the parent's Pell modes P+Qsqrt(2)=(1+sqrt(2))^j, k=(P-Q,Q).
Their joint integer in (2) is (-1)^j. For any large M choose the first
P^2>=M, giving sqrt(M)<=P<3sqrt(M). If a is nonzero on a compact input
patch, choose a fixed smooth unit-mass radial bump phi there and a fixed
compact slow bump Bslow. Set

    F_M=A (f_e a)^-1 phi Bslow
                   cos(2pi[k.y-M alpha U]).           (6)

The SAME drift cancels the radial phase. Nt^-1 changes the cosine to
sin/(2pi tau). Since |tau|<=P and Q>=P/2,

    N_J(F_M)<=|A|16^J C_J P^J,
    N_m[A_e(a Nt^-1 F_M)]>=c_* |A| P^(m-1).          (7)

Here fixed positive c_* includes the cutoff/output profile and harmless
2pi constants. Taking A=P^-J yields a ratio to M^-p N_J(F_M) growing
like P^(m+2p-J-1). Therefore

    J>=m+2p-1                                      (8)

is necessary for this uniform finite-jet estimate. Equality and removal
of the upper rule's extra three derivatives are not settled.

For e=0 and a=2V/R on a positive background patch, set the actual FAST
swirl v=-c^-1 Nt^-1 F_M, beta=gamma=Psi=0. This comparison is exactly
divergence free, satisfies both weighted Haar velocity moments, and
obeys the original fast formula. Its centrifugal linear source is
2Vv/R; the v^2 square has only zero/double Fourier modes and cannot
cancel its fundamental. The pressure's Haar variance is retained.
Thus the fast-map relation itself changes the derivative boundary by
ONE, not by an arbitrary amount. The prescribed F_M is not the actual
canonical prefast residual. Higher derivatives of the finite-J family
grow, so it is not a counterexample in the original all-fixed-jet class,
nor a counterexample to the construction or a lower bound on its residual.
Any fixed polynomial in log M cannot overcome a strict deficit in (8).
With amplitude epsilon^a P^-J and the original M flooring, the lower
epsilon power is a+kappa(J+1-m)/2, times fixed S powers and c^-1.
A strict deficit relative to a requested target beats every S polynomial;
equality is left open.

## 4. Rewrite the FULL accumulated centrifugal source, not one fast piece

Work at the FINAL stage after the latest fresh means. Put

    Fsum=-c Nt v, Gtheta=2(V+v)Fsum, Cm=2Vv+v^2.

The background V is auxiliary-independent. Every slow swirl increment
has Nt equal to zero. Every fast swirl is its original -c^-1 Nt^-1 row.
Within one common chart their residuals sum to Fsum. If old label/cover
representatives are used, keep the ORIGINAL intertwiners and temporal
normalizations; do not simply add differently normalized raw rows. An
alternative using only current coarse inputs is the exact definition
Fsum=-c Nt v, costing one current v derivative. No new desired fine-row
callback is supplied by that definition.

Since Nt Cm=-c^-1 Gtheta, and Nt^-1 Nt=Qcirc,

    Qcirc Cm=-c^-1 Nt^-1 Qcirc Gtheta.               (9)

v in Gtheta is the FULL accumulated mean, including all slow/fast parts.
Pulling it through Nt^-1 would omit the product commutator and the factor
two in Nt(v^2). Retain every old/new product. The exact Haar part is

    Haar Cm=2V Haar(v)+Haar(v^2),                    (10)

with the full fast variance and old/fast correlations. It is neither zero
nor generally 2V Haar(v)+Haar(v)^2. Formula (9) is valid ONLY for Qcirc;
using it for Cm would corrupt the current pressure and moments.

With the parent's S,H,K0 and SAME current stage, let

    Rrest=-tstar beta-S/R-Dz H+Wtt/R+epsilon Dz^2 beta.

Then the EXACT radial tail becomes

    A0 Qcirc gr=-c^-1 A0[R^-1 Nt^-1 Qcirc Gtheta]
                                      +A0 Qcirc Rrest.          (11)

The inherited coarse inputs imply Gtheta in M_B and Rrest in M_1:
the latter's groups have orders B+1,1,min(B+2,2B+2,2),1,B+4.
`split_envelope` gives their positive coefficient polynomials using
Fsum=-c Nt v with one extra v jet and beta through two extra jets.
It retains the original fixed background multiplier norms, rather than
asserting those backgrounds lie in a compact weighted class.

The complete pressure reconstruction is still

    r0=Haar(Cm)/R-c^-1 R^-1 Nt^-1 Qcirc Gtheta+Rrest,
    P=int Haar r0, pm=K0+T0(r0-rho P),
    K0=-S+epsilon D1 beta.                          (12)

Keep the moving unit-mass rho, full six-entry current wave covariance,
axial velocity/pressure moments, local K0 and all old pressure memory.
The balanced base's Pi1/Pi2 and assigned exterior residual are not removed.
No new slow derivative is commuted past a moving cutoff. In particular
(11) is not an improved estimate for pm or the entire residual.

## 5. Specific requests and the unchanged actual-input frontier

For target C0 or the parent's physical-route target
(N+3/2+2h+2m)/h, choose strict integer counts separately:

    pC=floor((target-B)/kappa)+1,
    pR=floor((target-1)/kappa)+1,
    JG=m+2pC+2, JR=m+2pR+3.                         (13)

All displayed targets have pC=pR+1>=1. At kappa=1/100000:

| Request, output m=0 | Gtheta raw jet | Rrest raw jet | Previous full gr raw jet |
|---|---:|---:|---:|
| Normalized target C0 | 200002 | 200001 | 200003 |
| Physical route h=1/200,N=0 | 60200006 | 60200005 | 60200007 |

These refer to DIFFERENT functions, not a free reduction of total original
source complexity. If Fsum is bounded generically through Nt v, Gtheta at
JG needs v through JG+1. Rrest at JR needs beta through JR+2. Their largest
LOCAL current field request is exactly the previous full-gr raw index,
where its old conservative program requested fields three orders further.
No older source/root/normal DAG has been expanded or a total minimal
original-profile index proved. Keep the ENTIRE older graph conservatively.
Using an independently derived same-prefix prefast-row sum may save costs,
but that sum's own norms and chart intertwiners must first be proved.

Gtheta and Rrest can be embedded at final pole 12q under the retained
unweighted fixed-multiplier contracts. The original A0 raw source
unweighting/output-edge floor and archived shifts remain; (4) additionally
retains the f0/R reciprocal-profile derivatives. If a fixed multiplier has
an edge pole, add it at each product before choosing q. All new local
requests fit within the older conservative graph; no lower-pole operator
monotonicity or new physical band is assumed.

The still-required inequalities are the ACTUAL all-band weighted/raw
Gtheta and Rrest norms at JG,JR, with complete field identification,
whole-hull smooth extension and all original profile/base/support bounds.
Manufactured norm data do not prove them. Equation (8), or its epsilon
version J>=ceil(m-1+2(target-a)/kappa) for a row of order a, also rules
out replacing those costs by arbitrary low finite row norms. It does not
prove failure of the actual canonical inequalities.

Next validate these actual reduced inputs and their support/pressure
identities, then advance the next generic signed cycle. The parent's
direct A1 order B+1, full regular order E and negative physical exponent
are unchanged. A controlled infinite diagonal, localization/forcing
regularity and complete velocity growth remain open before blow-up claims.

Exact checks use the joint integer form, signed Fourier/time normalization,
radial integration by parts, both resonant directions, the full nonlinear
Haar split and pressure variance, and independent epsilon norm budgets.
Finite polynomial fixtures are not smooth canonical fields or continuum
enclosures. Reproduce with `python test_audit.py -v` and `python audit.py
--verify-record evidence.json --output FRESH_PATH.json`.
