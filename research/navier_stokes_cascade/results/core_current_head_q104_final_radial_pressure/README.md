# Q104: complete final means and full shifted radial mean pressure

Parent: `f5c73ee64ef04ff45cf016745ccf94afc4f5df5b` (Q103).
Its CI runs 41479760558 and 41479760847 were still running when inspected
on 2 November 2026. Q102 field-comparison run 41465339945 had passed.
This checkpoint composes the FINAL six mean contributions and total wave
covariance with the frozen compact-source quotient calculus. It derives
conditional full mean-pressure/A0 coefficients and identifies the remaining
radial component obstruction. The quotient identities are inherited from
core_compact_source_quotient; they are not presented as new discoveries.

## 1. Assemble the actual final perturbation state before its source

Use exactly Q103.CurrentPostslow and its Q102/Q101 factories. The FINAL
perturbation means are the sums of

    initial fast + initial slow + first-real fast + first-real slow
                                      + sourced fast + sourced slow. (1)

All three components beta,v,gamma are retained in each summand. Old and
new axial fast fields include their actual gd-A1 gd velocities; slow fields
include their entire potential hull, including gaps. They are increments,
not replacements of older means. Unlike Q103's slow scalar reset, the
pointwise radial source needs FULL fields, not only auxiliary Haar means.

Take weighted initial-fast coefficients from the same budget._fast,
initial-slow coefficients from initial.slow_fields(weighted=True),
first-real fast/slow from CurrentScalars, and sourced fast/slow from
Q101/Q102. Embed EACH from its own pole r into p=max(all mean/wave poles)
using B(r,p)=max(1,Lhi)^(p-r), then add with safe 243 normalization padding.
The old operator recipes remain at their own proved poles; only coefficient
jets are embedded. All input bounds refer to the same original positive
norm hull, physical cover, root and flat width. No operator monotonicity is
inferred from the pole comparison.

Q97 supplies old and sourced complete-wave raw l1 coefficients A_j,V_j,
including +/-1 and +/-2 and all cut curls. Form U_j=A_j+V_j and

    W_n=sum_(j<=n) 5^j (U*U)_j/j!,                       (2)

where * is the raw binomial derivative convolution. This bounds EVERY
entry of the TOTAL symmetric angular covariance, not just Delta W.
Both complex cross orientations, all squares and the old covariance
remain. Its coarse grade is 1, from the old grade-1/2 wave; the higher
cross/sourced-square grades are conservatively embedded into it. The same
Gaussian and common-pole hypotheses are Q97's; P_i^2<=1 is used as there,
not inferred from unrelated overlapping phases. The wave is held fixed
during the mean updates.

## 2. Exact final radial source and local pressure memory

On this ONE final state, with Dr=partial_R+M d R^(d-1)vr.partial_y,
De=Dr+e/R, Dz=epsilon partial_Z and tstar=-epsilon partial_T+c Nt, put

    S=2b beta+beta^2+Wrr,
    H=b gamma+G beta+beta gamma+Wrz,
    C=2V v+v^2+Wtheta,theta,
    gr=-tstar beta-D1 S-Dz H+C/R
                             +epsilon(Dr D1 beta+Dz^2 beta). (3)

G beta remains on the FULL hull; the G=0 slow reserved-hull premise is
not extended to the whole source. The vector viscosity is really

    Dr D1 beta=Dr^2 beta+R^-1 Dr beta-R^-2 beta.           (4)

It must not be swapped for D1 Dr beta. Define

    r0=-tstar beta-S/R-Dz H+C/R+epsilon Dz^2 beta,
    K0=-S+epsilon D1 beta, gr=r0+Dr K0.                 (5)

Under the retained WHOLE compact smooth-extension and mixed-trace
identities, the SAME moving T0,A0 give

    P=int Haar(gr) dR=int Haar(r0) dR,
    pm=K0+T0(r0-rho P),
    Dr pm-gr=-rho P-A0(r0-rho P).                       (6)

The local -Wrr, nonlinear -S and viscous epsilon D1 beta pressure all
remain. A0 Dr K0=0 does NOT imply T0 Dr K0=0; the latter is K0.
No slow derivative is commuted through the moving inverse. P is the
Q103 postslow scalar of THIS state, not the earlier wave-only defect.

pm is the full angular-zero compact pressure response relative to the
retained exactly balanced base. Base Pi1/Pi2 and nonzero-angular wave
pressures keep their separate original ledgers; they are not declared
compact, zero or replaced by (6). Assigned base/exterior residuals remain
additive. The radial quantity in (6) is this pressure-balance component,
not an asserted formula for every term of the complete physical residual.

If only piecewise classical derivatives are known, let J be the sum of
transported jumps of K0, including exterior zero pieces, and Jless their
cumulative sum before the output point. Then T0 Dr K0=K0-Jless+chi J,
P=int Haar(r0)-Haar(J), and the last line of (6) acquires +chi_R J.
These signs follow the PLUS Dr K0 in (5), unlike Q98's minus-Dr-F split.
Quantitative bounds on every moving mixed trace are then required. Local
patch compactness or a scalar Haar check does not certify (6).

## 3. Explicit current-state coefficient program

Let X,U,Z bound the weighted final beta,v,gamma through N=n+2, W bound
the total covariance there, and wp=zeta delta^-p on [a,b]. Inflate each
mean coefficient by E=1+Eedge(p,width). Every mean-source monomial has
at least one compact mean factor; wp^q<=E^(q-1)wp makes feeding EX,EU,EZ
to the positive polynomial conservative, including squares. W is already
a SINGLE full weighted covariance and is not squared again. All fixed
b,V,G,c and graph multipliers require their original UNWEIGHTED whole-hull
jet bounds; safe 243 padding is retained. No weight inequality is
differentiated to claim derivative bounds.

Put d=n+1, dd=(n+1)(n+2), g=1+2 graph_N and
r=sum_(j<=n) a^(-j-1). For this paragraph X,U,Z denote the inflated
coefficients. S*,H*,C* are the positive polynomials of (3), with W for
each relevant tensor entry. A sufficient envelope is

    time=d(1+4c)X,
    visc=[dd(g^2+1)+r d g+r^2]X,
    Gr*=time+(d g+r)S*+d H*+r C*+visc,
    r0*=time+r S*+d H*+r C*+dd X,
    K0*=S*+(d g+r)X.                                   (7)

The original bounded Nt components are covered by the conservative factor
four; labels are held fixed. Through-N norms and normalized product rules
cover two coefficient derivatives and derivative jets of the graph.
Actual epsilon/graph grades are kept in the analytical ledger, then
beneficial powers are dropped in these positive coefficient bounds.

Rebuild the original e=0 compact recipe at EXACT p and output n. Let
Rho=243 B(0,p) weighted_rho(n), and let Q103 supply P*. Then

    Comp*=Rho P*, F*=r0*+Comp*,
    Pm*=K0*+KT0(n,p) F*,
    A0*=KA0(n,p) F*, Er*=Comp*+A0*.                     (8)

No ordinary Haar primitive or old lower-pole operator is substituted.
This is a same-state conditional normalized bound for the full compact
mean pressure and both radial response terms. All old pressure memory,
cutoffs and viscosity survive in the fields/source/operator programs.
The constants remain symbolic exponential-polynomials, not measured
canonical continuum suprema or uniform-in-stage bounds.

## 4. Honest coarse rate and the remaining radial inequalities

The complete final tangential/axial means retain coarse B=1-kappa and
radial B+1; later increments have higher grades. Total W has grade 1.
Thus gr and r0 have only the sufficient coarse grade B, due to Dr Wrr
and/or the centripetal 2Vv/R term. The local K0 has grade 1, but (8)'s
full mean pressure and direct A0 envelope are still only B. Q103 supplies
the stronger Comp grade C+1=5/2-kappa. It does NOT upgrade A0.

For h=r0-rho P, its whole Haar integral is zero, hence Haar(A0 h)=0.
The pressure-balance component Er in (6) therefore satisfies exactly

    Haar(Er)=-rho P, int Haar(Er) dR=-P.                 (9)

An improved grade-t bound on THIS component requires the corresponding
grade-t scalar/compensator bound, as well as control of A0 h. The
available upper bound at C+1 is not a lower bound or a failure theorem;
it simply does not prove a larger t. For the earlier illustrative complete
source target t=3-6kappa at kappa=1/100000, the gaps are

    t-B=2-5kappa, t-(C+1)=1/2-5kappa>0.                (10)

One still needs a justified inequality for the actual final shifted A0 h
at the chosen target, with its full high-jet/trace/domain program, and a
stronger scalar estimate or an accounted cancellation in the COMPLETE
residual ledger. A zero P alone cannot eliminate a nonzero-Haar A0 term.
The exact finite Fourier control exhibits that distinction; it is not an
actual canonical source lower bound. No complete residual improvement,
absolute physical decay, unconditional pulse bound or blow-up follows.

## 5. Requests and next obligation

Radial output n needs means/covariance through n+2 AND Q103 scalar n.
The latter and the fresh slow norm together require covariance n+11,
Q101 rows n+9, Q94 mean factory n+14, conservative target/profile
n+28,n+41 and existing new cutoff shape n+11. This new final-pressure
operator itself needs cutoff through n+1. Higher pressure gradients
require rebuilding the finite reserve; no extra derivative is silently
borrowed. Original compact tails and primitive DAGs remain unexpanded.

Next compose the FINAL angular/axial mean differences and nonzero-angular
regular-source ledger on this state, then address (9)--(10) and the actual
compact A0 bound. Native/trace/support premises, all-scale pulse estimates,
complete viscous/cutoff/interaction residual closure and infinite realization
remain unproved. Root, fields, width C_setup^-20, accepted band and frozen
historical records are unchanged. These are conditional norm recipes,
not actual field values or a mutation of the incoming state object.

Exact controls independently differentiate (3)--(5), retain the vector
connection and local viscous/covariance pressure, integrate a shifted
Fourier source with moving cutoff, verify (6) and (9), and detect a surviving
A0 at zero Haar moment. Polynomial fixtures have only finite endpoint
regularity; the analytical argument retains the stated smooth extensions.
Manufactured lower norm leaves verify all six mean stages, total old/new
covariance and exact operator domains, not actual canonical continuum data.
Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`.
