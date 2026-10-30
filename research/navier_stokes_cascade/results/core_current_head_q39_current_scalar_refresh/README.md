# Q39: current postfast scalars and the same slow response

Parent: `4d10116d38f4fd193eb42c5d94078b466767c0a8` (Q38).
Q38's focused CI 41134964277 passed; broad CI 41134964307 was still
running at inspection on 30 October 2026. This checkpoint connects its
complete current pressure/scalar state to the next fast scalar change and
slow five-map. It derives the OLD full mean norms from the same initial
constructor. The NEW actual fast-field estimate remains an explicit lower
input. It is not replaced by an assumed current scalar or a fine source.
Historical frozen evidence and the physical construction are unchanged.

## 1. Reconstruct the full old means on the whole hull

The real signed wave changes covariance and pressure, but no mean velocity.
Thus the old means are precisely the Q35 initial fast plus slow fields.
For normalized through-order n norms, unweight the old fast field at its
own original finite pole p0 using Eedge(p0,w), and add the already unweighted
slow field:

    O_k,n=Eedge(p0,w) F_initial,k,n + S_initial,k,n.            (1)

The actual slow labels X,U,Z map to radial, tangential and axial x,u,z.
The slow field is evaluated on its entire potential hull. The fast field
includes the actual axial gd-A1 gd term; it is not just gd and it is not
replaced by its Haar average. Its zero Haar average does not make it zero.
All original fixed-map support and pressure premises remain required.

Let X,U,Z bound the NEW complete fast increment through n+1, with declared
finite pole v and the actual weighted source contract. Unweight these at
Eedge(v,w), separately from (1). No monotonicity shortcut between old and
new operator constants is used. For the full common radial norm hull [a,b],

    M=(b-a)(a^-1+1+b+b^2)                                    (2)

bounds all required unweighted radial integrals (weights -1,0,1,2), with
normalized Haar measure. Both inputs must obey their norms on this entire
hull. The interval is not reduced to the small reserved bumps.

## 2. Current pressure already contains the signed change

Q38 returns Theta_signed=Theta_initial+Delta_wave for P,Jtheta,Jz. Add only
the actual fast change. For old complete means O_x,O_u,O_z and new X,U,Z,
the historical positive functional is

    S=2 O_x X+X^2,
    T=2 O_u U+U^2,
    H=O_x Z+O_z X+XZ,
    Delta_fast P <= M[S+T+(n+1)H],
    Delta_fast Jtheta <= M[O_z U+O_u Z+UZ],
    Delta_fast Jz <= M[2 O_z Z+Z^2]+(M/2)[S+T+(n+1)H].        (3)

Inputs through n+1 cover the derivative of the radial-axial flux. These
are sufficient norm bounds, with signs replaced by positive terms. Exact
scalar addition retains its original signs. The output is

    Theta_postfast = Theta_signed + Delta_fast.               (4)

Adding a second wave_change to (4) would double count Q38's update. Using
only the Haar old means in (3) would instead omit correlations: for old
cos(t) and new 2cos(t), the new field has zero Haar, but the averaged square
change is 4, whereas its own square averages only to 2. This elementary
control shows why the full old fast field must be retained. Exact zero
Haar cancels the prescribed auxiliary-independent base linear terms, not
products with old nonzero-Haar fields.

For the first real update, B=1-k, alpha=2-3k and D=2-4k. Assuming the new
fast tangential/radial orders D,D+1, (3) has order min(B+D,2D)=3-5k.
Combining the initial scalar order 2 and signed change order alpha gives
postfast order min(2,alpha,B+D,2D)=alpha. This is conditional on that
actual new fast estimate; the code does not establish it from a declared
pole or from successful polynomial arithmetic.

## 3. Fresh slow response and derivative ledger

At n and n+1 feed (4) into the SAME retained five-map coefficients:

    u_s=H_n(Kp_n Theta_P,n+Kz_n Theta_z,n),
    z_s=H_n Kg_n Theta_theta,n,
    x_s=(n+1)H_(n+1) Kpsi_(n+1) Theta_theta,n+1.              (5)

All coefficients are obtained from the existing initial budget's slow
recipe. There is no new current-scalar or slow-coefficient callback.
The optional weighted result applies its original whole-potential-hull
floor once. The response has conditional orders alpha,alpha+1. It is an
increment, not a replacement of old means. In particular computing (5)
does not yet prove the postslow reset or final row difference estimates.

For slow outputs through m, reserve Q38 pressure/scalars through m+1,
new fast fields through m+2, Q36 signed target through m+4 and its profile
and primary covariance through m+17. This finite low-order interface is
not a replacement for the full compact-tail request graph or an expanded
canonical primitive ceiling. Higher requested outputs require rebuilding
the finite ledger, not bypassing its guards.

## 4. Precise remaining estimate

The new fast input must be the actual compact/time inverse applied to the
CURRENT post-signed nonzero-Haar rows. It must satisfy, on the same bands
and whole hull, a bound of the form

    ||Fnew_k||_through j <= epsilon^(D+1_{k=x})
                            zeta delta^-v A_k,j(S),          (6)

including the axial A1 velocity. The callback supplies only the positive
norm recipe and its derivative/pole declaration; those do not verify (6).
The historical fast_envelope requires current rows at j+4 and j+5 and
new-domain KT1/KA1/time-floor estimates. Q38's radial maps do not supply
the nonzero-Haar current-row source or the oscillatory cutoff estimate.
The full A0/A1 high-order graph from core_real_refresh_obligation_join
therefore remains in force, as does its zero strict fast-difference margin.

Tests run Q38's real norm chain with manufactured lower inputs, compare
(3)--(5), retain the old/new pole distinction, detect removal of full old
correlations and check next-scalar usage. They are not actual-profile or
continuum source evaluations. The parent frozen closure is replayed.

Next derive the new fast fields from the current full rows, then apply the
fresh postslow scalar reset and complete row differences, retaining old
means, all pressure archives and viscosity/cutoff interactions. Actual
source closure, uniform pulse bounds, complete residual improvement and
blow-up remain unproved.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json` here.
