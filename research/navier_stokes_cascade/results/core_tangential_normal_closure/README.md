# Close the extra low-order inputs to the normal-floor recipe

Parent: `c86d5ad47289d9484f060801ea4727247d61c51d`.
This checkpoint removes three independent inputs from the sufficient
normal-floor test: the coefficient-one tuple norm, Borel cutoff K2, and
leading second jets. No new analytical premise replaces them. The SAME
complete normalized Euler diagonal and the actual leading first jets,
chart coverage, value gaps and original selector data remain premises.
This closes an input dependency, **not** the actual unchanged-band floor.
No unconditional pulse, complete residual improvement or blow-up follows.

## 1. Use the retained diagonal at its natural low-order rate

Write E=q partial_q at fixed X,eta. Retain the COMPLETE normalized tuple
from `core_background_finite_head_budget`, equation (1), including

    q^(2nh)[chi W_n^r-(E chi)J_n^r]

in its radial component. Coefficients are the same globally finalized
ones with all viscosity, products, pressure memory and moment repairs.
For every n>=1, the retained diagonal gives, through total Euler/X/eta
order m<=1 and on the entire cutoff support,

    |jet increment_n| <= C_(n,m) q^(2nh)
                       <= 2^-n q^(nh).                         (1)

Outside that support the increment and its jets vanish. Hence for ALL
0<q<=qstar<=1, putting t=q^h,

    sum_(n>=1) 2^-n t^n = t/(2-t) <= t,
    |background correction tuple|_(m<=1) <= q^h.                (2)

Uniform majorants justify the differentiated sum at these orders. This is
an analytic infinite-series argument, not a finite-prefix certificate.
It is conditional on the SAME diagonal for every n, including n=1.
The coefficient and cutoff norms are still used to DEFINE/justify that
diagonal. Their identities and bounds are not being inferred from (2).

The previous q^(2h) finite-head estimate is correct and remains useful at
higher precision/orders. Its coefficient-one head is unnecessary at the
weaker q^h rate already allowed by the phase normal. Merely dividing (1)
by q^(2h) at n=1 would still be invalid. At m>=2 the n=1 term has no
diagonal jet guarantee; (2) does not extend to those orders.

On the original physical chart 1/4<=s=q/Q<=4, epsilon=Q^h, the original
first-order composition constants CF1,CG1 yield

    Jfull1 = (CF1+CG1)(Lead1+1),
    |F-Flead,G-Glead|_1 <= 2(CF1+CG1) epsilon.                  (3)

Here Lead1 is a bound for the complete LEADING native tuple through one,
and 4^h<2 for the construction's 0<h<1/100. Composition differentiates
the actual increment, never a pointwise envelope. The native geometry and
one-sided boundary restrictions are exactly those of the parent. This
does not certify frozen representatives outside the physical covered set.

## 2. Bound the tangential normal with its complete axial transport

The actual phase normal is

    n=(x0-v H_R, p/R, pz-epsilon v H_Z), H=pF+pzG.               (4)

For an inverse-normal floor it suffices to bound its last TWO components
against Bs Kunit, with |Kunit|=1. Radial closeness to its tilted reference
is unnecessary for that particular implication. This is not a component
floor on p/R: p/R can be arbitrarily small. The complete axial component
and BOTH unrounded tangential tilt errors are retained.

Use the unchanged geometry Llow S<=L<=Lhigh S, S>=1, and the same bounds
U>=u, lambda_max>=lambda0, |g0|>=gmin>0, Rmin<=R<=Rmax. Put

    B=max(1,lambda_max), qz=B[1+U/(Llow gmin)],
    p0=Rmax qz, pplus=p0+1,
    ct=2 B U/(Llow gmin), ca=1/Rmin,
    cb=p0 Dbox/Rmin^2,
    cz=Lhigh(pplus+qz) Jfull1,
    Ctan=ct+ca+cb+cz.                                          (5)

The original nearest-NONZERO carrier rule costs 1/k, not 1/(2k).
The exact frozen carrier relation in the parent's equation (9), box
diameter Dbox S^-3, (3), and (4) give the Euclidean tangential error
via the component-sum upper bound

    |n_tan-Bs Kunit| <= ct/S+ca/k+cb/S^3+cz S epsilon.           (6)

In particular cz bounds the FULL H_Z of the summed background; no leading
substitution, local U=0 pressure reset or omission of axial velocity is
made. The leading second jets used to cancel the RADIAL normal are not
needed in (6). Those jets and radial errors remain relevant to other
phase/eikonal/cone estimates; this result only replaces the floor recipe.

Retain the original whole-path slow gate

    S^2(epsilon+epsilon^2+k^-1)<=1.                             (7)

Each of the four contributions in (6) is bounded by its OWN coefficient
divided by S. Thus their sum is <=Ctan/S, with no extra factor 2 or 3.
Carrier calibration gives Bs^2>=beta=lambda_min/[4(1+U^2)^2]. If

    4 Ctan^2 <= beta S^2,                                      (8)

the triangle inequality gives |n_tan|>=sqrt(beta)/2, hence
|n|^-2<=4/beta. This applies to the full moving normal used in the sourced
inverse. It does not remove that inverse's normal-source pressure term,
propagator hypotheses, supports, or complete covariance corrections.

## 3. What is closed, and what must still be checked

The previous recipe used 36 Cphi^2<=beta S^2, where Cphi included radial
background comparison and leading second jets. That remains sufficient,
but failure of its coefficient-one/K2-dependent margin is NOT an
obstruction to (8). For its SAME input geometry and leading norms, the
old finite head Bcorr1>=3/2 gives old Jfull1>=(CF1+CG1)(Lead1+3/2).
Consequently old Cphi strictly exceeds Ctan; no choice of K2, first repair
constant or radial second jet enters (5). A toy test where the old recipe
fails and (8) passes is a logical separation of sufficient tests, not a
fluid counterexample or certification of an actual band.

The unchanged selector supplies S=ell^2 and ell^2>=Lreq/(2r0).
A sufficient compatibility check on that SAME length input is now

    beta Lreq^2 >= 16 r0^2 Ctan^2.                             (9)

It replaces 144 r0^2 Cphi^2 for THIS floor only. Alternatively check (8)
directly at the actual existing minimum ell; its left-to-right margin
then improves for every later band. Neither operation reselects ell,
Lreq, r0, carriers, supports, qstar or the unchanged width C^-20.

Removed assumptions: independent coefficient-one, K2 and leading-second
input enclosures **for the floor**. Added assumptions: none.
Leading tuple first jets still require the radial/angular profile inputs
in the parent's reconstruction; this does not turn a first-tuple jet
into a bound using only first derivatives of each primitive profile.
Retained unresolved premises: identify/justify the same full infinite
Euler diagonal and cutoff; verify the actual chart/value neighborhood,
physical representatives, leading first norms and positive gaps; evaluate
(7)--(8) or (9) with those inputs. The complete phase/residual still needs
its other radial, pressure, viscosity, cutoff and interaction estimates.

Next obligation: discharge the original diagonal's prescription/identity
and substitute the actual leading first-jet and gap bounds into (8) at
the existing band. A numerical infinite list is not required if the SAME
prescription is proved to satisfy the diagonal by induction; a finite
fixture alone cannot establish it. Do not demand first-coefficient/K2
enclosures independently a second time for the floor once that is done.

## 4. Reproduction and scope

`python test_audit.py -v` and `python audit.py --verify-record evidence.json
--output fresh.json` check exact geometric-tail identities, rate/order
limits, the four transport terms, original cutoff dependence of the
diagonal, strict comparison with the old sufficient recipe, unchanged
selector algebra, and frozen source bindings. Finite rational controls
are implementation checks, not a continuum proof of the inherited
profiles, infinite schedule, accepted bands, covariance or PDE residual.
