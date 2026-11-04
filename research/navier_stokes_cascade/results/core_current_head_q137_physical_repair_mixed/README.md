# Q137: ordinary physical repair derivatives with a closed Gevrey-3 bound

Parent: Q136, 31ef905aa7183a9b87ba34ccda23e4c3eebdf0fb.
A fresh fetch agreed with this head; both parent workflows were in progress
on 4 November 2026. Current Q136 notes/evidence, Q121, Q132, Q63 and Q70
were read first. Q136 completes a CONDITIONAL pure-angular leading F0
assembly; it does not control arbitrary ordinary-X derivatives. This
checkpoint resolves the separable leading exact-repair contribution to
that next mixed-derivative obligation. It is not the full modulated profile.
The C^-20 activation baseline remains local leading core/first collar.
Original roots, N, gauges, bumps, physical scales and all historical frozen
records remain unchanged.

## 1. Retain the full Stirling conversion without an extra factorial

For y=log(X/rho), the exact identity is

    partial_X^r = X^-r sum_(k=0)^r s(r,k) partial_y^k.       (1)

The coefficients are SIGNED first-kind Stirling numbers. Replacing them
by second-kind numbers or retaining only the top derivative is wrong.
Q63/Q70 correctly retain this identity. Their conservative estimate using
sum |s(r,k)|=r! remains valid; it is not a sharp factorial classification.
Here is a refinement for the SAME Q63 bump sequence

    B_k=A D^k((k+1)!)^3, A=162*32600, D=32600.

For r>=1 and 1<=k<=r, the cycle decomposition of permutations gives

    |s(r,k)|/r! = (1/k!) sum_(n1+...+nk=r; ni>=1) 1/(n1...nk)
                <= binom(r-1,k-1)/k! <= binom(r,k)/k!.

Consequently, using (k+1)!^3/k!=(k+1)^3(k!)²,

    sum_k |s(r,k)| D^k((k+1)!)^3
      <= (r+1)^3 (r!)^3 sum_k binom(r,k)D^k
      <= [8(1+D)]^r (r!)^3.                             (2)

The last step uses r+1<=2^r. At r=0, s(0,0)=1 and (2) is exactly one.
This is an all-order analytical inequality, not a fit to finite tables.
No cancellation between signed terms was needed for this upper bound.
It improves the coarse radial factorial cost while leaving the old valid
frozen bounds in place.

## 2. The SAME leading five-row repair, not moving modulation

Q132 supplies, for P=2e_star max(cstar,gamma Z) and K=4+Kroot,

    |partial_y^r partial_eta^k Delta(E,U)|
           <= P K^k(k!)^3 B_r.                          (3)

The coefficients depend only on eta. The original disjoint width-.1
unit-mass bumps depend only on y; their supports and normalizations are
unchanged. The root is the same accepted SMALL C1 root of the normalized
leading quadratic five-row map, with all five baseline deficits and all
five incoming modulation memories. In particular the pressure-square row,
Pi_eta dependence and inverse-lambda normalization are not removed.
This argument does not apply to the earlier eta-dependent quadratic map
or to positive-order linear repairs without their own binding estimates.

Let xmin>0 be a valid lower bound on the ORIGINAL repair supports. Define
G=8(1+D)=260808. Combining (1)--(3) gives, for every r,k>=0,

    |partial_X^r partial_eta^k Delta(E,U)|
       <= A P (G/xmin)^r K^k (r! k!)^3.                 (4)

The flat zero extensions make this a bound on the repair contribution on
the whole X interval: outside its supports all its local derivatives
vanish. This does not say its cumulative moments vanish there. Actual
support membership and accepted-root/positive-field conditions remain
premises, not conclusions of entering numerical constants in a helper.
The inverse radius is retained; (4) is not uniform in physical scales.

For Delta F=Delta E/sqrt(2X), the multiplier has raw derivative bounds

    |partial_X^j (2X)^-1/2|
       <= max(1,xmin^-1) xmin^-j j!,  X>=xmin,

because (1/2)_j<=j! and (2xmin)^-1/2<=max(1,xmin^-1).
Padding j! to (j!)^3 and applying Leibniz gives

    |partial_X^r partial_eta^k Delta F|
       <= A P max(1,xmin^-1) [(G+1)/xmin]^r K^k(r!k!)^3. (5)

No division by sqrt(X) is used at the axis; this contribution is identically
zero on a neighborhood of it. No radial analyticity of a flat bump is
claimed. Only the physical derivative majorant is refined.

## 3. Raw total orders and the axis primitive

For a component with prefactor Pkind and radial base Rkind in (4),(5),
put B=max(1,Rkind,K). Since r!k!<=m! when r+k<=m,

    max_(r+k<=m) |partial_X^r partial_eta^k component|
          <= Pkind B^m(m!)^3.                           (6)

This is a RAW mixed maximum, not a factorial angular norm; no division by
k! is silently made. Its binary budget is
ceilLog(max(1,Pkind))+m ceilLog(B)+3m ceilLog(m), last term zero at m=0.
The separate (r,k) budget retains both factorial weights and both bases.
Both work at order 100005 without expanding factorials.

Q121's leading tuple reconstruction is linear in F,U and the ENTIRE
axis-to-X primitive M=int_0^X U. Thus (6) at m for Delta F and m+1 for
Delta U can be inserted into its SAME n=0 recipe to bound the repair's
contribution to (E,U,W,J), including its primitive. The helper does this
at moderate orders only. It does NOT supply the full background tuple.
Even where Delta U is zero after the last bump, its primitive can be
nonzero. Neither a velocity plateau nor a local derivative bound resets
pressure or either reserved mean. Q132's full target/exit identities remain
required for the completed field; a norm alone proves no pressure support.

## 4. Precise remaining obstruction

No estimate of arbitrary radial derivatives follows from Q136's angular
norm. For a manufactured family f_M(X,eta)=a sin(MX), independent of eta,
all its angular factorial norms equal a on a fixed interval containing a
maximum, while sup|partial_X f_M|=a M is unbounded. This is an obstruction
to that inference, NOT a counterexample to the selected Navier--Stokes
profile or to the existence of its radial estimates.

For the actual compact modulation, E_N=E exp(Acal/N), even its first
ordinary radial derivative contains

    (E_N)_X=exp(Acal/N)[E_X+E Acal_X/N].                 (7)

The real value gate |Acal|<=N only bounds the exponential. It gives no
bound on Acal_X/N. Q131's pure-angular estimate and (4) therefore do not
supply the whole mixed tuple. The next obligation is to bind actual mixed
X,eta derivatives of the moving-phase primitives and PREMOD fields,
retain the ORIGINAL N and inverse widths, then assemble all radial regions
on Q121's whole hull (U through m+1 on [0,Xhi], F through m on the annulus).
The low-order state/cone and same physical-scale gates still need checking.

Exact tests cover cycle/Stirling coefficients, (2), the signed physical
chain rule, E-to-F products, factorial/binary bounds, the missing radial
inference and a nonzero post-bump primitive. The parent evidence is replayed
and source hashes checked; a fresh record must match byte for byte. These
checks are not a continuum proof or an evaluation of the enormous profile.
Uniform physical pulses, complete covariance mean corrections including
axial velocity/pressure, the full viscous/cutoff/interaction residual and
blow-up remain unproved. Run test_audit.py and audit.py --verify-record
evidence.json --output a fresh path. Historical evidence stays frozen.
