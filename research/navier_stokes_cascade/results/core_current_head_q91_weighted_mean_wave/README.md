# Q91: weighted mean-wave source from the corrected unweighted means

Parent: `aee602e9eeb10cbfb864af8eae10286cceddcd67` (Q90).
The task is Q90 (5), the actual first-real full-domain source bound at
B=1-kappa. We isolate and compose BOTH mean-wave orientations. For this
part it is unnecessary to prove a new full-edge bound for the final means:
the wave supplies the whole required half-edge/Gaussian factor. This does
not remove the weighted wave or linear-source hypotheses. No field,
pressure, endpoint, width, domain or frozen evidence is changed.

## 1. Exact cylindrical expression with the complete wave

For angular-independent m=(beta,v,gamma) and complete Fourier amplitude
w=a exp(i omega Phi), omega=kmu with fixed nonzero integer mu, exactly

    exp(-i omega Phi)[N(m,w)+N(w,m)]
      = beta Dr a + gamma Dz a + ar Dr m + az Dz m
        + i omega (m.n) a
        + R^-1(-2 v atheta, v ar + beta atheta, 0).       (1)

n=(Dr Phi,ptheta/R,Dz Phi). The last vector contains the connections
from BOTH orientations. Complete curl amplitudes need not be transverse;
no n.a=0 substitution is made. In particular incompressibility does not
cancel mean advection of the phase. Cutoff curl terms stay in a. The mean
in (1) is the FINAL corrected first-real mean, not the incoming projected
row or the post-signed/pre-new-fast mean. Its axial component retains
gd-A1 gd and its scalars retain full shifted-pressure memory.

## 2. One weighted factor suffices

Work on Q26's fixed-q domain and each whole enlarged lifted rectangle,
using the same R coefficient chart and original Gaussian P. At a finite
source order n, take all factorial through-order bounds through N=n+1.
For the current final mean, suppose its component bounds, after their
separate original unweightings, are

    ||beta||_N <= epsilon^(B+1) Mx(S),
    ||v||_N, ||gamma||_N <= epsilon^B Mu(S), Mz(S).       (2)

The meaning of the last display is one component per corresponding
coefficient Mu,Mz. Put M=Mx+Mu+Mz and drop the beneficial radial epsilon.
These are unweighted norms over the ENTIRE hull and path; bounds only on
a physical evaluation or a set of separate response bumps do not suffice.
Q86 supplies their conditional coefficients from the corrected Q47 factory:

    (Mx,Mu,Mz)=243[old_full(N)+fast_unweighted(N)
                                      +slow(N,False)]. (3)

The edge costs of old and new means are already paid separately in (3).
The filled reserved potential hull and pressure-dependent slow update are
retained. Nothing multiplies those unweighted means by zeta again.

Let the COMPLETE, actually cut first-real wave have a component-max
factorial coefficient bound U(S) through N with weight

    epsilon^(1/2) sqrt(zeta) delta^(-p) P.               (4)

Here a norm for all raw coefficient derivatives divided by the common
weight is meant, as in Q20; it is not a norm after deleting derivatives
of P. The finite p must cover this N. Q27's two regular/cut norm slots
are distinct: (4) uses its cut slot. Any conversion to the current R chart
and exact target/root state must be supplied before using U. The wave
coefficient is not manufactured by evaluating the cut field at one point.

Let g,Nnormal,r bound, through N on that same domain, respectively the
normalized graph multiplier of Dr, the l1 normal, and R^-1. Thus Dr
costs at most epsilon^-kappa g times one ordinary derivative (g can be
padded to include its identity part); Dz=epsilon partial_Z. All multiplier
jet costs are included. With factorial through norms, one ordinary
coefficient derivative costs at most (n+1) times the through-N norm.
Multiplication is submultiplicative by the multiindex Leibniz formula,
including multiplication of an UNWEIGHTED function and weighted raw jets.
The wave's component-max norm converts to a vector l1 norm with factor 3.

For 0<epsilon<=1, k=ceil(epsilon^-1/2)<=2 epsilon^-1/2 and
0<kappa<1/4, (1)--(4) consequently imply

    ||N(m,w)+N(w,m)||_(n,weighted)
      <= epsilon^B Cmw_n(S),
    Cmw_n=3[2(n+1)(g+1)+2|mu| Nnormal+2r] M U.          (5)

The weighted norm here removes sqrt(zeta) delta^-p P, and no additional
mean edge pole is introduced. The two ordinary-advection groups have
orders at least B+1/2-kappa; the phase group has order B; connections have
order B+1/2. Dropping only positive epsilon reserves gives (5). The true
radial companion improves some terms but is not needed for this bound.
For a finite Fourier sum use its l1 coefficient mass and the largest
original |mu| and multiplier bounds. Different labels can be separated
only by the retained CLOSED enlarged support theorem.

This is a direct proof for the specified term, not promotion of an arbitrary
unweighted source to a weighted one. For example sqrt(zeta) P sin(y) has
the required value weight, while sin(y) need not, even though both have
finite unweighted norms. Conversely an unweighted mean bounded by one
may multiply the former without spoiling its value weight. Higher raw
jets require the through-N versions of these hypotheses, not just values.

## 3. Executable composition and finite reserve

MeanWave requires the corrected Q47 Postslow state through Q86 FinalMeans;
it rejects a different wave-state object. It derives M from (3) and g,r
from the original background archive. The remaining U and Nnormal slots
are explicit LOWER inequalities on that identical state/domain, not
physical certificates supplied by the class constructor. No independent
final-mean norm callback remains. The program uses positive EP algebra
and preserves its exact coefficients across historical module instances.

For source order n it requests final mean and wave order n+1. Q86's
factory reserve is therefore n+1; its conservative incoming-row, covariance,
signed-target and profile requests become n+11,n+12,n+15,n+28. This is
before expanding the complete initial-mean/target/background dependency
DAG. It neither raises nor shrinks the original accepted q-domain. Bounds
at a single derivative order cannot be extrapolated to arbitrary physical
flatness or to every later sourced state.

## 4. Remaining full-source inequality

Write f_R=f_lin+f_ww+f_mw, retaining the exact regular split of Q89/Q90.
Conditional on the actual common-domain (2)--(4) and fixed multipliers,
(5) removes the mean-wave contribution from the new independent input list:

    |D^I f_R| <= n! epsilon^B sqrt(zeta) delta^-p P
                     [Clin_n+Cww_n+Cmw_n], |I|=n.       (6)

The other two coefficients require their OWN actual weighted operator/wave
estimates with a common pole enlarged to cover all three terms. In particular
unweighted Q88 radial-residual output cannot stand in for (4) or Clin.
Pressure's homogeneous moving-normal term and every regular spatial cutoff
and viscous term remain in Clin; neither temporal commutators nor compact
mean tails enter the forwarded source. Full Pi1/Pi2 memory remains in the
means of (3). No conclusion about (6) follows from named rates alone.

Next derive the common original-root weighted regular/cut wave and complete
linear-source coefficient envelopes on this current factory's R chart,
including their background normal and pressure derivative requests. Join
the wave-wave estimate and (5), then apply Q20--Q21 and Q90 to the full
source. The final actual primitive inequalities, uniform actual pulses,
complete residual improvement, infinite realization and blow-up remain open.

Exact tests independently expand (1), verify its missing-connection and
mean-phase controls, check the norm polynomial and derivative shift, and
exercise the typed corrected-state composition with manufactured leaves.
They do not evaluate the actual physical coefficient constants. Run
`python test_audit.py -v` and `python audit.py --verify-record evidence.json
--output FRESH.json`; Q90's frozen evidence graph is retained.
