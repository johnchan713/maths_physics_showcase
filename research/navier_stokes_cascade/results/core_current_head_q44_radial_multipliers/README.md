# Q44: derive the native hull and radial compact multipliers

Parent `6ac762d58d76a714b8dc579a81efedd62e9e86c5` (Q43).
On 31 October 2026 focused CI 41161995153 passed; broad CI 41161995175
was still running. Historical evidence and all physical choices are unchanged.

Q43 left four positive radial enclosures as inputs. Two suffice: on the
original fixed hull 0<r<=R<=s and for 1/3<d<5/2, U=R^d obeys

    u=min(1,r^3) <= U <= v=max(1,s^3).                    (1)

For R<=1 use R^d>=R^3 and R^d<=1; for R>=1 use 1<=R^d<=R^3.
Thus (1) is proved from the existing R hull, not guessed from four positive
numbers or chosen by altering a support. Q44 uses these SAME enclosures
for the native raw multiplier bounds and both Q43 chart transfers.
The underlying assertion that the physical support lies in [r,s] remains
an actual input premise. No narrower support is inferred.

## Raw radial multipliers at every fixed order

For a fixed real exponent b, let m=ceil(|b|) and

    H_b=max(1,r^(-m),s^m).

Splitting at R=1 gives R^b<=H_b. Since R^b=U^(b/d), exactly

    partial_U^j R^b = (b/d)_j R^b U^(-j),

where (a)_j=a(a-1)...(a-j+1), (a)_0=1. Therefore the RAW through-n bound

    B_n(b)=H_b max_(0<=j<=n) |(b/d)_j| u^(-j)            (2)

holds uniformly on the whole hull. Other coefficient derivatives vanish:
these multipliers depend only on U, with the fixed exponent d and labels
held constant. This statement does not apply to the moving cutoff or
pressure profile, whose slow derivatives must still be supplied.
For channels e=0,1,2, the original Q29 multipliers are

    f_e=R^(e-d+1)/d,  g_e=R^(-e),  k_e=d R^(d-1-e).

Their bounds are B_n(e-d+1)/d, B_n(-e), and d B_n(d-1-e).
Also dU/dx=d U/2 for the logarithmic integral, so Jmax=d v/2.
These four inputs to Q42 are now derived. The falling factorial in (2)
is essential; neither a normalized Taylor coefficient nor a value-only
bound is a raw higher-derivative bound. For rational d,r,s the recipe is
exact rational arithmetic, with no evaluation of irrational powers.
For each fixed n these constants are finite and band-independent on the
fixed hull. They are not uniform as n increases or r approaches zero.

## Bind the current provider without dropping the unresolved terms

`build` creates Q42's raw provider, then Q43's normalized R-chart wrapper,
using the same d and hull throughout. It accepts only the seven remaining
bounds chi,1-chi,chi_U,length,edge,kernel,interior, with Q42's order/pole/
cutoff coverage declarations, plus the original time-inverse/S bound.
It rejects attempts to override f,g,k or the Jacobian independently.
The axial fast map uses e=1; e=0 and e=2 are separate moment channels and
are not substituted for e=1. No monotonicity in channel or pole is used.

Q29's directed kernel and transition weight floor still require the actual
fixed cutoff plateaux, their uniform interior gap, and their raw jets.
This checkpoint does not infer these from (1), does not identify the
logarithmic length with a native length, and does not discard the complete
shifted A1 moment or any pressure contribution. Q40 still uses gd-A1 gd,
and Q41 retains the complete moving pressure product.

Exact finite controls compare (2) with differentiated fractional-power
multipliers on manufactured rational points, verify fractional exponents
on both sides of R=1, and check the Q42/Q43 composition and channel guards.
Frozen Q43 evidence and all transitive bindings are replayed. The analytical
inequalities above establish the conditional bounds; test points alone
are not a continuum proof or a certificate of the selected support.

Next supply the actual cutoff/kernel/weight data and the initial A1 tail's
high-order source estimate. Boundedness of A1 still does not supply its
M^(-p) integration-by-parts gain. Uniform actual pulse bounds, complete
source/residual improvement and blow-up remain unproved.
