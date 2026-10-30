# Q37: Q36 complete real covariance and pressure memory

Parent: `fe8ff50754441ef89ab24a09db18ebca7a2d4876` (Q36).
This connects the freshly recomputed Q36 target to the existing first-real
signed covariance and pressure identities. It retains their conditional
scope. No historical frozen evidence, profile, root, width or band is changed.
On 30 October 2026 parent focused CI 41127207604 passed; broad CI
41127207607 was still running at inspection.

## 1. All six entries, including the signed square

For the first phase-locked REAL update write the complete old wave as
A cos(theta)+B sin(theta), and its signed increment as
L cos(theta)+S sin(theta). The coefficients here include their physical
powers and full cutoff curls. Normalized angular averaging gives exactly

    Delta W = Lambda + E,
    Lambda = (A tensor L + L tensor A)/2,
    E = (B tensor S + S tensor B + L tensor L + S tensor S)/2.   (1)

This is an identity for every symmetric tensor entry, not only the two
selected fluxes. The selected Lambda entries match the chosen target under
the original phase/root/partition identities. The other four need not vanish.
The signed square L tensor L cannot be dropped. Later complex waves have
additional cosine-quadrature crosses, so (1) is not a later-prefix rule.

Use the inherited common-chart conversion on each primary and Q36 signed
coefficient jet f_j:

    fbar_n=16*188^n*n! sum_{j<=n} f_j((25/16) Sscale).           (2)

These bound raw maxima in the common native U chart under the original
atlas/label normalization contracts. The physical target transfer is the
same one justified in core_symbolic_real_covariance_join: with r=Q*/Qi,
r^(h*alpha-2*Aexp)=r^(-1-3h*kappa)<=1. It is not a new assumption equating
unscaled targets. The same physical zeta is retained.

Let a,b,l,s denote (2). Applying Q34's established raw-U to through-order-R
transfer T_n ONCE gives sufficient bounds for all entries:

    Lambda_n = T_n 2^(n+1) a_n l_n,
    E_n      = T_n 2^n (2 b_n s_n+l_n^2+s_n^2).                (3)

The complete norm is their sum. T_n already bounds the sum of raw and
normalized R derivatives; no second derivative-count factor is needed.
The unchanged common wave and covariance poles are

    q_n=max(6n+9,6n+d+6),  covariance pole=2q_n, d=2p.          (4)

Q36 supplies p=12*(m_target+13)+18. Derivatives of flat factors and complete
curls are already in those inherited wave allowances; they are not removed
by treating the principal target as constant.

With kappa=1/100000 and alpha=2-3kappa, the three terms in E have orders
alpha+1-2kappa, 2alpha-1 and 2alpha-2kappa. Thus the sufficient error order
is nu=3-6kappa, while Lambda and the complete change have only order alpha.
These are conditional norm estimates, not a complete residual improvement.

## 2. Fresh pressure changes, old scalars retained

Hold the existing mean velocities fixed during this wave substep. For the
Haar-averaged tensor change D, the original identities read

    dg=-D1 Drr-epsilon partial_Z Drz+Dtt/R,
    dP=int dg, dp=I0(dg-rho dP),
    dJtheta=int R^2 Dtheta,z,
    dJz=int R Dzz-(1/2)int R^2 dg.                            (5)

The helper uses the existing complete pressure_target_budget with (3) at
n,n+1,n+2. No axial covariance, pressure compensator, radial derivative or
cutoff-curl component is deleted. For every scalar it returns

    current = Q36 initial postslow scalar + signed-wave change. (6)

This is addition of positive bounds for the exact addition of fields, not
a replacement of the old pressure by its change. It computes both

    epsilon partial_Z(dJz+c_rho dP),
    epsilon partial_Z(Jz_current+c_rho P_current),              (7)

with the same Q36 fixed profile, including c_rho,Z P. Background Pi1/Pi2
and previous mean/wave pressures remain in the initial state premises.

The density and c_rho norms are obtained directly from Q36's profile.
For whole-hull radial weights -1,0,1,2 a sufficient moment constant at
new pole v is

    Kmom(v)=(Rmax-Rmin)(Rmin^-1+1+Rmax+Rmax^2) Eedge(v,width). (8)

Indeed zeta delta^-v<=Eedge on the entire original radial norm hull; each
radial weight is bounded by the corresponding displayed term. Normalized
Haar averaging has norm one. This avoids applying the reserved-support
weight floor outside its actual moving support. No independent density,
moving moment or radial moment callback is accepted here.

## 3. Exact remaining new-domain inequality and finite ledger

Pressure-target output n requires covariance through n+2 and signed target
through n+3. To request pressure outputs through m, construct Q36 through
m+3, hence reserve its profile and primary covariance through m+16. This
is a safe composition ceiling, not the expanded canonical primitive ceiling.

The new pressure maps MUST cover order n+1 and pole

    v=2*max(6*(n+2)+9,6*(n+2)+d+6).                           (9)

In particular v>d. The old initial pressure contract at d is insufficient.
The callback supplies only Kp, primitive and bump, with declared order/pole.
Its declarations are checked, but the actual inequalities remain necessary:

    ||I0(g-rho int g)||_(n+1,v) <= Kp ||g||_(n+1,v),

and the corresponding oriented target primitive and normalized projector
bounds on the SAME domain and supports. The precise norm/edge allowances
are those of the historical pressure_target_budget contract. A finite
number entered in a callback is not a proof of these inequalities. The
existing symbolic two-edge kernels offer a route at this new finite pole;
this checkpoint does not silently reuse an older-domain certificate.

The complete Haar pressure/scalar change has order alpha; the pointwise
pressure may have only alpha-kappa. The updated projected target has the
inherited order nu only under all exact selected-flux and pressure moment
identities. No full pointwise residual rate is promoted.

## 4. Checks and next obligation

Tests independently verify all six tensor identities, detect an omitted
signed square and a later-complex cross, check the limiting exponent and
strictly enlarged pressure pole, and run the Q36-to-pressure chain with
manufactured lower bounds. They check old scalar addition, distinct change
and total moving-pressure moments, chart transfer, full curl and finite
ledger rejection. These are algebra/interface controls, not continuum
suprema of the actual selected field. Q36's frozen closure is replayed.

Next discharge the new-domain map inputs and compose fresh fast/slow means
on this CURRENT post-signed state, retaining old means and both pressure
archives, before binding the first complete unaveraged source. Actual
root/map/support premises, uniform pulse estimates, full viscosity and
cutoff interactions, complete residual improvement and blow-up remain open.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json` here.
