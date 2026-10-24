# Sharp initial scalar cancellation joined to symbolic primary and fast bounds

Parent: `d5ef67446a724e53913590a682f868f545e73cfb`.
On 24 October 2026 its fresh-mean CI step passed; the full run was
still in progress without failed steps. The parent correctly rejected
promotion of its COARSE scalar order 1-kappa to order1. The stronger
INITIAL identities already exist in `core_initial_scalar_seed_budget`.
This checkpoint joins those identities to the explicit charted primary,
actual fast-state norm recipe and symbolic coefficient ring. It does not
claim the old S2 cancellation as a new theorem. Historical evidence,
profiles, widths, supports, physical bands and gauges are unchanged.

## 1. Estimate after the exact Haar and radial identities

This is ONLY the first real-primary prefix: mean perturbations are zero
before the fast correction, and all three actual fast velocities have
zero auxiliary Haar mean. Keep the complete angular covariance W,
including its axial/radial entries and real curl squares. Integration
by parts of the exact averaged radial source gives

    P0=int (Wtt-Wrr)/R - epsilon partial_Z int Wrz,
    Jtheta0=int R^2 Wtheta,z,
    Jz0=int R[Wzz-(Wrr+Wtt)/2]
                            +(epsilon/2)partial_Z int R^2 Wrz. (1)

All W entries in these integrals are Haar averages. Haar kills graph
torus derivatives exactly; it does not kill correlations. Boundary
terms vanish by the ORIGINAL smooth compact supports on the whole
fixed positive hull. These operations precede taking absolute values.
Consequently no epsilon^-kappa graph cost occurs in (1), unlike the
parent's conservative full-row estimate.

Let [a,d] contain the ENTIRE primary and actual fast potential/cutoff
supports, including A1. This is not assumed to be the reserved I4.
Put M_j=(d-a)max(a^j,d^j). All derivatives retain those supports.
The explicit R-chart covariance coefficient W_n has full weight zeta
and a common pole p, with epsilon order1. The retained edge absorption
Eedge(p,w) bounds sqrt(zeta)delta^-p, hence also zeta delta^-p.
Set Omega_n=Eedge W_n (no torus-count multiplier). Then (1) gives

    IP_n=2M_-1 Omega_n+(n+1)M_0 Omega_(n+1),
    Itheta_n=M_2 Omega_n,
    Iz_n=2M_1 Omega_n+((n+1)/2)M_2 Omega_(n+1).        (2)

The derived covariance already dominates normalized through-order norms;
no old two-variable raw-count conversion is silently reused for five
variables. The factor n+1 safely pays the one slow derivative. The
whole-hull norm and support premises remain required.

## 2. Fresh fast scalar changes retain every quadratic term

Write the ACTUAL zero-Haar fast velocities as beta_f,v_f,gamma_f,
where gamma_f=gd-A1 gd, not gd. Exact averaging cancels their linear
time, base and viscous terms, including radial velocity, but leaves

    Delta gr=-D1 bar(beta_f^2)
             -epsilon partial_Z bar(beta_f gamma_f)+bar(v_f^2)/R.

Thus the three scalar changes are

    Delta P=int [bar(v_f^2)-bar(beta_f^2)]/R
                          -epsilon partial_Z int bar(beta_f gamma_f),
    Delta Jtheta=int R^2 bar(gamma_f v_f),
    Delta Jz=int R[bar(gamma_f^2)-(bar(beta_f^2)+bar(v_f^2))/2]
                          +(epsilon/2)partial_Z int R^2 bar(beta_f gamma_f).

These are exactly the old zero-old-Haar fast identities. Use the newly
constructed fast fields through n+1, remove their full edge weight with
Eedge, and call the historical `fast_budget(n,M,0,0,T,X)` plus
M_2 T^2 for the angular change. The tangential/axial rates are mu=1-kappa,
the radial rate mu+1. Every retained quadratic term has order at least
2mu>=1. Adding (2) therefore proves FRESH postfast P,Jtheta,Jz have
epsilon order1 with explicit positive EP coefficients. The scalar bound
does not require the false inequality epsilon^mu<=C epsilon.

Pressure's zero-Haar projection may simplify the fast input, but its
WHOLE moving scalar P remains here. No compensator or moving product
is deleted from subsequent pressure and axial residual identities.

## 3. Join the original slow solve and its C=1 remainder

Feed these fresh scalar coefficients into the SAME fixed slow inverse
at orders n and n+1. Its original distinct physical normalization
powers, reciprocal a floor, fixed lambda loss and total mixed-profile
jets remain lower inputs. The output field bounds now have slow
tangential/axial order1 and radial order2. For the sharper weighted
Haar interface, multiply these differentiated-field coefficients by
the parent's SYMBOLIC exp(E) floor on the actual moving potential
hull, with pole0. For scalar integrals use the unweighted response
coefficients and the WHOLE reserved union hull separately.

On that original reserved hull require the existing exact base identity
G=0, the prescribed V profile and all five moment-response identities.
They are hypotheses of this join, not consequences of norm bounds.
The actual base radial velocity b=epsilon b_N is retained. Old fast
times slow products vanish ONLY after Haar averaging because the old
fast means are zero. The historical C=1 `slow_budget` with r=a=0
therefore bounds the freshly recomputed nonlinear defects at order2.
This is the actual postslow residual of the five-map, not a reset of
old mean velocities or an assertion that the full PDE residual is zero.

The complete axial residual moment is bounded at order3 by

    (n+1)[Theta_z,n+1+C_rho,n+1 Theta_P,n+1],          (3)

where C_rho bounds the WHOLE moving c_rho in the normalized through
norm. This retains c_rho,Z P. The scalar input order for (3) is n+1;
the code never substitutes Theta_n. The geometry comparison and these
initial cancellations were proved historically; the new result is their
explicit substitution with the recent primary/fast symbolic recipes.

With this conservative route, scalar output n uses primary covariance
through n+12 and the whole axial moment through n+13. Fast output j
uses covariance through j+8; postfast scalar n uses fast fields n+1;
slow fields add one scalar jet for the radial response; the C=1 scalar
remainder uses slow fields through n+2. These are safe finite counts,
not sharper than the historical separately optimized scalar route.
Further root/background costs inside covariance remain unchanged.

## 4. Scope and next obligation

The exact historical pressure and fast/slow identity controls are replayed.
New controls compare the symbolic quadratic and C=1 remainder recipes,
check the common scalar and moving axial-memory orders, and trace the
explicit primary covariance route. Pipeline traces use a marked
manufactured covariance subclass, not the selected giant physical field.
Passing tests is not a continuum proof.

The parent's coarse bound remains valid; it was insufficient by itself.
The initial sharper order is now conditionally supplied using additional
EXACT zero-Haar, compact-boundary, reserved-profile and five-map identities.
It must not be extended to a later nonzero old-Haar/complex state. Full
background/phase/value gates, actual whole-hull map/profile bounds and
support identities remain premises, including the declared fast hull.
No new cutoff, first band, tilt or radius is selected to make them hold.

Next substitute these sharp initial slow/Haar coefficients and fresh
order2 scalar defects into the existing weighted initial-row and
original signed-target/forced-source budgets, retaining complete
pressure memory and later complex cross terms. The complete all-scale
viscous, cutoff and nonlinear residual, infinite realization and blow-up
remain unproved; no unconditional pulse or residual-improvement claim
is made.
