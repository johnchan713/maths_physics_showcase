# Q113: compact-shape identity and fixed-order uniformity scope

Parent: `7d98b37952bbb282a44c00478e647a8a1ce3e1e9` (Q112).
On 3 November 2026 its research CI 41557585205 and field CI 41557585148
were running; Q111 field CI 41549497511 had passed. The inspected remote
head was Q112. This checkpoint applies an EXISTING Q30 qualitative fact
to Q112's current raw profiles and corrects the interpretation of the
remaining shape obligation. It does not select a cutoff or reconstruct
completed mean rows. Historical frozen evidence, physical root, accepted
cover and width C_setup^-20 are unchanged.

## 1. A fixed smooth compact shape already has finite constants

Q30 section 1, following compact_mean_corrections (4), fixes ONE original
physical profile hat_chi with

    Xa<c0<c1<Xb, hat_chi=0 near the inner side,
    hat_chi=1 near the outer side, chi=hat_chi(X), X=R^2/(2s).

The profile is C-infinity and independent of the running band parameter.
For each fixed finite N, compactness therefore gives

    sigma_j=sup_[c0,c1] |hat_chi^(j)|/j! < infinity,
    Sigma_N=sum_(j=0)^N sigma_j < infinity.                 (1)

The plateau values account for the outside of the transition interval;
in particular sigma_0>=1. The derivative orders, original shape and hull
are fixed BEFORE varying S/epsilon. Q30 already states that numerical
values are not needed for this existence conclusion. Q111/Q112's words
"actual shape controls" must not be read as demanding new numerical
values before this qualitative fixed-order profile lemma can be used.
The numerical certificate remains absent; qualitative finiteness does not.

On the SAME positive normalized hull a<=R<=b, put u=min(1,a^3),
H2=max(1,a^-2,b^2), H=512(1+Xb)/u. Under the original Q46 coordinate
inequalities, Q112 gives the following EXISTENTIAL uniform profile bounds:

    N_N(chi)<=N! H^N Sigma_N,
    A_J=N_J(f0)<=3 H2 J! (3/u)^J,
    B_n=N_n(k0 chi_U)<=
        2^n [3 H2 n! (3/u)^n] [(n+1)! H^(n+1) Sigma_(n+1)]. (2)

Every right-hand side is finite and independent of the running band for
each fixed n,J and the original fixed hull/shape/coordinate controls.
No analyticity of hat_chi or convergence of its infinite Taylor series is
needed: Sigma_N is a FINITE sum. There is no bound uniform in derivative
order, shape choice, vanishing collar width or arbitrary physical hulls.
The slow Z/T motion of X is retained in Q112's coordinate majorant.

Thus the raw kernel MULTIPLIER finiteness slot in the qualitative compact
lemma needs no new shape construction. This is a binding and scope
correction of existing results, not a new proof of Q30. It does not bound
the actual high-jet SOURCE merely because that source is smooth for each
epsilon. A family of smooth sources can have unbounded normalized jets.

## 2. Three cutoff roles do not identify one another

The repository supplies explicit formulas for other profiles. Their
derivative certificates cannot be routed to (1) without proving equality
of the profiles and paying the argument conversion:

| Role | Argument | Existing provenance |
|---|---|---|
| Compact mean gluing | X=R^2/(2s) | Q30 and compact_mean_corrections (4) |
| Borel summation | a_n q | core_borel_cutoff_certificate |
| Modulation variance | 8(v-2)/delta-1 | Q66 cutoff amplitude |

The Borel profile uses a squared-flat exponential step. Q66 selects an
integrated bump and proves its own factorial derivative bound. Neither
statement identifies the original compact hat_chi(X) with that function.
The compact proof only needs its original fixed smooth choice; replacing
it by an explicit profile for convenient numbers would change the data.

`CompactShapeScope` binds Q112.KernelProfiles to that compact-X role and
the identical current source state. Its qualitative flag is explicitly
conditional on the original fixed C-infinity choice and coordinate
hypotheses. It rejects Borel/variance role substitution and keeps numerical
shape certification false. These checks enforce routing, not regularity
proofs. The existing source archive, exact common pole, bundle lengths,
frozen width and geometry guards are delegated without relaxation.

## 3. Why plateau geometry alone cannot supply numerical derivative bounds

Here is an analytical COMPARISON class, not a replacement actual cutoff.
Fix c0<c1 and any nonnegative C-infinity bump g supported in [-1,1],
positive on (-1,1), of mass I>0. Set

    S(t)=I^-1 int_(-infinity)^t g(v) dv,
    m=(c0+c1)/2, w_L=(c1-c0)/(4L), L=1,2,...,
    chi_L(X)=S((X-m)/w_L).

Each member is smooth, monotone, between zero and one, identically zero
for X<=m-w_L and one for X>=m+w_L. It has the SAME outer plateau
witnesses c0,c1. The fundamental theorem of calculus gives

    int_(m-w_L)^(m+w_L) chi_L'(X) dX=1,
    sup |chi_L'| >=1/(2w_L)=2L/(c1-c0) -> infinity.        (3)

Consequently there is no shape-class-uniform numerical first-derivative
bound from just smoothness, range and those plateau witnesses. This is
consistent with (1): each SELECTED member has finite constants, which
depend on that member. The relevant quantifiers are

    for each fixed hat_chi and N, exists finite C_(N,hat_chi),

not a single C_N for every smooth cutoff with those plateaus. This family
does not say the original fixed shape is ill behaved, nor is it a
canonical-source or Navier--Stokes counterexample. It changes no profile
in the current construction.

## 4. Precise source inequality still needed for the radial cone

Q110's exact current radial identity remains

    gr=r0+Dr K0, pm=K0+T0(r0-rho P),
    Er=-rho P-A0 Qcirc r0, P=int Haar r0 dr,
    K0=-S+epsilon D1 beta.                              (4)

All six final means, total old-plus-sourced covariance, incoming histories,
moving rho derivatives and WHOLE pressure memory remain. Removing Dr K0
under A0 does not remove K0 from pressure or its upstream field effects.
At target C=3/2-kappa the scalar compensator is already sufficient. The
remaining estimate is Q110 (3) for the actual A0 Qcirc r0, with
J=n+2m+3, B=1-kappa and m chosen so B+m kappa>C.

After the legitimate own-pole unweighting and R-to-U conversions, one must
still prove a band-uniform native coefficient bound for

    N_J(Qcirc r0/epsilon^B),                             (5)

and its common high-pole/archive/whole-hull/mixed-trace hypotheses. Equations
(1)--(2) make the finite raw profile factors in that inequality available
qualitatively, under their inherited hypotheses. They do not prove (5).
For kappa=1/100000, n=0, the current request remains m=50001, J=100005;
the full source program still requests profile jets through 100046.
The three conservative source-side J! instructions remain. Huge but
fixed factorial constants do not themselves obstruct a fixed-order
asymptotic estimate; the missing issue is the uniform ACTUAL source bound
and legitimate proof domain, not their numerical evaluation. A low-order
archive still cannot certify the request by changing its maximum label.

If the compact profile or normalized geometry instead depended on the
band, (1) alone would no longer give a band-uniform constant. That is not
the fixed-choice premise of Q30. Uniform estimates across all physical
scales still require the original normalization and transfer losses.

## 5. Verification and next obligation

The exact comparison checks verify fixed plateau witnesses, shrinking
transition length and (3). An independent squared-flat smooth step has
midpoint derivative 8 in its unit argument; its rescaling checks derivative
growth without a numerical quadrature. A polynomial unit-mass control
has only C^2 endpoint regularity and is labelled accordingly; it is not
the C-infinity analytical witness. Typed current-provider tests retain
the large request, role rejection, same-state and archive guards, and
separate qualitative finiteness from numerical certification. Frozen
replay verifies unchanged Q112 history and source bindings. None of these
finite diagnostics proves the continuum source inequality (5).

Next prove (5) on the whole consistent current high-jet/common-pole graph,
or obtain a sharper actual A0 cancellation, then join the radial estimate
to Q108/Q109 and the inner/outer cones. Quantitative shape values are needed
for a numerical certificate, not for the fixed-order existence lemma.
Native coefficient/reserved-mean supports, uniform actual-pulse bounds
on every physical scale, complete viscous/cutoff/interaction residual
improvement, infinite realization and blow-up remain unproved.
