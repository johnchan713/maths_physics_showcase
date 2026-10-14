# All-order background and divergence-preserving asymptotic summation

Status: `conditional-all-order-background-with-flat-augmented-residual`.

This checkpoint conditionally extends the finalized first background
coefficient to every order and sums those coefficients as smooth fields for
q>0. The cutoff sum preserves incompressibility, the heat exterior, the
reserved mean-patch velocity, and the background's growing core swirl.
Its **augmented** residual, after retaining the signed annular stress,
vanishes faster than every fixed power of q with all fixed physical
derivatives. The bare momentum residual still contains that stress's
divergence. Physical waves, admissible forcing, complete-field energy and
blow-up therefore remain unverified.

This is a conditional analytic audit, not a numerically evaluated huge-N
field, independent peer review, or a proof-assistant certificate. Its
foundations remain the inherited leading profile, analytic-eta collar,
exact moments, continuum matrix inverse and heat-edge estimates. The
manuscript already proposes this background construction; the work here
derives its general coefficient sources, tests their omitted terms, and
makes the summation losses and choices explicit.

Parent: `305b9288a86877ba994514e514835a0c8e0f036d`. Source: the supplied
[Navier-Stokes manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(5.1)-(5.6), Lemmas 5.1-5.4 and (5.41)-(5.46), especially the shifted
radial axial-viscosity term and the curl-cutoff commutator. SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The proposed theorem is not used as a premise. Earlier audit folders and
GOAL.md remain unchanged; the leading frequency stays `1+floor(H^32)`.

## 1. The general order has a new radial-viscosity source

Use the prior notation
`A=1/2+h`, `D=1/2-h`, `lambda_n=2nh`, `b=-A-1/2`, `c=-A`,
`d=1-eta^2`, `L=1-2h eta^2`, `F_n=E_n/sqrt(2X)`, and
`M_n=integral_0^X U_n`. The actual viscosity is one; the numerical
diagnostics also use other positive viscosities to expose missing factors.
The physical coefficient fields are

\[
u_{\theta n}=q^{-A+\lambda_n}RF_n,\quad
u_{zn}=q^{-A+\lambda_n}U_n,\quad
r u_{rn}=q^{\lambda_n}V_n,\quad
p_n=q^{-2A+\lambda_n}\Pi_n.
\]

The index n labels correction order and is unrelated to the fixed radial
oscillation frequency N. Define, for every lower order k,

\[
\overline\Omega_k=T_{\lambda_k}V_k-2\nu XV_{k,XX}
+\sum_{i+j=k}\{V_i(V_{j,X}-V_j/(2X))+U_iZ_{\lambda_j}V_j\},
\]
\[
\Omega_k=\overline\Omega_k
-\nu Z_{\lambda_{k-1}-D}Z_{\lambda_{k-1}}V_{k-1},
\]

where a negative index means zero. The pressure equation at order n is

\[
(\Pi_n)_X=2F_0F_n+P_n,\qquad
P_n=\sum_{i=1}^{n-1}F_iF_{n-i}-\Omega_{n-1}/(2X).
\]

Thus **order two already contains axial viscosity of V0** in its pressure
source. It must not be postponed to the third order. Omega_(n-1) includes
radial time differentiation, radial viscosity, all lower-order products,
and the shifted axial viscosity of order n-2.

All apparent axis divisions are regular. Put `v_j=V_j/X`. Directly,

\[
Z_k(Xv)=XZ_{k-1}v,\quad
X^{-1}Z_{k-D}Z_k(Xv)=Z_{k-1-D}Z_{k-1}v.
\]

For example the shifted term in Omega_(n-1)/X is exactly
`-nu Z_(lambda_(n-2)-1-D) Z_(lambda_(n-2)-1) v_(n-2)`.
The other terms use only smooth v_j, its first two derivatives, and U_j.
`hierarchy.py` evaluates the sources in this regular form, including X=0.

## 2. The same sparse linear system works at every order

With all orders below n already FINALIZED globally, define

\[
W_n=(F_n,U_n,M_n/X-U_n,\Pi_n,\partial_\xi F_n,\partial_\xi U_n),
\qquad \xi=\sqrt X.
\]

The system is

\[
\partial_\xi W_n+\xi^{-1}\operatorname{diag}(0,0,2,0,3,1)W_n
=A_{0,n}W_n+A_1\partial_\eta W_n+f_n.
\]

The entries of A0_n are the first-order entries with 2h replaced by
lambda_n in the CURRENT-order terms. The derivative matrix A1 is
independent of n: it maps the first four coordinates into the last two,
annihilates the last two, and has that same structure after every eta
derivative. Pressure and transport powers are shifted separately; a
blanket replacement in lower-order sources would be incorrect.

For an explicit source formula, put

\[
S_{\theta n}=\sum_{i=1}^{n-1}
\{V_i(F_{n-i,X}+F_{n-i}/X)+U_i Z_{b+\lambda_{n-i}}F_{n-i}\}
-\nu Z_{b+\lambda_{n-1}-D}Z_{b+\lambda_{n-1}}F_{n-1},
\]
\[
S_{zn}=\sum_{i=1}^{n-1}
\{V_i U_{n-i,X}+U_i Z_{c+\lambda_{n-i}}U_{n-i}\}
-\nu Z_{c+\lambda_{n-1}-D}Z_{c+\lambda_{n-1}}U_{n-1}.
\]

Then

\[
f_n=(0,0,0,2\xi P_n,2S_{\theta n}/\nu,
2[S_{zn}-2\eta X P_n/L]/\nu).
\]

The last source retains the pressure-gradient coupling. For arbitrary
manufactured lower and current profiles, the system defect equals

\[
(0,0,0,r_n/\xi,-2t_n/\nu,-2[z_n+\eta r_n/L]/\nu),
\]

where r_n,t_n,z_n are independently assembled full residual coefficients.
The comparison covers 216 cases at orders 1,2,3,5; it does NOT assume the
manufactured profiles solve any coefficient equation. The n=1 matrices
also reproduce the frozen first-order system. Omitting the new radial
axial-viscosity pressure source produces a nonzero defect.

## 3. Induction on one fixed collar, followed by the same global repair

Keep the same rectangle `0<=X<=a^2=Xa exp(tc/4)` and the same keep/cut
radii and reserved I3 patch from the
[annular checkpoint](../annular_first_correction/README.md).
Assume orders below n have been solved there, globally cut off and repaired.
Their radial cutoffs are eta-independent, and their moment bumps lie to
the right of this rectangle. All fixed lower-order derivatives are
therefore holomorphic on some common eta neighborhood for this finite
stage. That neighborhood may shrink with n. No uniform width or bound
over all n is asserted or required.

The leading profile remains fixed. On a smaller neighborhood there is a
finite bound C_n for A0_n,A1,f_n throughout the SAME radial interval.
The Volterra inverse has nonnegative diagonal weights, and adjacent
derivative factors vanish, even when one derivative hits the intervening
matrix. A length-k word has at most ceil(k/2) surviving eta derivatives.
The radial factorial therefore dominates the Cauchy loss, exactly as in
the first-order Picard proof. This constructs the unique analytic-eta
coefficient with zero axis traces on the same radial interval. Smoothness
at X=0 follows by the regular integral equation and parity, not by
estimating 1/xi on an annulus.

Cut off only F_n,U_n, add the five I3 bumps, and reconstruct M_n,V_n,Pi_n
forward from zero axis data. The five total moments are

\[
\int U_n,\quad \int RE_n,\quad \int(\Pi_n)_X,\quad
\int R\sum_{i+j=n}U_iE_j,\quad
\int[\sum_{i+j=n}U_iU_j-X(\Pi_n)_X].
\]

Every lower-order product and Omega_(n-1) is fixed while solving the
CURRENT five equations. Their variation is precisely the same linear B
map as before, because the current terms involve only U0=0 and the
power-law E0 on I3. Lower positive orders can be nonzero on I3; their
products go into the known discrepancy. Neither zero leading radial
memory nor zero lower-order profiles on I3 is assumed.

Invert B on the full finite incoming discrepancy, with the inherited
normalizations and inverse bound. This gives smooth bump coefficients
of any finite size, restores all five moments and keeps the zero axis
pressure. M_n=V_n=0 beyond the velocity patch; Pi_n can extend farther,
but vanishes beyond Xv because its radial source vanishes there and its
total third moment is zero. All positive-order fields have common compact
support within [0,Xv], and their velocities vanish on I4.

The conservative momentum identities close both total stress integrals.
At general n their non-viscous moment expressions are

\[
\frac{-k_\theta m_{n,2}+D\eta\partial_\eta m_{n,2}
+2l_\theta\eta m_{n,4}+d\partial_\eta m_{n,4}}{L},
\]
\[
\frac{-k_z m_{n,1}+D\eta\partial_\eta m_{n,1}
+2l_z\eta m_{n,5}+d\partial_\eta m_{n,5}}{L},
\]

with `(k_theta,l_theta,k_z,l_z)=(3/2-A+lambda_n,
3/2-2A+lambda_n,1-A+lambda_n,1-2A+lambda_n)`.
The lower viscous moments vanish by the previous order's exact moments;
at n=1 use the convergent SUBTRACTED leading angular moment. The tests
retain those lower viscous terms on fixtures where they are nonzero,
rather than arranging all test terms to cancel trivially.

The first stress T1 reaches the flat outer heat collar. For n>=2 its
residual source is zero beyond Xv, including its axial viscosity of the
positive order n-1. Its backward primitive thus has support in
`[Xkeep,Xv]`, strictly inside the leading active annulus. Consequently
every fixed derivative of T_n/zeta is bounded for n>=2. The first stress
keeps its already derived `zeta*delta^(-N_I)` bounds.

This completes the CONDITIONAL induction. It supplies finite constants
at each order and derivative, not a numerically materialized coefficient
sequence or effective uniform constants as n tends to infinity.

## 4. Finite residual powers improve without growing derivative losses

Let U^[J] denote the physical velocity, pressure and retained stress
through correction order J. The exact pressure equations cancel every
retained radial coefficient. The stress primitives cancel every retained
tangential coefficient. All remaining products and shifted viscosity
terms have increment at least `2h(J+1)`. On a fixed compact profile range,

\[
|\mathcal F_{\rm slow}(U^{[J]})|_m
\le C_{J,m}q^{2h(J+1)-K_m},\qquad K_m=3/2+2h+m.
\]

Here Fslow adds ONLY the prescribed tangential stress-divergence operator
to the momentum residual. One transverse derivative loses at most q^-1/2,
an axial derivative q^-D, and a time derivative q^-1. The coarser loss
q^-m suffices. Differentiation generates powers of the finite index J in
the constants, not an additional J-dependent loss in q. Cartesian
representatives handle the axis; inverse-radius stress coefficients occur
only on their common annular support.

For h=.005 and m=0, the first bound has exponent -1.49. Positive decay
in this finite augmented bound first requires J>=151. This is an exponent
example, not a claim that 151 actual coefficients were evaluated.

## 5. Sum streamfunctions, not velocity cutoffs

The raw coefficient series need not converge: C_n,m can grow arbitrarily
fast. Use smooth cutoffs chi(a_n q), equal to one for a_n q<=1/2 and
zero for a_n q>=1, with a_(n+1)>=2a_n. The physical Stokes streamfunction
and its Cartesian vector potential are

\[
S_n=q^{D+\lambda_n}M_n,\qquad
\mathcal A_n=\tfrac12q^{-A+\lambda_n}(M_n/X)(-x_2,x_1,0).
\]

Sum `curl(chi(a_n q) A_n)` for the radial and axial velocity. Sum the
azimuthal velocity, pressure and stress multiplied directly by chi.
The additional radial velocity is ESSENTIAL:

\[
r u_{rn}^{\rm cut}=q^{\lambda_n}
\left[\chi V_n-\frac{2\eta}{L}(a_nq)\chi'(a_nq)M_n\right].
\]

The other components are `uz_cut=chi uz_n`, `utheta_cut=chi utheta_n`.
The extra term is the curl product rule. Dropping it or reversing its
sign produces a nonzero divergence in the manufactured checks. Separate
derivatives of the Cartesian potential verify both velocity values and
their gradients, including at the axis and on both cutoff plateaus.

Every positive derivative of chi(aq) is supported where 1/2<=aq<=1.
Its factors a^k are compensated by the q^k from derivatives of q:
`|partial_z^alpha partial_t^beta chi(aq)|<=C q^(-D alpha-beta)`,
with constants independent of a. The tuple of coefficient potentials,
azimuthal velocity, pressure and stress has the sufficient loss
`ell_m=2A+m`. Including a curl and its cutoff product rule gives the
sufficient loss `ell'_m=2A+m+1` for the realized tuple.

### A scale choice controlling whole intervals, not endpoint samples

At step n collect the finite bounds needed through derivative order n,
including normalized profile seminorms and, for n>=2, weighted T_n/zeta
seminorms. A generic requirement is

\[
C(1+|\log q|)^P q^\gamma\le2^{-n},\qquad \gamma>0,
\quad 0<q\le a_n^{-1}.
\]

Pad C>=1 and P>=0. Put s=-log q. The elementary inequality
`log(1+s)<=sqrt(s)` follows because the derivative of
`u-log(1+u^2)` is `(u-1)^2/(1+u^2)>=0`. Young's inequality then gives

\[
\log C+P\log(1+s)-\gamma s
\le\log C+P^2/(2\gamma)-\gamma s/2.
\]

It suffices to choose log(a_n) larger than

\[
\max\left(\log(a_{n-1})+\log2,\ \log(q_0^{-1})+\log2,\ 1,
\max_{\rm requirements}\left[
\frac{2(\log C+n\log2)}\gamma+\frac{P^2}{\gamma^2}\right]\right).
\]

Only finitely many requirements occur at each stage. The actual background
uses gamma=nh for the half-order tail bounds, plus the finite normalized
requirements. The implementation keeps these scales as logarithms and
adds harmless strict padding. It does NOT evaluate actual Navier-Stokes
scales without the actual coefficient bounds. The analytic inequality
controls every smaller q; the finite high-precision samples merely check
the implementation of this choice.

For fixed m and J>=max(1,m), this yields

\[
|U-U^{[J]}|_m\le2^{-J}q^{h(J+1)-(2A+m+1)},
\qquad 0<q<(2a_J)^{-1}.
\]

On any compact region where q is bounded below, only finitely many terms
are nonzero. The sum is thus smooth for q>0, divergence-free and supported
as stated. No convergence of the uncut formal power series is claimed.

## 6. Compare with one FIXED truncation to obtain flat augmented error

In Cartesian fields Fslow is a fixed differential polynomial of order
two and degree two. Its stress coefficients have derivative bounds on
their common annular support, with sufficient loss `H_m=m+2` independent
of J. Subtracting its products one factor at a time gives

\[
|\mathcal F(U)-\mathcal F(U^{[J]})|_m
\le C_m q^{-H_m}|e|_{m+2}
(1+|U^{[J]}|_{m+2}+|e|_{m+2}),\quad e=U-U^{[J]}.
\]

For any desired derivative order m and decay power M, choose ONE finite
J>=m+2 so that

\[
2h(J+1)-(3/2+2h+m)\ge M+1,\qquad
h(J+1)\ge M+4A+3m+10.
\]

These deliberately conservative sufficient margins are checked with exact
rational arithmetic. Then keep J fixed and shrink the q-neighborhood.
The constant C_J,m in the partial sum can be absorbed by one extra power
q^-1. The finite residual and the tail comparison are both O(q^M).
The order of choices is (m,M), then J, then the neighborhood. Letting
J vary with q inside an unproved uniform bound would be invalid.

It follows conditionally that every fixed physical derivative of the
augmented background error is O(q^M) for every M. The finite constants
and cutoff scales can be enormous; no practical runtime for evaluating
this sequence is inferred from the existence argument.

## 7. What survives the whole background sum, and what does not follow

The extra finite normalized seminorm requirements give, on each fixed
enlarged compact annulus, all fixed derivatives in q*d_q, d_X and d_eta,

\[
q^A u_{\theta B}-E_0=O(q^{2h}),\quad
q^A u_{zB}-U_0=O(q^{2h}),\quad
q^{1/2}u_{rB}-V_0/R=O(q^{2h}),
\]
\[
|D^I(q^{A+1/2}T_{\rm phys}-T_0)|
\le C_I q^{2h}\zeta\delta^{-N_I}.
\]

For normalized expansion order K, separate a finite block through at
least 2K+2; its original powers are retained near q=0, while the half-order
tail supplies the full O(q^(2h(K+1))) asymptotic remainder. All correction
velocities and their cutoff commutators vanish on I4, and all positive-order
velocity and pressure terms vanish beyond Xv. The exact heat exterior is
preserved. Pressure is not asserted to vanish on I4.

Coefficient-by-coefficient moment cancellation does not, by itself, certify
exact nonlinear moment identities for the realized field after DIFFERENT
q cutoffs are applied to different orders. The augmented residual comparison
controls the resulting errors asymptotically. Later physical mean corrections
must check their constraints on the summed field itself.

At `X*=1/Lambda, eta=0`, the inherited bound `F0>.264/C` and the normalized
summation estimate give a finite B_F>=1 with
`|q^A utheta_B/R*-F0|<=B_F q^(2h)`. Therefore for

\[
0<q\le\min(1,[.132/(CB_F)]^{1/(2h)}),\qquad
u_{\theta B}\ge\frac{.132\sqrt{2/\Lambda}}{C}q^{-A}.
\]

This lower bound now survives ALL background corrections, not just the
first partial sum. It still excludes physical waves and later localization.
The bare momentum residual is

\[
\mathcal R(u_B,p_B)=-\operatorname{Div}_{\rm tangential}T_{\rm phys}
+\mathcal E_B,
\]

with E_B flat as above. The stress is still present and can be singular.
This is not a solution with an admissible smooth force. The two-component
operator is not the full divergence of a symmetric positive wave covariance.
The weighted edge estimates also do not prove uniform cone stability of
T0+small corrections: delta^-N can grow at a flat boundary.

## Reproduce and continue

The scalar test `q^2*f'-f+q+f^2=0` has positive formal coefficients
`1,2,8,44,...`, with `a_n>=(n-1)a_(n-1)`. Its uncut series diverges at
every positive q. The cutoff sum is locally finite and has improving
residual orders. Equal q powers are grouped before numerical evaluation
to avoid subtracting tiny almost-equal series values. This example tests
the summation mechanism and its cutoff derivative; it is not a fluid
simulation or evidence for Navier-Stokes blow-up.

```sh
python research/navier_stokes_cascade/results/background_summation/test_audit.py -v
python research/navier_stokes_cascade/results/background_summation/audit.py \
  --verify-record research/navier_stokes_cascade/results/background_summation/evidence.json \
  --output background_summation_reproduced.json
```

Next: verify the physical wave construction against this actual conditional
background, including the full covariance divergence, its signed stress
corrections, and all mean and cross interactions. Preserve radial pressure
memory and weighted edge costs at that interface. Smooth forcing through
the singular time, finite energy and a lower bound for the SAME complete
field remain separate unresolved obligations.
