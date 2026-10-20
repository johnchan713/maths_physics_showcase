# Global first-coefficient C3 boundary and the missing leading fifth jet

Parent: `df9bd482678a6c9cf502233bbc080c455c902630`.
Status: `local-first-C3-and-exact-global-leading-C5-repair-boundary`.

The explicit coefficient cutoff closes the shape part of the first and
second coefficient estimates.  This checkpoint audits the next stated input,
the finalized global first axial coefficient norm `K1,3=|U1|Ceta3`.  It
extends the existing complex-core/collar estimate from angular order two to
three and follows that estimate through the selected eta-independent cutoff.
It then retains the complete signed leading pressure memory and the actual
linear inverse-800 repair.

The audit finds a one-derivative gap before that repair can be quantified:
the signed leading source moments in `C_eta^3` require the finalized leading
axial field on its **whole support** in `C_eta^5`.  The current chain gives a
PREMOD native sixth-angular bound, a FINAL compact fifth-angular bound, and a
later axial/reference bound only through fourth angular order.  Those are not
a whole-support fifth-order certificate.  This is an input boundary, not
evidence that the selected profile fails.

No width, cutoff, support, root, frequency, pressure gauge or historical
evidence is changed.  The estimates are analytical implications of the
retained continuum construction.  The exact tests below check algebra and
bookkeeping; they do not prove the inherited construction or a Navier--Stokes
singularity.

## 1. The whole-inner first coefficient already has a local C3 bound

Keep the notation of `axis_first_coefficient`: `delta=rho/(800S)`, and let
`Mcore` and `Mcollar` be its two distinct holomorphic solution bounds.  The
existing proof bounds both solutions on a `3delta/2` strip and uses disks of
radius `delta/2`.  The same Cauchy argument at one higher angular order gives

\[
 c_3=1+{2\over\delta}+{4\over\delta^2}+{8\over\delta^3}.
\]

If `(p_r,q_r)` are the established positive-radius collar coefficients,
then, for `0<=r<=3`,

\[
 \sup |\partial_X^r(F_1,U_1)|_{C_\eta^3}
 \le c_3\{r!(2\Lambda)^rM_{core}+p_rM_{collar}+q_r\}.       \tag{1}
\]

This uses exactly the same radial cover as the earlier C2 result.  It does
not assert radial analyticity across either flat join.  Multiplication by the
selected coefficient cutoff is eta-independent, so with its physical jets
`k_j` the exact product estimate is

\[
 \widehat J_{1,r}^{(3)}\le
 \sum_{j=0}^r {r\choose j}k_jJ_{1,r-j}^{(3)},\qquad r\le3.   \tag{2}
\]

The explicit-cutoff checkpoint supplies every `k_j`; in particular the
third row retains `k_3 J_{1,0}^{(3)}`.  Equations (1)--(2) therefore close
the local/cutoff part of the `C_eta^3` estimate.  The disjoint I3 repair is
the remaining part.

## 2. Why the exact global repair needs leading U0 in C5

Write the regular leading radial velocity as

\[
 v={2\eta U-2D\eta\bar U-(1-\eta^2)\bar U_\eta\over L},
 \qquad L=1-2h\eta^2,
\]

and use the full-support signed moment formulas already proved in
`signed_pressure_source_moments`:

\[
 I_0=2\nu v(0)+{D\eta(a_0)_\eta+d(b_0)_\eta-2A\eta b_0\over L}
       +{c_0\over2},                                      \tag{3}
\]
\[
 I_1={-a_1+D\eta(a_1)_\eta+d(b_1)_\eta+2D\eta b_1\over L}
       -{c_1\over2}.                                      \tag{4}
\]

Here `a_j=int X^j v`, `b_j=int X^j Uv`, and `c_j=int X^j v^2` through
the true common support `Xv`.  To bound (3)--(4) in `C_eta^k`, their
explicit eta derivatives require `a_j,b_j` in `C_eta^(k+1)`.  Therefore
`v` is needed in `C_eta^(k+1)`.  The last term in the numerator defining
`v` contains `bar U_eta`, so this requires `U` in `C_eta^(k+2)`:

\[
 \boxed{|(I_0,I_1)|_{C_\eta^k}\quad\hbox{requires}\quad
        \sup_{0\le X\le X_v}|U_0(X,\cdot)|_{C_\eta^{k+2}}.} \tag{5}
\]

For the desired first repair `k=3`; hence the missing datum is

\[
 K_{0,5}:=\sup_{0\le X\le X_v}|U_0(X,\cdot)|_{C_\eta^5}<\infty
 \quad\hbox{with an explicit selected-parameter majorant}. \tag{6}
\]

This corrects the prior roadmap's suggestion that the newly joined fourth
leading derivative was already enough to derive `K1,3`.  Fourth leading
order is enough for the later second-source calculation once `K1,3` is
known; it is not enough to construct the first repair itself in third order.

## 3. Complete conditional bound, with pressure memory retained

The following makes (6) the precise boundary rather than an unspecified
regularity request.  In every factorial norm of order at least two,

`|2h eta^2| <= 8h`, so for `0<h<=1/100`,

\[
 R:=|L^{-1}|\le(1-8h)^{-1}\le25/23.                    \tag{7}
\]

If (6) holds with majorant `K`, the regular-average formula gives

\[
 |v|_{C_\eta^4}\le R(4+4D+20)K<29K=:V.               \tag{8}
\]

Product algebra, (3)--(4), and positive radial integration then give the
explicit full-support signed bounds

\[
 Q_{03}=2\nu V+RX_v\{8DV+(16+4A)KV\}+{X_vV^2\over2}, \tag{9}
\]
\[
 Q_{13}=R{X_v^2\over2}\{(1+8D)V+(16+4D)KV\}
                         +{X_v^2V^2\over4}.           \tag{10}
\]

These are bounds for the signed totals after integration.  They are not
bounds for `int |w0|`, and (9) keeps the nonzero axis-viscosity term.  The
weighted radial-viscosity boundary cancels in (10), exactly as before.

Let `J` bound the tentatively cut first fields in `C_eta^3`, and let `f,u`
bound the leading `F0,U0` values in the same norm.  The five ordinary I3
discrepancies satisfy

\[
 (\mu_1,\ldots,\mu_5)\le
 (X_cJ,X_c^2J,2fX_cJ+Q_{03}/2,
  X_c^2(u+f)J,2uX_cJ+fX_c^2J+Q_{13}/2).              \tag{11}
\]

For angular order three the actual reciprocal-amplitude cost is
`W3=5/e_*`.  The unchanged linear repair therefore gives

\[
 |c_1|_3\le800\max\left\{
 {\mu_1\over\rho},
 {\mu_1/\rho+K_gW_3\mu_4\over\lambda},
 K_g\mu_2,{W_3\mu_5\over\rho},W_3\mu_3\right\}.      \tag{12}
\]

Adding the fixed I3 bump `C_eta^3` operator norm to (2) and (12) produces
the desired finalized `K1,3`.  There is no high-order nonlinear root here:
the positive-order solve is the same exact linear map.  Complete future
pressure memory enters only through (9)--(11), and is not reset at the
cutoff, I3 or I4.

## 4. The current fourth-order statement cannot imply (6)

As an abstract norm control, put

\[
 f_n(\eta)=n^{-4}\sin(n\eta).
\]

For all `n>=1`, `|f_n|Ceta4<=65/24`, while
`|f_n^(5)(0)|/5!=n/120`.  Thus no bound depending only on a `C_eta^4`
majorant controls the missing fifth derivative.  This is not a family of
selected PDE profiles and does not show that their fifth derivative is
large.  It proves only that the missing row must come from the defining
maps (or another structural identity), not from norm monotonicity.

The retained all-order smoothness theorem does give qualitative finiteness
at every fixed order.  That is enough for an existential later diagonal
choice, but not for the selected-parameter quantitative repair (9)--(12),
the historical first band, or a stated finite `H2,3` constant.

## 5. Status and next obligation

Established here:

- the local whole-inner first coefficient has an explicit `X^3 C_eta^3`
  bound, and the selected cutoff propagates it without angular loss;
- the complete signed pressure-memory derivative ladder is `k -> k+2`;
- equations (7)--(12) reduce the global first repair to one explicit
  whole-support leading input `K0,5`, together with already named geometry,
  leading swirl and fixed bump constants.

Not established here:

- a selected-parameter whole-support bound for `K0,5`;
- the resulting finalized `K1,3`, second signed moments `Q10,Q11`, second
  repair or `H2,3`;
- the original raw-deficit first-band inequality, uniform physical pulses,
  full covariance/mean corrections, the complete residual, or blow-up.

The next obligation is to extend the later finalized leading axial/reference
stage from its current fourth-angular estimate to fifth order and join it to
the already available native and FINAL compact fifth-order rows.  Then insert
that whole-support `K0,5` in (9)--(12), close `K1,3`, and feed it with the
leading fourth-order input to the shifted-viscosity `Q10,Q11` calculation.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/first_coefficient_c3_boundary/test_audit.py -v
python3 research/navier_stokes_cascade/results/first_coefficient_c3_boundary/audit.py \
  --verify-record research/navier_stokes_cascade/results/first_coefficient_c3_boundary/evidence.json \
  --output first_coefficient_c3_boundary_fresh.json
```

The output path must be fresh.  Frozen evidence is never overwritten.
