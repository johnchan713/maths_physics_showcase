# Explicit coefficient cutoff, all fixed jets, and the pulse-gate boundary

Parent: `61ee5de6493e359c9a4e6056d5923ef9dc794d89`.
Status: `explicit-coefficient-cutoff-all-fixed-jets-and-pulse-gate-boundary`.

The compact/global reconciliation fixes the variance cutoff used in the
moving modulation loop.  The next coefficient and pulse estimates use a
different cutoff: the radial coefficient cutoff from Lemma 5.2, Step 1.
That cutoff had fixed plateaus and width but no formula.  This checkpoint
corrects the repeated implication that a quantitative shape was already
fixed, and specializes the construction to one explicit admissible member.

The same `C^-20` collar, `Xkeep`, `Xcut`, roots, gauges, supports, moments and
pressure memory are retained.  The selection supplies every fixed cutoff
jet, including the third jet needed by the second coefficient.  It does not
evaluate the original first physical band, finalize `H2,3`, or prove an
all-family pulse, covariance/mean, residual or blow-up theorem.  Historical
frozen evidence is preserved.

Source: Lemma 5.2, Step 1 of the retained
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA-256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The source's final theorem is not used as a premise.

## 1. Select the previously qualitative coefficient cutoff

Let

\[
 \phi(s)=\exp[-s^{-1}-(1-s)^{-1}],\qquad 0<s<1,
\]

extended by zero, and put

\[
 K(s)=1-\frac{\int_0^s\phi(r)\,dr}{\int_0^1\phi(r)\,dr},\qquad
 \kappa(X)=K\!\left(\frac{\log(X/X_{\rm keep})}{\Delta}\right),
 \quad \Delta=\frac{t_c}{16},\quad t_c=C^{-20}.             \tag{1}
\]

Thus `kappa=1` through `Xkeep`, `kappa=0` from
`Xcut=Xkeep exp(Delta)`, `0<=kappa<=1`, and it is smooth, decreasing and
flat at both ends.  This is a legal specialization of the original
qualitative choice, not an identification with the already selected
variance cutoff.  The two cutoffs may use the same normalized bump formula
without becoming the same function: their arguments, windows and roles are
different.

## 2. Exact normalized derivative bounds at every fixed order

Set `u=1/s`, `v=1/(1-s)`.  If
`phi^(n)=P_n(u,v) phi`, direct differentiation gives a finite integer
polynomial recurrence.  Since

\[
 u^a e^{-u}\le a!,\qquad v^b e^{-v}\le b!,
\]

the sum of `|coefficient| a! b!` bounds `|phi^(n)|`.  On
`[1/3,2/3]`, `phi>=exp(-9/2)`, so

\[
 \left(\int_0^1\phi\right)^{-1}<3^6=729.              \tag{2}
\]

The exact integer majorants for `phi` through order seven are

\[
 1,4,80,3600,292320,37704960,7139819520,1867003084800.
\]

Consequently the selected cutoff has the strict normalized bounds

\[
 L_1<729,qquad L_2<2916,qquad L_3<58320,              \tag{3}
\]

and the evidence records explicit bounds through order eight.  The same
recurrence supplies a finite constant at every fixed order.  No constant
uniform in derivative order is asserted or needed by the fixed-order Borel
and pulse statements.

## 3. Physical radial jets keep every logarithmic Jacobian term

With `D=partial_(log X)`, the exact identity is

\[
 \partial_X^m\kappa
 =X^{-m}\prod_{i=0}^{m-1}(D-i)K
 =X^{-m}\sum_{j=1}^m s(m,j)\Delta^{-j}K^{(j)},         \tag{4}
\]

where `s(m,j)` are the signed first-kind Stirling numbers.  Replacing
`X^-m` by `Xkeep^-m` and taking absolute values gives a reproducible bound
for every fixed physical radial jet.  In particular

\[
 |\kappa_{XXX}|<X_{\rm keep}^{-3}
 \left(58320\Delta^{-3}+8748\Delta^{-2}+1458\Delta^{-1}\right). \tag{5}
\]

The last two terms in (5) are real logarithmic-coordinate costs; omitting
them would repeat the error isolated by the higher-jet interface.

Equations (3)--(5) discharge the cutoff-shape inputs `M,N,L_j` in the
first-coefficient, second-inner and finite-head frame estimates.  They do
not by themselves bound the inner coefficient, its global moment repair or
the physical pulse operators.

## 4. The first unresolved scale inequality is now shape-explicit

Let `Cstar(M,N)` be the complete affine coefficient bound from the sharp
cutoff/current-scale checkpoint, retaining pressure, shifted axial viscosity
and both derivative rows.  The safe selected value is

\[
 C_{\rm coeff}=C_\star(729,2916).                       \tag{6}
\]

For the same original tilt and a supplied raw-deficit lower bound `d0`, put

\[
 \gamma_B=\frac{\omega^2d_0}{9\kappa_u},\qquad
 H_\ell=\frac{\gamma_B}{8}-C_m\ell^{-6},\qquad
 S_\ell=2^{2h(\ell-2)}.
\]

The existing sufficient mixed certificate at the actual first band requires

\[
 H_\ell>0,qquad
 C_{\rm coeff}\le S_\ell\min\left\{
 \frac{r_b}{8},\frac{\omega^2d_0}{16L_\Psi},
 \frac{H_\ell}{8C_p}\right\},                         \tag{7}
\]

together with the independent chart and endpoint tests.  Smoothness or the
new derivative bounds do not supply `d0` or those endpoint constants.  Thus
(7), not an unspecified cutoff norm, is the precise first-scale obstruction.
For an existential later tail, finite `Ccoeff` can be absorbed by taking a
later band once the other positive inputs are established; that does not
retroactively verify the historical first band.

## 5. What remains before `H2,3` and uniform pulses

The nested second-core checkpoint gives pre-cut jets `J2,r`.  The selected
cutoff now makes the complete product row explicit:

\[
 \widehat J_{2,r}
 \le\sum_{j=0}^r {r\choose j}k_jJ_{2,r-j},\qquad r\le3, \tag{8}
\]

with `k_j` from (4).  The compact/global reconciliation also supplies the
leading fourth angular compact input relative to its retained premises.
The next missing coefficient datum is the finalized global first axial
coefficient norm `K1,3=||U1||_(C_eta^3)` on the full support.  It is needed
in the shifted-viscosity moments `Q10,Q11`; those totals then enter the exact
inverse-800 second five-row repair.  Only after adding the explicit repair
bump jets to (8) is the finite-head constant

\[
 H_{2,3}\ge\max_{0\le r\le3}
 \|\partial_X^r(E_2,U_2)_{\rm finalized}\|             \tag{9}
\]

available for the current-frame derivative estimates.  Local collar bounds
or exact leading moment cancellation cannot replace the full-support data.

Therefore this checkpoint establishes a concrete cutoff and all its fixed
jets, but does **not** promote the historical primary-pulse theorem to an
unconditional or complete pulse-family theorem.  The next obligation is to
derive `K1,3`, close the full second repair and `H2,3`, verify (7) or construct
a distinct certified later tail, and then replay the defining/inverse/gap/
support estimates before the complete covariance, axial-pressure mean and
viscous/cutoff/interaction residual joins.

The diagnostic suite checks the exact polynomial recurrence, the Stirling
chain rule through order eight, frozen dependency hashes, and scope guards.
These computations reproduce the stated bookkeeping; they are not numerical
evaluation of the enormous profile or a continuum proof of retained premises.

Reproduce with:

```sh
python3 test_audit.py -v
python3 audit.py --verify-record evidence.json --output fresh.json
```

The output path must be fresh; frozen evidence is never overwritten.
