# Q15: fixed-width flat-edge transfer to closed-support seminorms

Parent: `f80034adc47312dcddf0e09c147880b6244cf18e` (Q14).
Status: `fixed-width-weighted-to-unweighted-transfer-proved`.

The compact/slow reconstruction uses the unchanged logarithmic-shell width
w=C_setup^-20, with 0<w<=1, and the exact weight

    zeta=exp(-w^2/x^2-4/(L-x)^2), delta=min(1,x,L-x), 0<x<L.

Q5 supplies weighted source/operator statements at finite derivative and
pole requests. Q14 corrects how those requests reach the primitive domains.
This checkpoint proves the explicit edge-to-unweighted implication required
to use such statements in Q13's closed-support representative bound. It does
not assume a positive absolute lower bound on zeta; Q2 already ruled that out.

## 1. Uniform majorant for every finite pole

For any real p,s>=0, let n=ceil((p+s)/2). Then, for every L>0 and 0<x<L,

    zeta delta^(-p) <= n! w^(-2n) delta^s.                 (1)

If delta<1, either x=delta or L-x=delta. In the first case the inner
exponent is w^2/delta^2; in the second the outer exponent is
4/delta^2>=w^2/delta^2. Thus zeta<=exp(-w^2/delta^2).
For u=w^2/delta^2, the exponential series gives exp(u)>=u^n/n!, hence

    zeta delta^(-p-s)
      <= n! w^(-2n) delta^(2n-p-s) <= n! w^(-2n).

If delta=1, use zeta<=1 and n! w^(-2n)>=1 instead. It would be incorrect
to assert zeta<=exp(-w^2/delta^2) in that interior case. This separate
case makes (1) valid without an upper bound on the shell length L.
The case n=0 is included. No factorial or tiny root parameter needs numerical
evaluation: the expression in (1) is a finite symbolic constant at each
fixed request, even at Q14's enormous pole domains.

## 2. Consequence for the actual weighted jet contracts

Fix stage j and finite jet request R. Suppose the inherited weighted
estimate has the form, on all accepted bands ell,

    |D^alpha F_(j,ell)|
      <= H_(j,R) epsilon_ell^b (1+|log epsilon_ell|)^P
                    zeta delta^(-p_alpha),    |alpha|<=R,

with finite H, P and finite nonnegative p_alpha, independent of ell.
For any requested edge vanishing power s, (1) gives the explicit normalized
closed-support bound

    sup_(ell,x) epsilon_ell^(-b) (1+|log epsilon_ell|)^(-P)
                   delta^(-s) |D^alpha F_(j,ell)|
      <= H_(j,R) n_alpha! w^(-2*n_alpha),
    n_alpha=ceil((p_alpha+s)/2).                           (2)

Taking the maximum over the finite jet list is legitimate. The constants
need not be uniform in j,R,p or s. The fixed width is independent of physical
band; (2) therefore introduces no new band-dependent constant and spends
no epsilon/q exponent reserve. The logarithmic degree is retained, not
silently replaced by a pure-power bound. This applies to every weighted
term, including the complete axial/moving-pressure terms, when its exact
weighted hypothesis is supplied; it does not erase old pressure tails or
replace total pressure by a wave-only quantity.

If the all-fixed-order version of the hypothesis holds, all normal and
mixed jets tend to zero at both edges faster than every fixed power.
Their zero extensions are smooth in the normalized support charts: induct
on derivative order, using the vanishing traces and difference quotients
across each boundary. This statement uses smooth support charts and the
jets in those charts; it is not a claim that arbitrary moving boundaries
or physical coordinate changes have bounded derivatives. Those transformations
must still use their original derivative and support contracts. For finite
R, conclude only the corresponding finite regularity, not C-infinity.

## 3. Why fixed width and actual weighted hypotheses matter

The constant cannot be uniform over widths tending to zero. Set L=2,
w_ell=1/ell, x_ell=w_ell, ell>=1. Then delta=x_ell and

    zeta delta^(-p) >= exp(-5) ell^p.

For p>0 this diverges. This counterfamily does not describe the selected
field: its width is fixed. It shows why the width quantifier in the provider
is indispensable. Likewise, finite seminorms separately on each band do not
supply the uniform H in (2). The all-band weighted source hypothesis must
be established before applying (1).

## 4. Join to Q13 and remaining obligation

For each complete representative Z_j=(A_j,b_j e_theta,pi_j), the missing
physical bound has the form

    |D^k Z_j| <= C_(j,k)(1+|log q|)^P_(j,k)
                         q^[h*delta_gain*j-(2A+2k)].

For its weighted compact pieces, (2) removes the flat-edge pole obstruction
once the original weighted estimates and physical coordinate transfer are
available at Q14's fully expanded requests. It also justifies their flat
zero extensions. It does not derive the q exponent, prove the phase/potential
transfer, or show that an unweighted background pressure component has the
weight zeta. Those terms must remain in their separate complete pressure
and representative estimates. In particular a velocity/curl estimate alone
still does not control the fixed-gauge potential used in the cutoff.

Next: give a term-by-term map from the actual sourced/signed wave potentials,
mean potentials, direct swirl and full pressure increments to their weighted
or unweighted physical estimates, with the same root h, fixed supports and
Q14 domains. This checkpoint removes a specific edge-transfer step, not the
remaining representative theorem. Endpoint closure, unconditional residual
improvement and blow-up are not claimed.

## 5. Reproduction

The proof is (1)--(2) and the all-order extension argument. Exact rational
checks exercise the exponent majorant, interior case, width degeneration and
symbolic large-order handling; tests are not a continuum proof. Q14 evidence
is regenerated and all inherited hashes are checked without editing frozen
records. On 29 October 2026 Q14's focused CI succeeded; its broad workflow
was still running at inspection.

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q15_flat_edge_transfer/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q15_flat_edge_transfer/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q15_flat_edge_transfer/evidence.json --output q15_fresh.json
```
