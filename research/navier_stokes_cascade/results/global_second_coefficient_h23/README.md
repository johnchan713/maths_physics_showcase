# Full-pressure-memory second repair and explicit global `H2,3`

Parent: `5b875cd821b89c02f23939346c7fa0440c2cf33a`.
Status: `shifted-viscosity-second-repair-and-global-H23-assembled`.

The parent finalized the whole-support first axial coefficient bound

\[
 K_{1,3}\le K_{1,3}^*=J_{1,3}+320C_{1,3}.
\]

This checkpoint inserts that bound and the already proved
`K0,5 < H^262` into the complete order-one signed pressure-source formulas,
assembles the unchanged second five-row inverse, and adds every radial jet
through order three of the actual repair bumps.  Relative to the retained
analytic-core, collar, cutoff, support and inverse premises, it gives a
finite explicit selected-parameter certificate `H2,3*`.

This is not an unconditional Navier--Stokes result.  The nested Picard
exponentials and all small geometric denominators remain visible rather than
being replaced by a guessed power of `H`.  No uniform physical-pulse family,
complete covariance mean, residual improvement or blow-up statement is made.
Historical frozen evidence is unchanged.

## 1. The complete full-support signed source

Let `K05=H^262` and retain the parent's explicit `K13*`.  Since the
factorial angular norms are nested, `K05` bounds every leading `C_eta^r` row
needed below for `r<=5`.  The regular-average estimate gives

\[
 v_{i,r}\le R_r[4+4(D+2hi)+4(r+1)]K_{i,r+1},\qquad r=2,3,
\]

where the exact positive rational expressions `R2,R3` are retained.
Consequently

\[
 u_0=K_{0,5},\quad u_1=K_{1,3}^*,\quad
 v_0=v_{0,2},\quad v_1=v_{1,2},\quad p=v_{0,3}
\]

are valid inputs to the original whole-support `k=1` moment theorem.  Put

\[
 B=u_0v_1+u_1v_0,\qquad C=2v_0v_1,
\]

and, for `j=0,1`,

\[
 A_{1j}=\nu R_1R_2(6+4|j-D|)(12+4j)
              {X_v^{j+1}\over j+1}\,p.
\]

Then the complete signed bounds are

\[
\begin{aligned}
 Q_{10}&=R_1X_v[(2h+4D)v_1+(6+4|2h-A|)B]
                    +X_vC/2+2\nu v_1+A_{10},\\
 Q_{11}&={R_1X_v^2\over2}[(1+2h+4D)v_1+(6+4|2h+D|)B]
                    +X_v^2C/4+A_{11}.                 \tag{1}
\end{aligned}
\]

Both shifted axial-viscosity terms `A10,A11` remain.  In particular the
weighted radial-viscosity cancellation does not erase `A11`.  Both ordered
cross terms `(0,1)` and `(1,0)` occur in `B`.  The integrals end at the true
common support `Xv`, not at the coefficient cutoff or either reserved
interval; all inherited pressure memory is therefore retained.

## 2. Do not identify the finalized `F1` and `U1` bounds

The tentative first fields share the parent's `J13`, but the exact repair
uses two additive `U` bumps and three relative `E` bumps.  If `r_p` is the
proved lower bound for `sqrt(2X)` on those supports, disjointness gives

\[
 u_1^*=J_{1,3}+320C_{1,3},\qquad
 f_1^*=J_{1,3}+320C_{1,3}\max(1,r_p^{-1}).             \tag{2}
\]

The second expression is needed in the fixed lower-order squares and
products.  Replacing it by `K13*` would silently discard the relative
`E`-to-`F` inverse-radius cost.

## 3. Close the second collar envelope, then cut with all Jacobian terms

There is no free `Tcollar` input.  Let `J1,r` be the explicit whole-inner
first-field `C_eta^3` rows from the first-axis/collar checkpoint.  Apply the
same cutoff product rule and the actual first repair jets (10)--(11) below to
obtain finalized `F1,U1` radial rows through order three.  If `H0,mix` is the
already proved leading whole-support mixed-`J3` bound (the retained result is
`H0,mix<H^70`), set

\[
 T_{col}=\max\{1,H_{0,mix},
        \max_{r\le3}(\|\partial_X^rF_1\|_3,
                     \|\partial_X^rU_1\|_3)\}.        \tag{3}
\]

This is a fixed selected-parameter expression using `C13`, not a qualitative
smoothness constant.  The nested second-core/collar checkpoint then gives

\[
 J_{2,r}=A_{c,r}M_{2c}
       +T_{col}^{p_r}(A_{a,r}M_{2a}+B_{a,r}),
 \qquad (p_0,p_1,p_2,p_3)=(0,0,2,4),                 \tag{4}
\]

where `M2c` and `M2a` are its displayed nested Picard exponentials.  They
remain explicit and enormous; neither is evaluated at `Tcol=1`.

For the selected coefficient cutoff

\[
 \kappa(X)=K\left({\log(X/X_{keep})\over\Delta}\right),
 \qquad \Delta=t_c/16,
\]

the fixed normalized bounds are `(L0,L1,L2,L3)=(1,729,2916,58320)`.
The physical bounds used here are

\[
\begin{aligned}
 k_0&=1,\\
 k_1&={L_1\over\Delta X_{keep}},\\
 k_2&={L_2\Delta^{-2}+L_1\Delta^{-1}\over X_{keep}^2},\\
 k_3&={L_3\Delta^{-3}+3L_2\Delta^{-2}+2L_1\Delta^{-1}
       \over X_{keep}^3}.                              \tag{5}
\end{aligned}
\]

Thus the tentatively cut jets are

\[
 P_{2,r}=\sum_{j=0}^r {r\choose j}k_jJ_{2,r-j},
 \qquad 0\le r\le3.                                  \tag{6}
\]

In particular neither the `k3 J2,0` term nor the lower logarithmic-Jacobian
terms in `k3` are omitted.

## 4. Five discrepancies and the exact second inverse

Using leading inner bounds `f0,u0`, the separate values (2), and (1), the
ordinary second-order discrepancies are

\[
\begin{aligned}
 \mu_1&=X_{cut}P_{2,0},\qquad \mu_2=X_{cut}^2P_{2,0},\\
 \mu_3&=2f_0X_{cut}P_{2,0}+X_v(f_1^*)^2+Q_{10}/2,\\
 \mu_4&=X_{cut}^2(u_0+f_0)P_{2,0}+X_v^2u_1^*f_1^*,\\
 \mu_5&=2u_0X_{cut}P_{2,0}+f_0X_{cut}^2P_{2,0}
       +X_v(u_1^*)^2+{X_v^2(f_1^*)^2\over2}+Q_{11}/2. \tag{7}
\end{aligned}
\]

The lower `F1^2`, `U1^2` and `U1 F1` terms are fixed inputs, not a quadratic
current-order repair.  The same linear moment map gives

\[
 C_{2,1}=800\max\left\{{\mu_1\over\rho_3},
 {\mu_1/\rho_3+K_g(4/e_*)\mu_4\over\lambda},
 K_g\mu_2,{(4/e_*)\mu_5\over\rho_3},{4\mu_3\over e_*}\right\}. \tag{8}
\]

There is no nonlinear root in (8).

## 5. Add the actual repair jets

For the unchanged width-`1/10` logarithmic bumps, the analytically proved
bounds through the required order are

\[
 (B_0,B_1,B_2,B_3)
 \le(320,121600,220000000,18870643200000).             \tag{9}
\]

With `p` a common lower `X` bound on their supports, the `U`-bump jets are

\[
 \left(B_0,{B_1\over p},{B_2+B_1\over p^2},
 {B_3+3B_2+2B_1\over p^3}\right),                    \tag{10}
\]

while the relative `E` bumps contribute to `F` the jets

\[
 {1\over r_p}\left(B_0,{B_1+B_0/2\over p},
 {B_2+2B_1+3B_0/4\over p^2},
 {B_3+(9/2)B_2+(23/4)B_1+(15/8)B_0\over p^3}\right). \tag{11}
\]

Disjoint support again avoids a factor of five.  First form

\[
 U_{2,r}^*=P_{2,r}+C_{2,1}U_r,\qquad
 F_{2,r}^*=P_{2,r}+C_{2,1}F_r.                       \tag{12}
\]

The finite-head norm is defined for `(E2,U2)`, not `(F2,U2)`.  On the
common positive support retain

\[
 0<R_{lo}\le\sqrt{2X}\le R_{hi},\qquad
 (R_0,R_1,R_2,R_3)=(R_{hi},R_{lo}^{-1},R_{lo}^{-3},3R_{lo}^{-5}).
\]

The product rule for `E2=sqrt(2X)F2` gives

\[
 E_{2,r}^*=\sum_{j=0}^r{r\choose j}R_jF_{2,r-j}^*.   \tag{13}
\]

Combining (6), (8), and (10)--(13) yields the explicit certificate

\[
 \boxed{H_{2,3}\le H_{2,3}^*:=
 \max_{0\le r\le3}\{U_{2,r}^*,E_{2,r}^*\}.}         \tag{14}
\]

Every term in (14) is a fixed selected-parameter expression supplied here
or by a frozen dependency.  It can be enormous; finiteness and correct
dependency closure, not a cosmetic `H` exponent, are the conclusion.

## 6. Limits and next obligation

Equation (14) closes the coefficient-level second repair relative to the
retained foundations.  It does not by itself verify the historical first
physical-band inequality involving the raw deficit `d0`, nor does it give
uniform actual-pulse estimates across every physical scale.

Next, insert `K05`, `K13*`, and `H23*` into the physical pulse maps and prove
uniform value, derivative and support bounds on every retained scale.  Then
audit the complete-wave covariance corrections, including axial velocity and
both exceptional pressure tails, before assembling viscosity, cutoff and all
interaction terms in the full residual.

Reproduce with:

```sh
python3 test_audit.py -v
python3 audit.py --verify-record evidence.json --output fresh.json
```

The output must be fresh.  Passing tests verifies the exact bookkeeping and
frozen dependency chain; it is not an independent continuum proof of the
retained analytical premises.
