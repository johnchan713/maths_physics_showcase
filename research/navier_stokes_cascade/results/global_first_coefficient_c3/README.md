# Finalized whole-support first-coefficient third-angular bound

Parent: `6ef0567c8ef3373735d956dc0fec59e864878205`.
Status: `full-pressure-memory-first-repair-and-global-K13-assembled`.

The parent supplied the last missing leading input,

\[
 K_{0,5}=\sup_{0\le X\le X_v}|U_0(X,\cdot)|_{C_\eta^5}<H^{262}.
\]

This checkpoint substitutes that bound into the complete signed leading
pressure-source totals, combines them with the already proved local
first-coefficient and selected cutoff bounds, applies the unchanged linear
inverse-800 moment repair, and obtains one explicit finite bound for the
finalized first axial coefficient:

\[
                    \boxed{K_{1,3}\le K_{1,3}^*:=J_{1,3}+320C_{1,3}.} \tag{1}
\]

All quantities on the right are displayed below in the unchanged selected
parameters.  This is an analytical implication of the retained continuum
foundations, not a numerical evaluation of their enormous constants.  The
tests check algebra and dependency bookkeeping; they are not a continuum
proof.  No second repair, uniform physical pulse family, covariance-mean
closure, complete residual improvement or blow-up is claimed.

## 1. The local/cutoff value constant

Let `delta=rho_axis/(800S)`.  The complex-core and activated-collar Picard
majorants retained by the parent are

\[
 M_c=(z_c+6z_c^2/\delta)e^{6z_c^2/\delta},\qquad
 M_a=(z_a+12z_a^2/\delta)e^{12z_a^2/\delta}.           \tag{2}
\]

Here `z_c=2C_core b_core` uses the explicit complex-core row/source envelope,
and `z_a=2C_col b_col` uses the explicit actual activated-collar envelope.
The first exponent uses the core strip loss `delta/2`; the second uses the
actual collar choice `Delta=delta/4` and `e<3`.  Cauchy on the remaining
strip and the already proved core/collar join give

\[
 J_{1,3}:=(1+2\delta^{-1}+4\delta^{-2}+8\delta^{-3})(M_c+M_a). \tag{3}
\]

This bounds both tentatively cut fields `F1,U1` in `C_eta^3`.  The selected
coefficient cutoff is eta-independent and between zero and one, so its
radial jets `729,2916,58320` are needed for mixed/radial estimates but do
not enlarge the value norm (3).  This distinction prevents a spurious
third radial-cutoff charge in `K1,3`.

## 2. Full pressure memory and the five repair targets

Use `K=H^262` for `K0,5`; it also bounds the leading `U0` `C_eta^3` value.
The previously joined leading swirl bound supplies

\[
 f_3=\max\{K_{axis},2(1+X_a^{-1})H^{70}\}.             \tag{4}
\]

For `h<=1/100`, put `R=(1-8h)^-1<=25/23`, `A=1/2+h`,
`D=1/2-h`, and

\[
 V_4=R(4+4D+20)K<29K.                                  \tag{5}
\]

Integration over the **entire unchanged support** `[0,Xv]` gives

\[
 Q_{03}=2\nu V_4+RX_v\{8DV_4+(16+4A)KV_4\}+X_vV_4^2/2, \tag{6}
\]
\[
 Q_{13}=RX_v^2\{(1+8D)V_4+(16+4D)KV_4\}/2+X_v^2V_4^2/4. \tag{7}
\]

The first term in (6) is the nonzero axis-viscosity boundary.  There is no
corresponding term in (7).  Neither integral stops at the coefficient cutoff,
`I3` or `I4`; consequently later pressure memory is not reset.

With `Xcut` the selected coefficient-cutoff endpoint, the five ordinary
discrepancy bounds are

\[
\begin{aligned}
\mu_1&=X_{cut}J_{1,3},&\mu_2&=X_{cut}^2J_{1,3},\\
\mu_3&=2f_3X_{cut}J_{1,3}+Q_{03}/2,&
\mu_4&=X_{cut}^2(K+f_3)J_{1,3},\\
\mu_5&=2KX_{cut}J_{1,3}+f_3X_{cut}^2J_{1,3}+Q_{13}/2.
\end{aligned}                                           \tag{8}
\]

## 3. Exact linear repair and the finalized bound

Keep the actual positive patch parameters `rho3`, `lambda_patch`, `e_*`
and `K_g=(sqrt(2)rho3^(3/2))^-1`.  Because
`e(eta)^-1=(1+eta^2)/e_*`, its complete third-angular factorial norm is
`W3=5/e_*`.  The established matrix estimate therefore gives

\[
 C_{1,3}:=800\max\left\{
 {\mu_1\over\rho_3},
 {{\mu_1/\rho_3+K_gW_3\mu_4}\over\lambda_{patch}},
 K_g\mu_2,{W_3\mu_5\over\rho_3},W_3\mu_3\right\}.     \tag{9}
\]

This is the exact positive-order linear solve; no quadratic leading-repair
term or small-root hypothesis is introduced.  The two `U1` unit-mass bumps
are disjoint and have the proved value bound `B0=320`, so their contribution
is at most `320C1,3`.  The three `E1` bumps and the `E1`-to-`F1` inverse-radius
costs do not enter the axial value norm.  Equations (1)--(9) thus define a
finite selected-parameter global `K1,3*` and complete the first repair.

The formula is deliberately not replaced by a guessed small power of `H`.
Its small denominators and symbolic Picard exponentials remain visible.  A
power-only simplification would require separately comparing every selected
`rho3^-1`, `e_*^-1`, `lambda_patch^-1`, radius and Picard exponent with the
`A,G,H` hierarchy; it is unnecessary for the next source calculation.

## 4. Limits and next obligation

The new bound is sufficient to enter the already derived order-two source
formulas together with the finalized leading `C_eta^4` row.  The next
checkpoint must evaluate the full-support shifted-viscosity totals
`Q10,Q11`, assemble the second five-moment repair, and quantify `H2,3`.

After that remain uniform actual-pulse size/derivative/support estimates at
every physical scale, complete-wave covariance mean corrections including
axial velocity and both exceptional pressure tails, and the complete
viscous/cutoff/interaction residual.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/global_first_coefficient_c3/test_audit.py -v
python3 research/navier_stokes_cascade/results/global_first_coefficient_c3/audit.py \
  --verify-record research/navier_stokes_cascade/results/global_first_coefficient_c3/evidence.json \
  --output global_first_coefficient_c3_fresh.json
```

The output path must be fresh.  Frozen evidence is never overwritten.
