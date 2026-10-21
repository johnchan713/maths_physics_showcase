# Coefficient support join and an unavoidable reserved-pressure error

Parent: `3651a61c308633530eafd9ceb728b3355177d3b8`.
Status: `same-profile-finite-coefficient-supports-and-reserved-pressure-obstruction-rechecked`.

The new global moment checkpoint supplies the same leading profile's five
functional identities. This checkpoint inserts them into the existing finite
coefficient construction and rechecks its ordinary row normalization, shifted
viscosity totals, support induction, and reserved-patch curl. Its new estimate
is a sharp pressure-error floor on the reserved mean interval. It also derives
the complete time and axial derivatives of both surviving pressure orders.

The finite inner Volterra existence and regularity lemmas remain analytical
dependencies of this join. Exact rational controls check the algebra; they
neither construct the selected enormous coefficient fields nor prove a
continuum supremum. All historical evidence, parameters, roots, supports,
frequency, pressure gauge, and width `C^-20` are preserved.

## 1. Same-profile inputs and the positive-order linear solve

Use `A=1/2+h`, `D=1/2-h`, `d=1-eta^2`, `L=1-2h eta^2`, and
`lambda_n=2nh`. The preceding checkpoint gives `M_inf=R=J_inf=S_inf=0`
and `Pi_axis+Cp_inf=0` on this leading profile. The global angular row is
renormalized at leading order; positive orders have ordinary compact moments.

On I3 put `X=rho3 exp(y)`, `0<y<5`. The same reserved leading plateau has

\[
 U_0=0,\quad E_0=e(\eta)e^{-(1/2+\lambda)y},\quad
 e(\eta)=e_* /(1+\eta^2)>0,\quad M_0=\mu\eta,\quad V_0=-\mu.       \tag{1}
\]

Here `mu>0`; zero axial velocity does not remove accumulated radial flux.
Finalize orders below `n` before forming the ordinary discrepancy vector
`(m1,m2,m3,m4,m5)`:

\[
 \left(\int U_n,\ \int\sqrt{2X}E_n,\ \int(\Pi_n)_X,
 \ \int\sqrt{2X}\sum_{i+j=n}U_iE_j,
 \ \int\{\sum_{i+j=n}U_iU_j-X(\Pi_n)_X\}\right).          \tag{2}
\]

Every integral includes its whole source support. At degree `n>=1`, a
current-order field occurs only in the pairs `(0,n),(n,0)`. Its square
belongs to degree `2n`. Since `U0=0` on I3, the changes from the five
additive bumps are exactly

\[
 \left(\int\delta U_n,\ \int\sqrt{2X}\delta E_n,
 \ \int E_0\delta E_n/X,\ \int\sqrt{2X}E_0\delta U_n,
 \ -\int E_0\delta E_n\right).                           \tag{3}
\]

Let `K=1/(sqrt(2)rho3^(3/2))`. Divide (2) by
`(rho3,1/K,e,e/K,rho3 e)` to obtain `(M,I,Cp,J,S)`, then use rows
`(M,(M-J)/lambda,I,-S,Cp)`. The correction is `c_n=-B^-1 z_n`, with the
same constant bump matrix and continuum inverse bound `<800`. The rational
positive-mass determinant argument in `matching_root_foundation` is replayed
directly; no sampled matrix inverse replaces that argument.

For the factorial angular norm, the reciprocal amplitude costs are

\[
 W_k=|e^{-1}|_k=e_*^{-1}\begin{cases}2&k=0,\\4&k=1,\\5&k\ge2.\end{cases}
\]

If `mu_i=|m_i|_k`, the bound is

\[
 |c_n|_k\le800\max\left\{
 \mu_1/\rho_3,\ (\mu_1/\rho_3+KW_k\mu_4)/\lambda,
 \ K\mu_2,\ W_k\mu_5/\rho_3,\ W_k\mu_3\right\}.          \tag{4}
\]

This keeps both the amplitude derivatives and the raw `lambda^-1` cost.
It is finite for every fixed order and selected parameter. It supplies no
uniform bound as the coefficient order, derivative order, or physical scale
varies. There is no small-target condition for this linear solve.

## 2. Finite support induction and its pressure exceptions

The inherited six-variable inner Volterra system constructs every finite
coefficient on the same radial collar from finalized lower sources. Its
analytic angular strip may shrink with order. Apply the fixed inner cutoff
to `F_n,U_n`, add the I3 bumps, and reconstruct `M_n,V_n,Pi_n` from zero
axis data. The bumps end by `y=2.55`; `Xvel=rho3 exp(3)` lies before I4,
whose range is `6<log(X/rho3)<11`.

The first exact total gives `M_n=0` beyond Xvel, including every angular
jet. The flux formula

\[
 V_n=\{2\eta XU_n-2(D+2nh)\eta M_n-d(M_n)_\eta\}/L       \tag{5}
\]

then gives zero positive-order velocity there. Multiplying the potential
by `chi(c_n q)` adds `-2eta(c_n q)chi'(c_n q)M_n/L` to the radial flux;
it also vanishes on I4. Thus any locally finite potential-cutoff sum agrees
with leading velocity on that interval, even inside cutoff transitions.

Pressure is reconstructed by

\[
 (\Pi_n)_X=\sum_{i+j=n}F_iF_j-\Omega_{n-1}/(2X).         \tag{6}
\]

After Xvel, its remaining sources are `-Omega0/(2X)` in order one and
`nu Z_-D Z_0 V0/(2X)` in order two. Higher orders have zero source.
The third exact total removes the final pressure constants. For stresses,
the conservative total identities also include the previous-order axial
viscosity: `-nu Z_(a_(n-1)-D)Z_(a_(n-1)) I_(n-1)` and
`-nu Z_(g_(n-1)-D)Z_(g_(n-1)) M_(n-1)`, where
`a_j=1-h+2jh`, `g_j=1-A+2jh`. The new leading renormalized `R=0` supplies
the `n=1` angular term. Current exact totals and induction remove the other
terms. The exact-rational residual integration controls keep this viscosity
active on manufactured profiles with nonzero lower moments.

The inherited stress-free inner collar and backward primitives therefore
give these right support envelopes on the same profile:

| Order | Velocity | Pressure | Stress |
|---|---|---|---|
| 1 | Xvel | Xv | Xb |
| 2 | Xvel | Xv | Xv |
| n>=3 | Xvel | Xvel | Xvel |

On I4, `V0=-mu`, `Omega0=-mu^2/(2X)` and `Z_-D Z_0 V0=0`, so

\[
 \Pi_1=K_1(\eta)-\mu^2/(4X),\quad
 \Pi_2=K_2(\eta),\quad \Pi_n=0\ (n\ge3).                \tag{7}
\]

K1 and K2 contain the later pressure source memory. A local zero derivative
does not determine either constant. The inherited manufactured smooth-drop
control has a nonzero order-two plateau `3nu mu` at `eta=0`; this is a
counterexample to deleting K2, not its value for the actual profile.

## 3. Sharp necessary inequality for pressure accuracy on I4

At fixed `(q,eta)`, define `P=q^(2A)p_B-Pi0`. Equation (7) gives

\[
 P(X)=\chi_1q^{2h}[K_1-\mu^2/(4X)]
             +\chi_2q^{4h}K_2,\qquad \chi_i=\chi(c_iq).             \tag{8}
\]

For any `Xa<Xb` in I4, both unknown angular constants cancel exactly:

\[
 P(X_b)-P(X_a)=\chi_1q^{2h}\frac{\mu^2}{4}
                    (X_a^{-1}-X_b^{-1}).                         \tag{9}
\]

The triangle inequality yields the sharp two-point lower bound

\[
 \sup_{X\in I4}|P(X)|\ge
 \frac{|\chi(c_1q)|q^{2h}\mu^2}{8\rho_3}
                    (e^{-6}-e^{-11}).                            \tag{10}
\]

Take the limit to the open interval's endpoints to obtain (10). The factor
`1/8` is optimal when the X-independent part of (8) is arbitrary: centering
its two endpoint values attains the bound. This does not assert the actual
K functions choose that center. The positive exponential difference is
verified with rational Taylor remainder enclosures.

Thus any proposed tolerance `epsilon(q,eta)` for replacing the full reserved
pressure by leading pressure must satisfy the right side of (10)
`<=epsilon(q,eta)`. It is strictly positive when `chi(c1q)!=0`, `q>0` and
`mu>0`. The construction's genuine gauge remains fixed; the impossibility
even allows arbitrary X-independent cancellation.

## 4. Complete derivative transfer and an exact stationary cancellation

For `q^b chi(cq)P(X,eta)`, the profile coordinate rules give

\[
 \partial_z=\frac{q^{b-D}}L\left[
 \chi\{2b\eta P-2\eta XP_X+dP_\eta\}
                       +2\eta\sigma\chi'(\sigma)P\right],       \tag{11}
\]
\[
 \partial_t=\frac{q^{b-1}}L\left[
 \chi\{-bP+XP_X+D\eta P_\eta\}
                            -\sigma\chi'(\sigma)P\right],       \tag{12}
\]

where `sigma=cq`. Use `b1=-2A+2h=-1` and `b2=-2A+4h=-1+2h` for the
two physical pressure orders. For `P=-mu^2/(4X)` in order one, both brace
terms vanish exactly. Indeed

\[
 q^{-1}[-\mu^2/(4X)]=-\mu^2/(2r^2).                       \tag{13}
\]

This piece is stationary in physical time and axial position while chi is
constant. Its radial derivative remains nonzero. During a cutoff transition,
the last terms of (11)--(12) restore time and, away from `eta=0`, axial
contributions. The K1 and K2 terms also have their own derivatives in
(11)--(12). A full mean/covariance calculation needs all of these terms;
zero positive-order velocity alone eliminates none of (8).

## 5. Reproduction and next obligation

Run `python3 test_audit.py -v` and
`python3 audit.py --verify-record evidence.json --output fresh.json` here.
Only the standard library is required for this checkpoint. The record binds
each dependency's declared source hashes and checks the inherited exact
conservative residual controls. Separate reproduction of the older quadrature
suites requires mpmath and their declared numerical dependencies. Those
diagnostics remain separate from continuum existence and support arguments.

Local validation on 2026-10-21: 191 tests passed across this checkpoint (24),
the global pressure/moment recheck (17), modulation/repair recheck (16),
coefficient foundation (31), historical support audit (34), annular/exterior
foundation (34), and weighted-stress foundation (35). Their frozen records
were verified in those suites; the historical support CLI also replayed to a
fresh output. Numerical suites used mpmath 1.3.0 in an isolated dependency
directory. These checks do not upgrade the retained analytical premises.

Next: discharge the remaining actual-pulse existence, inverse, gap, value,
derivative and support bounds uniformly over every physical scale, with the
selected cutoff and retained tilt deficit. Then join full wave covariance
and means using (8), (11)--(12), radial and axial velocity and curl remainders.
The complete viscosity/cutoff/interaction residual and blow-up remain open.
