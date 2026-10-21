# Independent same-profile modulation, repair, and regional cone recheck

Parent: `f9c45658e75a6f0a8122875e383185679f4bec55`.
Status: `same-profile-compact-modulation-repair-and-regional-cone-join-independently-rechecked`.

The parent fixes one actual continued profile, identifies both matching roots
inside the old uniqueness balls, and carries all five forward moments and
`Pi_axis+Cp` to the compact input.  This checkpoint performs the next
operation on that same profile: it rechecks the signed variance root, the
moving phase, the selected finite-frequency modulation, and application of
the same exact compact repair.  It then joins the resulting value cone to the
unchanged inner and annular/exterior regions.

No parameter, support, pressure gauge, root, frequency, or activation width is
changed.  The estimates below are analytical implications of the displayed
continuum bounds.  Exact-arithmetic tests check their algebra and scale
comparisons; they are not a proof-assistant certificate or a numerical
evaluation of the enormous selected field.  Global heat-tail pressure/moment
identities, positive coefficients, physical pulses, covariance means, the
complete residual, and blow-up are not claimed here.

## 1. The variance root is signed and regular through zero

For

\[
 q(z,\theta)=\frac{e^{z\sin\theta}}{I_0(z)},\qquad
 R(z,\theta)=\int_0^1q_z(tz,\theta)\,dt,
\]

write `z=mu p` and

\[
 F(\mu,p)=d_0\mu\sqrt{\langle R(\mu p,\theta)^2\rangle}.
                                                               \tag{1}
\]

This is the signed root, not the nonnegative square root on both sides of
zero.  At `z=0`, `R=sin(theta)`, so

\[
 \langle R^2\rangle=\tfrac12,qquad
 F_\mu(0,p)=d_0/\sqrt2,\qquad F_p(0,p)=0.                 \tag{2}
\]

Thus the ordinary implicit-function denominator stays nonzero at zero
variance.  Dividing by `V_mu` there would be singular and is not used.  The
explicit source-admissible flat cutoff gives a smooth right side

\[
 s=\zeta(v_s)\sqrt{(2+\delta_L/2-v_s)/a},                 \tag{3}
\]

extended by zero across its flat boundary.  The retained lower bound
`F_mu >= (d0/32) exp(-4|mu p|)` away from zero and (2) give one smooth root
on a real neighborhood of the closed angular interval.  This is the same
root family used in the accepted compact loop, now applied to the parent's
same continued input.

## 2. The auxiliary phase must move with the slow variables

With

\[
 w=\frac{a(1+t^2)}{2\pi v},\qquad
 \phi(x,\theta)=\int_0^\theta w(x,r)\,dr,
 \qquad \phi(x,\Theta(x,\varphi))=\varphi,               \tag{4}
\]

the variance equation gives one normalized period exactly and the retained
floor gives `w>0`.  Hence

\[
 \Theta_\varphi=w^{-1},\qquad
 \Theta_x=-\phi_x/w,
\]

and every pulled-back slow derivative contains

\[
 \partial_x[t(x,\Theta(x,\varphi))]
   =t_x+t_\theta\Theta_x.                                \tag{5}
\]

The second term is generally nonzero.  Freezing `Theta` would therefore
change the loop and its repair target.  The exact rational control in
`proof_checks.py` uses a normalized nonconstant phase map and rejects that
error.

Let `Acal,Bcal` be the zero-mean loop primitives and substitute
`varphi=Ny`.  Differentiating before substitution gives

\[
 \partial_y\frac{\mathcal A(y,\eta,Ny)}N
   =\mathcal A_\varphi+\frac{\mathcal A_y}N.              \tag{6}
\]

Thus a radial field derivative is not `O(N^-1)`.  The fast term is the
designed loop; the error relative to that loop is the slow remainder.  For
every fixed pair `(r,k)`, the exact chain rule yields

\[
 |\partial_y^r\partial_\eta^k\Delta(E,U)|
       \le C_{rk}N^{r-1},                                \tag{7}
\]

while integration of the five actual densities gives the extra radial gain

\[
 |\partial_y^r\partial_\eta^k
   \Delta(M,I,J,S,C_p,\Pi)|
       \le C_{rk}N^{\max(r-2,-1)}.                       \tag{8}
\]

Equation (8) keeps the incoming moments and `Pi=Pi_axis+Cp`; it uses no
oscillatory cancellation and no pressure reset.

## 3. The selected integer accepts the same exact repair

Keep

\[
 A=C^{4096}\ge2^{16},\quad G=e^{A^{32}},\quad H=e^{G^4},
 \quad N=1+\lfloor H^{32}\rfloor,\quad D=H^{16}.         \tag{9}
\]

The parent's compact map is

\[
 Bc+Q(c,c)=z_N,qquad
 T_c(M,I,J,S,C_p)=(M,(M-J)/\lambda,I,-S,C_p),            \tag{10}
\]

where `lambda>0`, `||B^-1||<=800`, `||Q||<=25000`, and
`||z_N||_C1<=D/N`.  Its exact inverse is

\[
 (r_1,r_2,r_3,r_4,r_5)\mapsto
 (r_1,r_3,r_1-\lambda r_2,-r_4,r_5).                    \tag{11}
\]

At `lambda=0`, (10) cannot encode independent `M,J` equations, so the
positive selected value is essential.  The fifth `Cp` row is equally
essential.

The contraction coefficient satisfies

\[
 8(800)^2(25000)=128000000000<2^{37},
\]

and, since `D/N<H^-16` and `H>2^16`,

\[
 8\beta^2qD/N<2^{37}H^{-16}<2^{-219}.                   \tag{12}
\]

The same root therefore has `C1` norm at most `1600H^16/N`.  It lies in the
already accepted uniqueness ball and in the repair neighborhood
`A^-128/2`.  For the latter it suffices that

\[
 3200A^{128}<H^{16}.                                     \tag{13}
\]

Indeed `log(3200)<12`, `log A<A`, and `G^4>5A^32`, so the logarithm of the
left side is below `12+128A<80A^32<16G^4`.

The loop state error is bounded by `H^16/N`; the repair-to-state constant is
`A^128`.  Hence the complete error relative to the intended loop state is

\[
 e_{state}\le\frac{H^{16}+1600A^{128}H^{16}}N.           \tag{14}
\]

To get `e_state<G^-64/2`, it suffices that

\[
 2H^{16}(1+1600A^{128})G^{64}<H^{32}.                    \tag{15}
\]

After logarithms, the non-`16G^4` part is below
`13+128A+64A^32<80A^32<16G^4`.  This proves (15) without
materializing the hierarchy.  No new frequency or high-norm contraction is
introduced.

## 4. Compact cone loss and the regional join

The exact loop has common gap `kappa=G^-4`, first-coordinate floor `G^-3`,
and lies in a box of radius `R<=4G^3`.  The four cone functions have common
Lipschitz bound `256R^10`.  With `epsilon=G^-64` and (14)--(15),

\[
 256R^{10}e_{state}
 <128\,4^{10}G^{-34}<\tfrac12G^{-4},                    \tag{16}
\]

because `G^30>=2^480>2^29=512*4^10`.  Also
`epsilon<G^-3`, so the first coordinate stays at least `G^-3/2`.  The
repaired compact gaps are therefore at least `G^-4/2`.

The parent's support ledger now identifies the whole chain on one profile:

* the modulation and repair do not touch the protected inner collar;
* the annular and compact coefficients are the same unique roots already
  checked on that profile;
* all five compact endpoint defects vanish as functions of `eta` at
  `Xrep=rho_p exp(5)=Xstar`;
* the axis pressure datum is unchanged, and the later heat edit starts at
  `Xstar` with the already audited annular/exterior value estimates.

Consequently the inner value cone, (16), and the annular/exterior value cone
cover the active leading profile with their shared endpoints and no root or
frequency reselection.  The activation width remains `C^-20`.  This is a
leading value/direction join; it is not a uniform physical-pulse derivative
or weighted positive-coefficient theorem.

The pressure row cannot be inferred from the other four.  An omitted defect
`p(eta)` contributes

\[
 \frac{X\{4(1/2+h)\eta p-(1-\eta^2)p_\eta\}}
        {1-2h\eta^2}                                    \tag{17}
\]

to the axial stress numerator.  With `p(0)=0,p_eta(0)=1`, (17) is `-X`.
The exact fifth-row endpoint equality avoids this defect but does not, by
itself, reprove the later infinite heat-tail identities.

## 5. Exact frontier

This checkpoint independently closes the fixed-profile compact modulation,
same exact repair, and leading regional cone obligation.  The next checkpoint
must retrace the global renormalized heat moment and complete pressure/moment
identities on this repaired profile.  After that remain finite coefficient
supports and reserved means, uniform physical pulses at every scale,
complete covariance means including axial velocity and both exceptional
pressure tails, and the viscosity/cutoff/interaction residual.  No
unconditional pulse bound, residual improvement, force, or blow-up statement
is made.

Reproduce from the repository root:

```bash
python3 research/navier_stokes_cascade/results/core_modulation_repair_join_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_modulation_repair_join_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_modulation_repair_join_recheck/evidence.json \
  --output core_modulation_repair_join_recheck_replay.json
```
