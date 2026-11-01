# Q77: the exterior axial profile pulse and its actual strength

This supplies mixed leading-field bounds on Xp<=X<=Xp exp(13/lambda).
This pulse is part of the slow profile entering the physical wave generator;
it is not the canonical high-frequency physical wave. Keep the exact
Md64 roots, supports and all five true incoming moments.

## Which energy selects the root?

Use the equations and constants of [the pulse moment audit](../pulse_moment_audit/README.md).
The first compact repair agrees in all moments after its exit. The subsequent
heat compensation changes E only on a U=0 stage before Xp. It leaves M,J
unchanged and the value Ep=e_b f unchanged at Xp, but generally changes
S(Xp). The terminal heat edit is after the axial stage and also has U=0.
The exact heat compensation equation for S says their two changes cancel
in S(infinity), identically in eta. Both edits are independent of the trial
axial strength. Thus the PRE-HEAT energy equation below selects the same
strength even after both edits. Substituting post-compensation s0 but a
pre-heat future mass would be wrong. The actual partial pressure and
energy on the pulse must instead start from Q76's actual traces.

## Quantitative angular bound for the post-pulse mass in that equation

Use ||v||_k=sum_(j<=k) sup|partial_eta^j v|/j!. The pulse audit's
two later angular (A.11) coefficients, denoted u, satisfy a constant
quadratic map with inverse linear norm <=5, bilinear norm <=32 and
accepted ||u||_infty<=7000 lambda^30. Its normalized target is
(-exp(-beta(Tu-3)) delta_I,0), beta=1-lambda. In the formula in its
section 3 the only angular scalar is a=log(2/(1+eta²)). On complex disks
of radius 1/2 at real centers, Re(1+eta²)>=3/4, |1+eta²|<=13/4.
Thus |a|<=4 (principal log), |exp(a v)|<=exp(4) for 0<=v<=1.
Using |delta0|<=1/beta, integral sigma'=1 and beta>0 gives the explicit
disk cap Z=5 exp(4) exp(-beta(Tu-3))/beta for the target. Cauchy
bounds its normalized kth derivative by 2^k Z. Only the target is
extended here; the radial cutoff is evaluated on its original real argument.

Take d0=7000 lambda^30 and
d_j=10(2^j Z+32 sum_(i=1)^(j-1) d_i d_(j-i)). Since
2·5·32 d0<1/2, the pointwise Jacobian inverse is <=10. Differentiating
the original equation shows sup|u^(j)|/j!<=d_j. The high derivatives
are not used to select a new small root.

Let R be the PRE-HEAT future integral of E² divided by Xend Eend².
Remove the two angular bumps temporarily, keeping the original schedule
and release/wait lengths, and use the unedited normalized density at
eta=0 as a positive integration measure. Its mass is <=512T: at eta=0
the restored squared multiplier is >=.99² everywhere and the original
audit bounds R(0)<=256T. Relative to this measure the unedited density
at eta is exactly (1+eta²)^(2(1-theta)), where theta decreases from
one to zero on interpolation and is zero thereafter. This follows by
dividing the swirl by its pulse-end factor f, including the normalization
on every later uniform stage. The scalar has disk modulus <=16 and
therefore angular norm <=16 sum_(j=0)^k 2^j. Restoring the two bumps
multiplies by (1+sum u_i beta_i)^2, with |beta_i|<=64. Hence

    ||R||_k <= 512T · 16(sum 2^j) (1+128 sum_(j=0)^k d_j)^2.     (1)

The sum in the first factor also runs from zero to k. Fixed endpoints
and this integrable majorant justify differentiation of the entire mass.
Both bump squares and the infinite release tail are retained.

## Angular jets of the axial strength and its end corrections

Write c_i(A,eta)=v_i A+b_i p(eta), p=eta(1+eta²). This is exact:
the two incoming normalized M,J are constant multiples of p, the main
pulse drives an eta-independent RHS linear in A, and the end-bump
matrix is constant with inverse <20/lambda at the retained positive
lambda. The pulse audit gives |v_i|<=lambda^20 and
|c_i(.9,eta)|<=lambda^20. Evaluation at eta=1, where p=2, gives
|b_i|<=lambda^20. No cancellation of the small determinant is presumed.
The pre-heat s0=s_*+s_1 p², |s_*|<=500T, |s_1|<=1, by the prefix
U² and E² decomposition (U² contributes <=1 on the whole interval).
All these polynomial factors have directly bounded normalized jets.

Use tau_i=int exp(-2lambda y) beta_i² dy<=20. The EXACT energy polynomial
(after the cancellation above) is

    F(A)=K_b A²-(1-exp(-26))/4 + lambda s0
          -(lambda/2) exp(-26) R + lambda sum tau_i (v_i A+b_i p)².

K_b<.246; the original real root A0 lies in [.9,1.2] and F_A>.35.
Here is a finite bound procedure not involving high jets of A as inputs.
At each center form P(x)=(1+x)(1+(1+x)^2), which dominates the absolute
Taylor coefficients of p. Use (1) for the normalized jth coefficient of
R (use k=j). Set Abar_0=1.2. For j>=1 form V=sum_(i<j) Abar_i x^i
and set Abar_j to 1/.35 times the degree-j coefficient in

    .246 V² + lambda (500T+P(x)²)
     + (lambda/2) exp(-26) sum_(i=0)^j bound(R,i) x^i
     + 40 lambda · lambda^40 (V+P(x))².                         (2)

Indeed the missing highest coefficient occurs in F only as F_A A_j,
and the real value inverse is <=1/.35. Induction bounds
|partial_eta^j A|/j! by Abar_j. The corresponding coefficient of either
actual c_i is <=lambda^20(Abar_j+[x^j]P). This uses the full trial-strength
dependence of the end corrections, even though their energy is tiny.

## Mixed radial derivatives, pressure and transfer

In y=log(X/Xp), E=e_b f exp(alpha y), alpha=-1/2-lambda,
U=E[A R0(lambda y)+sum c_i beta_i(y)]. Use the original
R0(xi)=phi(xi)(1-sigma(xi-10)), phi'=sigma(xi/.02).
On 0<=xi<=13, phi<=13; for n>=1 its derivative is bounded by
B_(n-1) .02^(-(n-1)), B_0=1, B_i=162·3260^i(i!)³ for i>=1.
The complementary factor has raw bounds B_i. Thus R0(lambda y) at
radial order r costs L_r=lambda^r sum_i binom(r,i) Phi_i B_(r-i)
with Phi_0=13 and Phi_i as above. The translated beta_i=sigma' at
width .3 instead cost H_r=B_(r+1) .3^(-r): these end terms do NOT
receive lambda^r. The flat extension at the two ends agrees with U=0.

For angular order k take A_k^*=sum_(j<=k) Abar_j,
C_k^*=lambda^20 sum_(j<=k)(Abar_j+[x^j]P) and f_k=2 sum_(j<=k)2^j.
A uniform raw radial bound in the angular norm is then

    ||partial_y^r U||_k <= e_b f_k sum_(i=0)^r binom(r,i)
                     |alpha|^(r-i) (A_k^* L_i+2 C_k^* H_i).

For E the corresponding bound is e_b f_k |alpha|^r. Q70 (3) with X0=Xp
converts to normalized physical X jets. These provide explicit finite
leading-field inputs to Q58/Q60 on the pulse region. Use the genuine
positive E floor e_b exp(alpha·13/lambda)/2 if a state reciprocal is
needed. Integrate the five actual Q67 densities through the interval;
M,J vanish at the right endpoint by the exact affine solve, not at
interior cuts. Neither the local pressure Pi nor its angular derivative
is obtained by making the unedited s0 substitution used only for total
energy. The construction does not change any reserved support.

Next bound the remaining post-pulse interpolation, uniform release/wait
and terminal heat collar through Xb, then verify the physical inverse and
two-time propagator gates. The pressure tails, all-band physical pulse,
complete residual and blow-up remain unproved. Run test_audit.py and
audit.py --verify-record evidence.json --output a fresh path.
