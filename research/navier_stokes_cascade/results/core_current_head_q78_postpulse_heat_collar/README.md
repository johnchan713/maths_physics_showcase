# Q78: the post-pulse leading background through the terminal heat edge

Continue Q77 on [Xend,Xb], Xend=Xp exp(13/lambda), using the
[frozen stage formulas](../outer_pressure_pilot/schedule.py), with the
selected parameters of the [moment audit](../pulse_moment_audit/README.md)
and [heat edit](../heat_exterior_audit/README.md). This bounds the
remaining leading-field jets in the axial test; it is not an estimate
on the actual physical propagator. All endpoints are independent of eta.
U=0 on these stages. This includes the angular corrections: those edit
E only. Keep the exact M=J=0 pulse exit and the inherited nonzero other
moments, including full pressure.

## Interpolation, two relative corrections, release and wait

Use local log-radius y from each stage entrance and E amplitude c equal
to its positive entrance value at eta=0. All c are known scalar products
obtained by evaluating the preceding stage formula at its right endpoint
and eta=0, starting at c=e_b exp((-1/2-lambda)13/lambda). No large
length or high derivative need be determined numerically. Write
B_0=1, B_i=162·3260^i(i!)³ for i>=1 (Q63 step sigma). For a log
field with positive radial derivative bounds H_j on |partial_y^j log E|/j!,
D_0=1 and D_n=sum_(j=1)^n j H_j D_(n-j)/n bound
|partial_y^n E|/n! divided by a cap for |E| at the center. The same
recurrence works for a complex angular argument when the bound is taken
in modulus; it is just Leibniz for E_y=E(log E)_y.

* On the Tf=1000 interpolation, b=-1/2-lambda, s=sigma(y/Tf),
  E=c exp(by) f exp[-s log(2f)], f=(1+eta²)^(-1).
  On eta disks of radius 1/2, |f|<=2, |log(2f)|<=4, using the original
  branch. Thus |E|<=2c exp(4), and
  H_j=(|b| 1_{j=1}+4 B_j/Tf^j)/j! suffices. The normalized mixed
  (n,k) coefficient is <=2c exp(4) D_n 2^k by angular Cauchy.
* On the Tu=30 log(1/lambda) uniform stage let u_i be the SAME two
  angular correction coefficients from Q77 and d_j its normalized
  derivative bounds. The field is c exp(by)(1+sum u_i beta_i(y)).
  The two relative bumps have width .3, centers Tu-3 and Tu-1 and
  NO unit-mass prefactor. Put V_0=64 and V_j=B_(j+1)/.3^j for j>=1.
  The raw y derivative of order n in Q77's angular algebra of order k
  is bounded by
  c[|b|^n + 2(sum_(j=0)^k d_j) sum_i binom(n,i)|b|^(n-i)V_i].
  This includes derivatives of the product. The angularly uniform base
  is restored exactly after the bumps, not the intermediate pressure.
* On the two unit release ramps, the H=sqrt(2X)E slope changes from
  -lambda to -1, then from -1 to -h, using sigma(y). The E log slope
  is the interpolated slope MINUS 1/2. Use H_1=3/2 and
  H_j=B_(j-1)/j! for j>=2 in the recurrence above; |E|<=c.
  On the intervening hold of length 4 log(1/h) and the final wait of
  length log(Q_before/Qp)/(1-h), the constant E slopes are respectively
  -3/2 and -1/2-h. Hence raw radial order n costs at most c(3/2)^n.
  These stages are eta-independent, and the retained exact waiting
  length is a scalar, so differentiation never moves a seam.

## Two terminal cutoffs and a real finite heat table

On 0<=y=log(X/Xtail)<=3, Xb=Xtail exp(3), write
 rho=co h, b=-1/2-h and

    E=c exp(by) [1-rho sigma((3-y)/2)]/(1-rho)
                [1+chi(y)(Hcal(Z)-1)],
    chi=sigma((y-.2)/.3),   Z=2(1-eta²)/X.

The first bracket divided by its entrance value is not chi. The
infinite tail has the continued power and Hcal, with both cutoff
factors already constant. For real Z>=0 the Q76 Gamma bound gives
|Hcal^(a)(Z)|/a!<=h(a+1)! for a>=1 and |Hcal-1|<=4h/X.
It does NOT give a complex disk at Z=0. For desired orders r,k use the
rectangle through (r,k) of finite normalized jets in y,eta; all positive
series below are truncated there. With formal increments p,q define

    J=2[(exp(p)-1)+exp(p)(2q+q²)],
    W=4+sum_(a=1)^(r+k) (a+1)! J^a.                         (1)

Indeed the increment of Z is
(2/X)[(1-eta²)(exp(-p)-1)-exp(-p)(2 eta q+q²)].
Thus J/X dominates it coefficientwise. Since X>=1, composition gives
normalized mixed coefficients of Hcal-1 bounded by (h/X)W.
The constant entry uses the mean value bound; the positive-degree
entries follow by finite differentiation on the interior and their
one-sided limits. Polynomial factors in the Gamma variable are
integrable at every fixed order, so those limits exist. The proof
keeps the angular quadratic increment and the radial derivatives hitting
it. In particular it does not extend the integral to Z<0.

For step derivatives set raw-jet polynomials
S_w(p)=sum_(i=0)^r B_i p^i/(i! w^i). A single coefficientwise
bound for the terminal normalized E table is

    c exp(|b|p) (1+rho S_2(p))/(1-rho)
                      [1+ S_.3(p) (h/X) W(p,q)].             (2)

Use X>=Xtail in the denominator for a uniform bound. All smooth joins
are at their original radii and are covered by the one-sided tables.
Terminal U is zero. This is an all-order formula at each FINITE pair,
not an assertion of radial or angular Taylor convergence.

## Consequence for the physical test

Use Q70's Stirling conversion on each positive interval (multiply the
log coefficient by its radial factorial before applying the raw
operator). Take maxima with Q69--Q77, including Q74's actual modulation
and accepted repair, to bound the mixed leading jets on the whole
native stress shell [Xa,Xb]. On any additional part of a fixed positive
chart hull beyond Xb, (1),(2) apply with constant cutoffs and the
decaying power; below Xa use the actual analytic core on the positive
image. One can therefore form Q58's leading weight constants and the
leading-profile derivative inputs in Q60 and Q62 at their requested
finite orders without an unknown post-pulse velocity seminorm. In
particular on this new interval G=0, hence w3=0; w4 retains both the
radial derivative of F and the moving pressure-cutoff moment. Nothing
sets that multiplier to zero just because the axial field vanishes.

For the unstripped target keep the five Q62 terminal identities and
integrate the true densities, or propagate the matched pulse exit
traces forward. In particular Pi_X=E²/(2X); chi and the reference
cutoff both enter this density. Exact cancellation of total moments
does not cancel Pi_eta locally. No estimate of a stress divided by
its vanishing edge factor follows from an unstripped derivative bound.

Next verify the actual lower-order full-background errors, normal and
column inverse floors, and the weighted two-time propagator bounds on
the SAME physical hull. Those are independent inputs of Q61; bounded
leading derivatives do not certify the axial rate, the distinct radial
pressure tail, all-scale physical waves, complete residual improvement
or blow-up. Run test_audit.py and audit.py --verify-record evidence.json
--output a fresh path. Manufactured algebra tests are not continuum proofs.
