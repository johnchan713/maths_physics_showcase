# Q71: the continuation driven by the radially flattened reference

Let delta=C^(-20), Xa=4/Lambda, B=Xa exp(2 delta) and Xi=110.
This gives quantitative bounds on the ACTUAL prescribed continuation on
[B,Xi] for arbitrary finite radial order, through six angular field
derivatives and five angular state derivatives. Q70 then applies from Xi
onward. We use the inherited exact core, native sixth-angular estimates
and earlier annular match. Nothing here bounds the remaining short interval
immediately to the left of B by extrapolating the flattened reference.

## Three radial powers, with all incoming constants retained

After the reference taper, F_r=E_r/sqrt(2X)=F(eta)>0 and U_r=u(eta).
The five **reference** moments on [B,Xi] are their actual reference values
at B plus the following increments, in the order M,I,J,S,Pi:

    u(X-B),  F(X²-B²),  Fu(X²-B²),
    u²(X-B)-F²(X²-B²)/2,  F²(X-B).                        (1)

These are the exact integrals of the five X densities; the full Pi(B)
contains the original axis datum. In particular flat F_r does not mean
constant E_r, zero pressure, or constant integrated source coordinates.
From the [complete reference map](../prescribed_continuation_foundation/README.md),
W=1-(2D eta M+d M_eta)/X, H=2XF, and
p1r=(-XW+N_Q/H)/L with N_Q linear in I,J and their eta derivatives.
Thus W has powers 0,-1 and N_Q has powers 0,2. Similarly

    nsr=[-Wu+(D(M-eta M_eta)+4h eta S-d S_eta)/X
                  +4A eta Pi-d Pi_eta]/L,

where d=1-eta², L=1-2h eta² and A,D=1/2±h respectively.
Consequently BOTH p1r and nsr belong to span{X^(-1),1,X} with
eta-dependent coefficients. This assertion is about nsr=N_r/L, NOT
p2r=X nsr/E_r. In particular it neither removes a row of memory nor
identifies p1r with the vanishing physical shear of the reference.

## Bound the modes without differentiating unknown trace coefficients

In the factorial C6 angular algebra of Q70, the [native sixth-order
source bounds](../native_sixth_angular/README.md), equation (3), give
Vp=10^8 Xi Lambda K0^26 and Vn=10^5 K0³ as uniform bounds for p1r and
nsr. For a three-power source v with such bound V, at each X write

    v(q^t X)=q^(-t) a_- + a_0 + q^t a_+,    t=0,1,2.      (2)

Choose q=2 if 4X<=Xi, else q=1/2. All three points lie in [B,Xi]:
Xi/B>=16 since Xa<1, 2delta<1/50. The row sums of absolute values of
the inverse of (2) for q=2 are 8,9,2; for q=1/2 interchange 8 and 2.
For example the q=2 inverse rows are (8/3,-4,4/3), (-2,5,-2),
(1/3,-1,2/3). This is interpolation of an exact three-dimensional
function space against an analytical supremum, not numerical sampling.
One never differentiates the choice of q. At the center, y=log(X/Xa),
positive derivatives of v have values (-1)^r a_- + a_+.
For Xv the powers instead are 0,1,2. Hence for all r>=0 sufficient
raw derivative bounds in that angular algebra are

    ||partial_y^r p1r|| <= 10 Vp,
    ||partial_y^r (X nsr)|| <= 19 Xi 2^r Vn.               (3)

At r=0 the original value bounds suffice. The constants do not carry
any power of the tiny taper width, since (3) concerns the already flat
reference region alone.

## Return to the actual prescribed slopes and memory

On this interval kappa=delta. Extend beta as one up to X=100, switch
it to zero over the original log width delta, then switch theta from
zero to one on the following width delta. Extend their plateaus to
the whole interval. The exact prescriptions are

    a=delta(1-theta)p1r+(4/5)theta,    (log E)_y=(1-a)/2,
    U_y=-(delta/2) beta X nsr.                               (4)

The last plateau has slope 1/10 for log E. These are not equations
replacing the actual fields by reference fields. Put T_0=1,
T_i=162·3260^i(i!)³ delta^(-i) for i>=1 from the actual Q63 step.
For r>=0 set

    A_r=delta sum_i binom(r,i) T_i 10Vp + (4/5)T_r,
    G_r=(delta/2) sum_i binom(r,i) T_i 19 Xi 2^(r-i) Vn,

summing i=0,...,r. Leibniz gives ||partial_y^r a||<=A_r,
||partial_y^(r+1) U||<=G_r. For j>=1 the normalized log-field
bound is H_j=(1_{j=1}+A_(j-1))/(2j!). Use Q70's exponential
recurrence D_0=1, D_j=sum_(k=1)^j k H_k D_(j-k)/j. Since
||log(CE)||_6<=J6 and ||U||_6<9, this gives explicitly

    ||partial_y^j E||_6 <= exp(J6) j! D_j/C,
    E>=exp(-J6)/C,  ||U||_6<9.                              (5)

The first bound at zero follows by exponentiation in the normed algebra;
for positive order differentiate E_y=E(log E)_y. Q70's falling-factorial
conversion (3) with X0=B now supplies every normalized mixed physical-X
coefficient with angular order <=6. The powers delta^(-i) in the
remaining switches cannot be set to one.

For the ACTUAL state take its five traces at Xi from Q70, using the exact
earlier annular match there via backwards integration from X_R. Subtract
the integral from X to Xi of the ACTUAL densities of Q67, formed from
(5) and the U bounds. This gives pure angular primitive rows <= the Xi
trace bound plus (Xi-B) times the density bound, and normalized positive
radial rows <= the preceding radial density coefficient divided by the
radial order. Q67/Q69 then bound p1,p2,a,b_s at arbitrary radial order
and angular order <=5 with both pressure terms. Formula (1) must NOT
be used for this last integral. These bounds need not meet a physical
smallness threshold.

Next quantify the native and taper inputs on the short interval ending
at B, and the higher angular orders of the shape and annular data when
needed. The normalized compact repair defect, uniform physical pulses,
pressure tails, full residual improvement and blow-up remain unproved.

Run test_audit.py and audit.py --verify-record evidence.json --output a fresh
path. The symbolic controls check the reduction and its indexing, not
the inherited continuum existence or smallness hypotheses.
