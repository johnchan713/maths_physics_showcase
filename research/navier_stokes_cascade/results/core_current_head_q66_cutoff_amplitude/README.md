# Q66: explicit all-order cutoff-amplitude envelope

Parent: `828450b93865cd97f6b9e7314432e6f6d460f801` (Q65).
This checkpoint removes the independent signed-amplitude jet input from
Q65, conditional on the original pre-modulation state jets. It preserves
the selected cutoff, pressure memory, physical widths, root and frequency.
The cutoff estimate below is analytical; manufactured algebra checks are
not a proof of the actual physical input bounds.

## The selected cutoff at every finite order

Use exactly the [frozen variance cutoff](../explicit_variance_cutoff/README.md):
phi(t)=exp(-1/t-1/(1-t)) on (0,1), zero elsewhere,
I=integral phi, K(t)=1-I^(-1) integral_0^t phi, with flat plateaus,
and zeta(v)=K(8(v-2)/delta-1). The earlier certificate supplies I^(-1)<729.
For u=1/t, w=1/(1-t), phi^(n)=P_n(u,w) exp(-u-w), with

    P_0=1,
    P_(n+1)=(-u² partial_u+w² partial_w)P_n+(u²-w²)P_n.

The degree is at most 2n. If L_n is the sum of absolute coefficients,
then L_(n+1) <= (2n+2)L_n, hence L_n <= 2^n n!.
For a monomial of degree i+j<=2n, u^i w^j exp(-u-w)<=i!j!<=(2n)!.
The elementary central binomial bound (2n)!<=4^n(n!)² gives

    ||phi^(n)|| <= 8^n(n!)³,
    ||K^(j)|| <= T_j := 729·8^(j-1)((j-1)!)³, j>=1; T_0=1. (1)

The flat extension satisfies the same inequalities. This is deliberately
coarser than the frozen first-five bounds, which remain unchanged. Unlike
a polynomial expansion to order 100000, (1) gives a symbolic bound at any
requested finite order without enumerating its monomials.

## Finite multivariate majorants for the actual amplitude

All jets below are normalized by multiindex factorials in the original
slow variables, before phase substitution. Let D be a finite downward-closed
index set, m its maximum total degree, and truncate every operation to D.
At a real base point let A_+(h), B_+(h) be nonnegative zero-constant
polynomials majorizing the increments of a and b_s. Assume a>=a_*>0
and |b_s|<=B_* on the region. These are explicit input requirements.
Use the original identities t_s=-b_s/a and v_s=a+b_s²/a. Define

    R(h)=sum_(k=0)^m (A_+(h)/a_*)^k,
    V_+(h)=A_+(h)+[(B_*+B_+(h))² R(h)-B_*²]/a_* .       (2)

Then every positive-degree v_s jet is bounded coefficientwise by V_+.
Indeed, reciprocal differentiation is the finite geometric expansion about
the actual a_0; replacing a_0 by its lower bound and |b_0| by B_* yields
(2). The subtraction removes only the constant coefficient. No convergence
of an infinite slow Taylor series is assumed.

Put g=delta/4 and H_*=2+delta/2. On the closed nonzero-cutoff region,
H=2+delta/2-v_s>=g; also 0<H<=H_* since v_s>=a>0. Define

    Z(h)=sum_(j=0)^m (T_j/j!) (8V_+(h)/delta)^j,
    W(h)=sum_(k=0)^m (V_+(h)/g)^k,
    S(h)=sqrt(H_*/a_*) Z(h) W(h) R(h).                  (3)

For the ACTUAL s=zeta(v_s)sqrt(H/a),

    |partial^alpha s|/alpha! <= [h^alpha] S(h), alpha in D. (4)

To prove it, compose the finite Taylor polynomial of K with the affine
rescaling of v_s, using (1). For the remaining factor, normalized Taylor
coefficients of (1+x)^(1/2) and (1+x)^(-1/2) have absolute value at most
one. The expansions about H_0,a_0 are thus dominated by W,R and the
value bound sqrt(H_*/a_*). Product coefficients give (4).
At v_s=2+delta/4 all derivatives of zeta vanish and H still has gap g;
the same product extends smoothly by zero. Beyond that threshold s and
all its jets vanish. At the lower plateau boundary K is likewise smooth.
One never evaluates sqrt(H) on the region where H may be negative.

Take S_+=S-S(0) and use it in Q65 equation (4), together with the
full pressure-input majorant P_+ for p=p2. This supplies an explicit
finite bound for the signed root from a,b_s,p2, rather than independent
s or mu jets. Downward-closed truncation preserves Q63's anisotropic
orders: no additional eta derivative is introduced by (2)--(4).
The retained compact floors allow a_*=A^(-1); |t_s|,v_s<=A and a<=v_s
also give |b_s|=a|t_s|<=A², so B_*=A² is sufficient. Higher jets of
a,b_s,p2 are NOT bounded by these value floors.

## Remaining inequality and scope

The next obligation is to enclose the positive-degree polynomials
A_+,B_+,P_+ on each original pre-modulation region and combine them with
the low-order normalized repair defect. In particular P_+ must come from
the full incoming five moments, including both pressure terms
4Aphys eta Pi and -(1-eta²)Pi_eta. No reset of Pi is allowed.
The derived constants may be huge; no smallness gate or uniform physical
pulse bound follows merely from finiteness. The axial-tail rate remains
conditional; the radial pressure tail, complete residual and blow-up are
unproved. Historical frozen evidence is not modified.

Run test_audit.py and audit.py --verify-record evidence.json --output a fresh
path. The replay checks provenance; it does not certify continuum estimates.
