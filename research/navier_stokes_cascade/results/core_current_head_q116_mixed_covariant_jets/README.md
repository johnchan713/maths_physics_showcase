# Q116: all fixed mixed covariant jets with harmonic residence

Parent: `8bbc5124f39769e58bdcdaa4c13663b310b0f322` (Q115).
On 3 November 2026 the inspected remote was Q115; research CI 41575057323
and field CI 41575057731 were running. Q114 field CI 41569012916 had passed.
The latest frozen evidence and core_growth_activation scope were reread.
The core result remains local to the leading core and first activation
collar at the unchanged width C_setup^-20. This checkpoint extends Q115's
projected first derivative to EVERY FIXED mixed order, including pulse-time
derivatives under explicit entrance compatibility. It does not establish
the actual source norm, a high-pole archive or a radial gain. All historical
frozen evidence, actual fields, cutoffs, roots and pressure stages remain.

## 1. Ordinary mixed derivatives and the exact triangular equation

Use Q115's SAME constrained normalized system on a common closed lifted
interval 0<=v<=L:

    u'=B_m u-Pi g, n.u=0, u(0)=0,
    N=n.n, h=n/N, Pi=I-h n^T,
    A0=-Pi K-h(n')^T, B_m=A0-(m^2 d+p)I,
    d=epsilon k^2 N, p=P'/P, g=f/P.                    (1)

Here n is real and bounded away from zero. Discrete labels, m,k, interval
and entrance are fixed before differentiating. D^alpha denotes ordinary
COMMUTING coordinate derivatives in parameters AND pulse time v, as in
Q95/Q96's raw whole-path norms. Derivatives of physical variable-coefficient
frames still require their original chart/connection conversions.

Set U_alpha=D^alpha u and, for alpha nonzero,

    Z_alpha=sum_(0<beta<=alpha) binom(alpha,beta)
                         (D^beta n).U_(alpha-beta),
    w_alpha=Pi U_alpha, U_alpha=w_alpha-h Z_alpha.     (2)

For alpha=0 put Z_0=0, w_0=u. This is just the differentiated constraint,
so Z_alpha involves only LOWER response orders. With Q=Pi(Kn+n')/N,

    w_alpha'=B_m w_alpha+F_alpha,
    F_alpha=Q Z_alpha
       +Pi sum_(0<beta<=alpha) binom(alpha,beta)
          [D^beta A0-(m^2 D^beta d+D^beta p)I] U_(alpha-beta)
       -Pi D^alpha(Pi g).                             (3)

Every forcing term in (3) is tangent. To prove it, differentiate (1),
then compute (Pi U_alpha)' using (2). The terms multiplying w_alpha
combine to B_m w_alpha, since its normal part is -h(n').w_alpha.
The remaining coefficient of Z_alpha is Pi(h'-A0 h)=Q. Scalar damping
and reference-rate terms on h disappear under Pi; they do not create an
uncontrolled m^2 Z_alpha term. This uses the actual moving projection,
not a frozen normal. No differentiated response at order alpha occurs
on the right except in the homogeneous B_m w_alpha.

At first parameter order, expanding the last two terms of (3) recovers
Q115's full pressure-numerator formula, including n.g. Keeping only
projected source jets before differentiating would invalidate this step.
At higher orders the full derivative of Pi g must likewise be retained.

## 2. Entrance hypotheses for pulse jets are essential

Pure parameter derivatives of u(0,a)=0 have zero entrance automatically.
For mixed derivatives containing r pulse derivatives, Q90's full
entrance recurrence implies zero entrance response jets through total
order J if the full mixed source jets through J-1 vanish at v=0, with
smooth coefficients on the entrance collar. This sufficient condition
also sets Z_alpha(0)=0 and w_alpha(0)=0 for |alpha|<=J.

These are hypotheses on the ACTUAL full unprojected source, not consequences
of sup norm tables. Q90 supplies their conditional cutoff-grammar route.
For zero extension of pressure through J, full source flatness through J
is additionally sufficient; velocity through J needs only J-1. Declaring
these orders in the program does not certify the physical endpoint data.
Flat exit source still does not imply a zero uncut exit response; the
original complete-potential/pressure cutoff and its residual remain.

Without entrance compatibility there is an exact obstruction to the
harmonic residence claim for pulse jets. Take n=e1,K=p=0,d=1 and constant
tangent g=e2. Then

    u(v)=-(1-exp(-m^2 v))e2/m^2, u'(0)=-e2.            (4)

The source jets are uniformly bounded in m, yet sup|u'|>=1 while ell_m
below tends to zero. One must retain a homogeneous entrance term or prove
compatibility. This is an operator comparison, not a canonical-source
counterexample. It leaves Q96's original ambient estimates intact.

## 3. A scalar majorant for every fixed mixed order

Retain the inherited normalized tangent two-time estimate

    ||V_m(v,s)||<=T exp[-(m^2-1)dmin(v-s)], T>=1,
    ell_m=L (|m|=1),
    ell_m=min(L,1/[(m^2-1)dmin]) (|m|>=2).              (5)

This comes from core_tangent_propagator_closure on the SAME background,
normal and Gaussian reference; its primitive validity is still required.
Let n_j,np_j,K_j,p_j,d_j,G_j bound every raw mixed jet of total order j
of n,n',K,p,d,g on the whole domain; q0 bounds 1/N. Thus np through J
requires the appropriate normal jets through J+1. Use binomial convolution
(a*b)_j=sum_i binom(j,i)a_i b_(j-i). A sufficient coefficient construction is

    D=n*n, inv_0=q0,
    inv_j=q0 sum_(i=1)^j binom(j,i)D_i inv_(j-i),
    h=n*inv, P_0=1, P_j=(h*n)_j for j>=1,
    A=P*K+h*np, Q0=(K_0 n_0+np_0)q0,
    C=L+4/(3 dmin), F=P*G.                            (6)

P_0=1 uses the norm of an orthogonal projection. The remaining P_j
bound derivatives of Pi; signs are discarded only in norm estimates.
For V_j and Zhat_j define the triangular recurrence

    Zhat_j=sum_(i=1)^j binom(j,i)n_i V_(j-i),
    V_j=T { F_j+L Q0 Zhat_j
            +L sum_(i=1)^j binom(j,i)(A_i+p_i)V_(j-i)
            +C sum_(i=1)^j binom(j,i)d_i V_(j-i) }
          +h_0 Zhat_j.                               (7)

Empty sums vanish, so V_0=T G_0. Duhamel in (3), induction in total
order and (2) give

    sup_|alpha|=j |D^alpha u| <= ell_m V_j.             (8)

Indeed Z_alpha is bounded by ell_m Zhat_j using lower orders;
ell_m<=L and m^2 ell_m<=max(L,4/[3 dmin])<=C absorb the two lower-order
sums. The reconstruction h Z_alpha supplies the last term of (7).
Multiindex Leibniz sums reduce to the displayed binomials because
sum_(|beta|=i,beta<=alpha) binom(alpha,beta)=binom(j,i).
There is no missing dimension factor in this Euclidean/operator max norm.

If the supplied native jets, T,L,q0 and dmin^-1 have band-uniform positive
polynomial bounds in S, (6)--(7) produce a finite such polynomial for each
fixed j. C is a positive SUM rather than a max, so this remains a positive
coefficient program even when L has polynomial dependence on S. There is
no exp(S) loss or numerical evaluation of giant factorials in this proof.
Harmonic independence of V_j requires harmonic-independent input bounds;
otherwise their actual m dependence remains. No arbitrary-source harmonic
summability follows. The code evaluates finite rational controls, not the
unverified native coefficient arrays at J=100005.

## 4. Full pressure and potential jets

Use the original exact formulas

    pi/P=i inv[(n^T K-(n')^T)u+n.g]/(km),
    potential/P=i h cross u/(km).                     (9)

With U_j=ell_m V_j, binomial convolution gives

    k |D^j(pi/P)| <= |m|^-1 {inv*[(n*K+np)*U+n*G]}_j,
    k |D^j(potential/P)| <= |m|^-1 (h*U)_j.             (10)

The pressure has the direct n*G term, with NO residence gain. For constant
n=e1 and g=a^2 e1, u=0 but the second normalized pressure parameter jet is
2i/(km). All mixed moving-normal np_j terms remain. Gaussian Bell
restoration, packet evaluation, full curl, spatial/temporal cutoffs and
physical derivative losses still follow Q95/Q96; (10) does not replace
them. Neither local pressure nor its derivatives reset Pi1/Pi2 history,
axial A1 velocity, the compensator or K0 in the current radial mean pressure.

## 5. Verification and frontier

Exact symbolic controls differentiate an independent moving-normal
manufactured solution through second parameter, mixed pulse/parameter,
and third pulse order, then check (2)--(3) at rational points. A separate
damped polynomial-source model has the required entrance jets and checks
(8) numerically at moderate harmonics. Exact controls exhibit (4), the
normal-source pressure term and multiindex combinatorics. They are finite
diagnostics, not enclosures of the canonical fields. The analytical proof
is (2)--(3), the entrance induction and Duhamel in (7)--(8). Frozen replay
retains Q115 and its historical evidence graph.

Next supply the actual whole-path native inputs and compatible high-pole
proof domain for this recurrence, then restore the complete packet/curl
bounds and rebuild covariance/means at their current orders. The radial
request still needs source order 100005 and profile jets through 100046;
none were materialized or certified here. The inner/outer cone join,
coefficient/reserved-mean supports, all-scale actual-pulse bounds, complete
viscous/cutoff/interaction residual improvement, infinite realization and
blow-up remain unproved.
