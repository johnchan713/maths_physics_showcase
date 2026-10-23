# Canonical nonzero entrance: explicit homogeneous mixed-jet closure

Parent: `9e8f33555745e8a5a07eeeeadb5c10d8798c5a74`.
Status: `conditional-canonical-entrance-and-homogeneous-jet-closure`.
The parent's hosted run 40016330531 was in progress when inspected on
23 October 2026; its reference-normal-floor step had passed. No historical
evidence, field, root, band, length or cutoff is changed.

This checkpoint replaces the independent homogeneous wave-jet input in
`core_primary_root_budget` by an explicit recipe from the retained full
phase/background jets and the original accepted-tail propagator. It makes
the qualitative fixed-jet existence statement in `canonical_primary_foundation`
quantitative for its ORIGINAL nonzero canonical entrance. It is not a new
existence theorem, evaluation of actual background suprema, or new value/cone
acceptance proof. The original canonical and inherited families stay distinct.

## 1. Retained domain and coefficients

Fix an original lifted label and its closed interval 0<=v<=L. Freeze the
discrete carriers, epsilon, lambda0, c0, tilt, L, representative and entrance
v=0 under continuous parameter derivatives. The profile is normalized by
P(L/2,x)=1 and p=P'/P. Work in the same fixed coefficient chart as the parent;
all input jets include its full common-chart/lift conversions. Endpoint or
moving-frame chain terms must be converted BEFORE using this statement.
Write S=ell^2 and let Length(S) bound L. All jets below are raw maxima over
mixed continuous derivatives of a fixed total order; x denotes ANY tuple of
continuous parameters, not just the single manufactured test parameter.

Use `core_phase_coefficient_budget`'s full normal n_j, moving-normal n'_j,
matrix K_j and profile-rate p_j. Its p_j already retains the actual reference
Gaussian and frozen-label semantics. A leading first-jet comparison alone
does not supply these full inputs. Put

    N=|n|^2, h=n/N, H=h n^T, d=epsilon*k^2*N,
    A0=-K+h(n^T K-nprime^T), A=A0-dI,
    Ahat=A+(d+p)H.                                      (1)

This is the original homogeneous primary |m|=1. Ahat agrees with A on
tangent solutions. Since nprime^T+n^T Ahat=p n^T, its normal response is
(n.t)'=p(n.t). `core_first_sourced_inverse_budget`'s lift construction gives

    ||Vhat(v,w)|| <= Kamb(S) P(v)/P(w), 0<=w<=v<=L,
    Kamb=T[1+Length Q n_0]+h_0 n_0,
    Q=(A0)_0 h_0+h_1.                                  (2)

Here T is the ORIGINAL tangent bound derived by
`core_tangent_propagator_closure`. The production factory derives it from
the original frame/length gates, rather than accepting a new T. It uses
`core_reference_scaled_normal_floor`'s smaller reciprocal without recomputing
the old CE, Cd, T or length. Recurrences for h, Ahat and Kamb are unchanged.

## 2. Canonical entrance without identifying the inherited column

Let tau=n_tan, J=|tau|^2, nr=n_r and rot(tau)=(tau_z,-tau_theta).
The manuscript's Br=U Mr gives the physical entrance vector

    b0=Br(0)eplus
      =er-nr*tau/J+c0*wr(0)*rot(tau)/sqrt(J),
    t(0,x)=P(0,x)b0(x),
    wr(0)=sqrt(1+u_actual^2/4).                         (3)

The actual fixed tilt satisfies |u_actual|<=u_upper (the code's `u`), so
|c0 wr(0)|<=Centry=cmax(1+u_upper/2). Both are frozen-label factors.
Replacing wr(0) by sqrt(1+(nr/|tau|)^2) would construct the inherited
entrance instead. This checkpoint never makes that replacement.

The preceding reference-floor PROOF bounds both tangential and full norms:

    1/J <= qtan=min(4/bmin^2, 1/[(1-theta)^2 beta]).      (4)

A bound for 1/N alone would NOT imply (4); n=(1,delta,0) is a counterexample
as delta tends to zero. This use depends on the preceding tangential error
proof, not on relabeling an ambient reciprocal. Let D_j=(n*n)_j, which bounds
the raw jets of J as well as N. With raw convolution
(a*b)_j=sum_i binom(j,i)a_i b_(j-i), define

    R_0=qtan,
    R_j=qtan sum_(i=1)^j binom(j,i) D_i R_(j-i),
    W_0=max(1,qtan),
    W_j=(1/2)(W*R*D_shift)_(j-1), D_shift_i=D_(i+1),
    B_j=1_(j=0)+(n*n*R)_j+Centry(n*W)_j.                (5)

These bound J^-1, J^-1/2 and D^j b0. The inverse-square-root recurrence
comes from D(J^-1/2)=-(1/2)J^-1/2 J^-1 DJ. The same raw multibinomial
grouping applies to every mixed derivative. No tangential component floor
or actual-wave jet is introduced.

Do not discard entrance Gaussian derivatives. Because the endpoints and
L are frozen,

    log P(0,x)=-integral_0^(L/2) p(v,x) dv,
    ell_j=(Length/2)p_j, j>=1,
    G_0=1,
    G_j=sum_(i=1)^j binom(j-1,i-1) ell_i G_(j-i),
    I_j=(G*B)_j.                                       (6)

I_j bounds |D_x^j[P(0)b0]|/P(0). In strictly fixed reference coordinates
the parameter derivatives of p can vanish, but setting them to zero in a
different representation requires a proof. Formula (6) is conservative
and retains them. There is no P(0)^-1 majorant: the common P(0) cancels
against (2). For moving entrances, missing endpoint chain terms would make
(6) inapplicable; this theorem explicitly uses the original fixed entrance.

## 3. Nonzero-data all-mixed derivative induction

Let Ahat_j be the parent's positive generator-jet recurrence. Define

    U_0=Kamb I_0,
    E_j=Kamb[I_j+Length sum_(i=1)^j binom(j,i) Ahat_i U_(j-i)],
    V_j=sum_(i=0)^(j-1) binom(j-1,i) Ahat_i U_(j-1-i),
    U_j=E_j+V_j, j>=1.                                (7)

Parameter-only derivatives satisfy the differentiated homogeneous ODE
with initial data D_x^j[P(0)b0]. Duhamel with (2), (6) gives E_j: every
integrand's P(w) cancels the propagator denominator. For any derivative
containing a path derivative, differentiate t'=Ahat t j-1 times and use
the lower total-order bounds; this gives V_j. Their sum bounds the maximum
over both cases. Induction gives |D^I t(v,x)|<=P(v,x)U_j(S), |I|=j,
on the WHOLE retained interval. Each fixed U_j is a finite positive scale
polynomial. No exponential in S, Gaussian reciprocal, or independent W_j
is introduced. Fixed-order coefficient sizes may be enormous.

Homogeneous pressure and potential are retained exactly:

    k*pi=i(n^T K-nprime^T)t/N,
    k*Avec=i(n cross t)/N.                             (8)

The inherited `pressure`/`potential` recurrences now use (7). The source
is exactly zero for this homogeneous primary only. This does not authorize
projecting away n.f in any forced mode, dropping the moving-normal term,
or truncating any later global pressure primitive, moment or mean memory.

## 4. Join to the ORIGINAL root norm slot

For the original product cutoff chi and coefficient-chart dimension d=5,

    Wroot_j=(d^j/j!)(chi*U)_j.                          (9)

Indeed sum_(|I|=j)1/I!=d^j/j!. Equation (9) bounds J_j(chi*t)/P, NOT
J_j(t/P). This distinction retains derivatives of P. `root_wave_join`
replaces this NORM slot only, keeping H, target, curl, loghalf callbacks
and K0, beta, R, mass and width identical. Thus the existing root, signed
difference and full curl norm recurrences can use (9) on their ORIGINAL
canonical value domain. No physical root is rebuilt, no root/Haar gap is
proved anew, and no edge constant is applied a second time.

## 5. Evidence, frontier and reproduction

`controls.py` checks a manufactured moving normal n=(v+x,1,0) and exact
homogeneous solution t=P(1,-v-x,1) on [0,1]^2, with a nonzero pressure
k*pi/i=P(1+v+x). Its exact rational coefficient jets and every mixed
derivative through order3 are enclosed on the whole box by numerator l1
and proved positive denominator floors. These are manufactured exact
controls, not sampled diagnostics or actual background enclosures. The
tests detect loss of moving-normal pressure, Gaussian entrance derivatives,
cutoff derivatives and the nonzero initial datum; a separate varying
tangential norm checks the inverse-square-root chain. Production-factory
controls retain the original accepted gates and root-budget objects.
Frozen historical source hashes and eight parent records are bound anew.

An initial investigation of edge-weight interpolation did NOT yield a new
obstruction: `core_initial_haar_row_budget` already derives the initial
fine weighted row from the full mean equations and lower primary/map
inputs. That established conditional implication is retained, together
with `core_signed_prefix_source_budget`'s generated-source recipe. This
checkpoint closes a LOWER homogeneous wave input quantitatively; it does
not reopen those completed joins.

The next lower obligation is to make the canonical Haar-column H_j jet
input equally explicit from (7), the actual covariance integrals, all
moving maps/cutoffs and the retained Gaussian/edge normalization, or to
identify the precise inequality blocking that composition. The validity
of full actual phase/background jets and same physical domains remains
a premise. Complete generated-source/weighted-target estimates still
depend on all their retained lower inputs. Global pressure/moment and
reserved-mean supports, fresh axial/full-memory pressure means, and every
viscosity, cutoff and interaction residual must remain in the full ledger.
No unconditional pulse bound, full residual improvement or blow-up follows.

Reproduce with `python test_audit.py -v`, then
`python audit.py --verify-record evidence.json --output /fresh/path.json`.
Install the declared SymPy requirement if necessary. Passing tests check
the exact controls/provenance; the analytical induction above supplies the
conditional continuum statement.
