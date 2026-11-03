# Q115: covariant first parameter jet and harmonic residence

Parent: `c534ecc17e4c4a75f407f2c7cd9798036aa55103` (Q114).
On 3 November 2026 the inspected remote was Q114; research CI 41569012859
and field CI 41569012916 were running. Q114 bounds coefficient algebra,
not native source leaves. This checkpoint advances one such interface:
the first parameter derivative of Q96's normalized constrained response.
It retains the original root, widths, pressure, profiles and accepted domain.
Historical frozen evidence is unchanged. This is an analytical identity
and conditional norm estimate, not a computed canonical-field enclosure.

## 1. Use the existing two-time tangent estimate

core_tangent_propagator_closure and core_full_path_ode_reconstruction
already establish, conditional on their retained primitive/frame bounds,

    ||V_m(v,s)|| <= T P(v)/P(s) exp[-(m^2-1)dmin(v-s)],
    0<=s<=v<=L, m a nonzero integer, dmin>0.              (1)

This is the SAME normal, K and reference P as Q96, on every subinterval;
it is not a bound merely for one growing column. Labels, m, carrier k,
entrance and the common lifted interval are fixed under the parameter
derivative below. On the normalized tangent equation the P ratio cancels.
Define the sufficient residence bound

    ell_m=L for |m|=1,
    ell_m=min(L,1/[(m^2-1)dmin]) for |m|>=2.             (2)

It bounds the integral of the decay factor in (1). The first derivative
estimate below retains (2). No new value-domain selector, first band or
Gaussian minimum is introduced. Actual validity of the primitive bounds
in (1) remains an input; this is not a replacement proof of that theorem.

## 2. Differentiate the constraint before applying a tangent inverse

Let N=n.n, h=n/N, Pi=I-h n^T, d=epsilon k^2 N, p=P'/P,
A0=-Pi K-h(n')^T and B_m=A0-(m^2 d+p)I. Q96's actual normalized response
is

    u'=B_m u-Pi g, n.u=0, u(0)=0, g=f/P.               (3)

Prime means path time. Let a be one continuous slow/transverse parameter,
with a common smooth lifted domain. Set z=n_a.u and w=Pi u_a. Exactly

    u_a=w-h z, n.w=0, w(0)=0.                           (4)

Thus u_a need not be tangent. Define the tangent vectors

    q_a=Pi n_a/N, Q=Pi(Kn+n')/N,
    J=(n^T K-(n')^T)u+n.g.

Direct differentiation of (3), followed by projection and (4), gives

    w'=B_m w+F_a,
    F_a=q_a J-Pi K_a u-(m^2 d_a+p_a)u-Pi g_a+Q z.     (5)

Every term of F_a is tangent. To check the moving-plane terms, use
Pi h_a=q_a and Pi(h'-A0 h)=Q (the latter is also Q79 (10)).
The projected derivative of A0 contributes
q_a[(n^T K-(n')^T)u]-Pi K_a u. The derivative of Pi g contributes
q_a(n.g)-Pi g_a. These yield (5), including both signs and the Qz term.
There is no n'_a term in this PROJECTED velocity identity; it is killed
by Pi h=0. That does not delete n'_a from differentiated pressure.

J is precisely the full normalized pressure numerator. Keeping only
projected forcing in this calculation would incorrectly discard n.g.
For a purely normal source, its contribution cancels the corresponding
part of Pi g_a in the VELOCITY equation. Those two terms must both be
present for this cancellation. Local pressure can still be nonzero.

## 3. Explicit first-jet coefficient bound

Use whole-path Euclidean/operator upper bounds n0,np0,na0,K0,Ka0,
pa0,da0,G0,Ga0 for n,n',n_a,K,K_a,p_a,d_a,g,g_a and q0>=1/N.
All are nonnegative, band-uniform norm coefficient premises on the SAME
proof domain, not arbitrary callbacks certified by smoothness. Put

    R_m=T ell_m, U=R_m G0,
    qa=na0 q0, Q0=(K0 n0+np0)q0, h0=n0 q0,
    F1=Ga0+qa n0 G0
         +[qa(n0 K0+np0)+Ka0+pa0+m^2 da0+Q0 na0]U.

Duhamel applied to the TANGENT equation (5), then (4), proves

    sup |u|<=U,
    sup |w|<=R_m F1,
    sup |u_a|<=R_m F1+h0 na0 U.                         (6)

In particular

    m^2 ell_m <= max(L,4/[3 dmin])                     (7)

for every nonzero integer m. For |m|>=2 this follows from
m^2/(m^2-1)<=4/3; for |m|=1 it is L. Hence if the displayed source and
coefficient bounds are uniform in m, (6) retains an ell_m factor even
in the term with differentiated m^2 d. This is a CONDITIONAL first-jet
harmonic improvement over using Q96's undamped ambient length everywhere.
If the actual g_m norms grow with m, their dependence remains in G0,Ga0;
no summability of the actual harmonic source is inferred from (7).

An ambient residence substitution would be false. For n=e1,K=0,p=0,d=1,
the ambient normalized generator is diag(0,-m^2,-m^2). Its normal mode
has no damping, and constant normal forcing e1 gives y(v)=v e1.
The ambient inverse therefore cannot have the decaying ell_m bound.
This is an exact operator comparison, not a canonical-source counterexample.
Q96's original ambient estimate remains valid. Projection plus constraint
reconstruction, rather than simply replacing L by ell_m, is essential.

## 4. Full normalized pressure derivative

The exact formula is pi/P=i J/(k m N). Let npa0 bound n'_a and let Ua
be (6). With J0=(n0 K0+np0)U+n0 G0, differentiation gives

    |(pi/P)_a| <= (1/[k |m|]) {
        2 q0^2 n0 na0 J0
        +q0[(n0 Ka0+K0 na0+npa0)U
                +(n0 K0+np0)Ua+na0 G0+n0 Ga0] }.       (8)

Here |(1/N)_a|<=2q0^2 n0 na0. The moving-normal mixed jet n'_a and full
normal-source derivative are indispensable. A constant normal with
g=a n has u=0 yet (pi/P)_a=i/(km). The code returns k times the sufficient
bound (8); k is the original frozen carrier. Q96's k^-1<=sqrt(epsilon)
transfer and full Gaussian Bell restoration still apply. Derivatives of
raw pi, rather than pi/P, require those extra reference-profile terms.

The local pressure in (8) does not reset cumulative mean pressure, Pi1/Pi2,
axial A1 velocity, scalar compensation, or the K0 memory in Q110/Q114.
No pressure stage is substituted for another in the reduced radial source.

## 5. Scope, checks and next obligation

Exact symbolic controls construct a moving normal, nontrivial K, varying
damping and reference rate, and a tangent manufactured response. They
verify (4)--(5) at exact rational points with and without normal forcing.
Separate controls show the undamped ambient obstruction and the nonzero
normal pressure derivative. Rational tests verify (7), negative harmonics,
all required norm leaves and pressure contributions. These diagnostics
are not continuum enclosures of the selected fields. The analytical
proof is projection/Leibniz in (4)--(5) and Duhamel in (6).

Next extend the covariant recurrence to the remaining mixed jets or supply
the actual whole-path source/coefficient bounds needed in (6)--(8). This
first-jet result does not meet the radial source request J=100005, rebuild
its common pole, prove all harmonic summability or close the cone join.
Native coefficient/reserved-mean supports, uniform actual-pulse bounds on
all physical scales, complete viscous/cutoff/interaction residual
improvement, infinite realization and blow-up remain unproved.
