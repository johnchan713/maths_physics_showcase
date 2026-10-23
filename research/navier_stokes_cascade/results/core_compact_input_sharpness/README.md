# A sharp finite-jet boundary for the actual compact operators

Parent: `df3b93b4f07489e4f2806f359200dc34f5537d84`.
Its fresh-mean/current-tail CI check passed when inspected on 23 October
2026; the complete workflow was still running without a failed step.

The parent requests tens of millions of source derivatives for its physical
A0/A1 tail estimates and leaves their ACTUAL bounds unverified. This
checkpoint audits whether those requests can be removed simply by improving
the scalar compact-operator inequality. There is a sharp obstruction:

    N_m(A_e f) <= C M^(-p) N_J(f), uniformly in M,
    for every smooth compact zero-Haar f, REQUIRES J >= m+2p.       (1)

The inherited sufficient integer rule J=m+2p+3 is within THREE derivatives
of this necessary boundary. Even allowing any fixed polynomial in the
original S does not remove a strict power deficit. Thus the huge requests
cannot be dismissed as a loose isotropic finite-jet estimate. A proof using
special structure of the SAME actual sources may still do better.

This is an analytical scalar-operator obstruction, not a counterexample to
the canonical source, its all-fixed-jet class, a physical residual estimate,
or Navier--Stokes. We construct comparison sources explicitly; none is
inserted into the actual field. Historical evidence, roots, carriers,
width C_setup^-20, cutoffs, h, pressure, physical scales and accepted band
are unchanged. Actual input validity, complete residual decay, infinite
realization and blow-up remain unproved. The next generic signed cycle is
not advanced here: its incoming high-jet obligation has not been discharged.

## 1. Retain the original moving compact kernel

Use U=R^d, d=d_r>0, vr=(1,-b), b=sqrt(2)-1. At fixed M,c the
original compact error, for e=0 or 1, is

    A_e f(U,y)=k_e(U) chi_U(U,Z,T)
          int_all f_e(U+a) f(U+a,y+Ma vr,Z,T) da,
    f_e=R^(e-d+1)/d, k_e=d R^(d-1-e).                 (2)

These are the SAME profiles and moving cutoff used in the parent's raw
kernel recipe and `core_fast_inverse_budget`, not a replacement unweighted
operator. On a positive fixed radial hull f_e,k_e are smooth and nonzero.
Since chi goes from zero to one, at any fixed slow parameter its derivative
is nonzero somewhere. Fix one such point U*, within its existing output
support. Write c*=abs(k_e(U*)chi_U(U*))>0. It may be enormous or tiny;
only positivity and independence from M matter for the obstruction.

Choose a fixed nonnegative C-infinity radial bump phi, supported strictly
inside the original positive input hull, with integral phi=1. Choose a
fixed compact slow bump B(Z,T) equal to one at the chosen slow point.
No choice of a new physical chi is involved. All derivatives of
f_e^-1 phi B are finite; their actual numerical values are not needed.
On this compact interior input support any retained flat shell weight has
a positive lower bound, so its fixed finite weighted norm changes these
constants only. The proof does not assert a uniform edge constant on a
larger hull or replace the parent's weighted source theorem.

## 2. Pell modes cancel the radial phase exactly

For j>=1 put

    P_j+Q_j sqrt(2)=(1+sqrt(2))^j,
    k_j=(P_j-Q_j,Q_j),
    omega_j=vr.k_j=P_j-sqrt(2)Q_j
             =(-1)^j/(P_j+sqrt(2)Q_j).               (3)

Here |k_j|_1=P_j, P_j/2<=Q_j<=P_j, omega_j!=0 and
abs(omega_j)<=1/P_j. For ANY M>=1 select the first P_j with P_j^2>=M.
The recurrence P_next=P+2Q<=3P gives

    sqrt(M)<=P_j<3sqrt(M), abs(M omega_j)<=P_j.        (4)

Thus no arithmetic matching between Pell numbers and rounded carriers is
assumed. In particular the comparison works at each unchanged actual band
value of M once it is >=1. This selection is a per-band scalar comparison,
not a construction of one coherently glued smooth canonical field.

For amplitude A define the real zero-Haar source

    f_M(U,y,Z,T)=A f_e(U)^(-1) phi(U) B(Z,T)
                     cos(2pi[k_j.y-M omega_j U]).     (5)

It has smooth whole-line zero extension, compact slow/radial support and
zero auxiliary Haar at EVERY radius and slow point. In particular every
integrated Haar radial moment vanishes. It meets the scalar A1 zero-Haar
condition and the compensated A0 zero-mode condition; the compensator
does not eliminate this oscillatory source.

In (2) the phase becomes

    k_j.(y+Ma vr)-M omega_j(U+a)
                 =k_j.y-M omega_j U.                (6)

Consequently the exact output, including the SAME moving cutoff, is

    A_e f_M=A k_e chi_U B cos(2pi[k_j.y-M omega_j U]). (7)

There is no small integral left: the radial oscillation of the SOURCE
cancels the drift oscillation. Deleting the source's radial derivatives
from the input cost would miss precisely this resonance. The profile
f_e^-1 in (5) is essential; it cancels f_e in the actual kernel. The real
conjugate pair, rather than just a complex non-real mode, gives the same
conclusion. No total wave covariance or pressure term is removed.

## 3. Prove the finite-jet lower boundary

Let N_J be the maximum RAW coefficient derivative through total order J
in (U,Z,T,y1,y2), as in the inherited compact proof. Derivatives hold M,c
fixed. Set A_J=N_J(f_e^-1 phi B), enlarged to at least one if necessary.
All phase frequencies in (5) are <=2pi P_j. Since P_j>=1 and 2pi<8,
Leibniz gives the sufficient bound

    N_J(f_M) <= abs(A) 16^J A_J P_j^J.                (8)

At the chosen slow point and U*, take m derivatives in y2, then a torus
supremum. Equation (7), Q_j>=P_j/2 and 2pi>=2 give

    N_m(A_e f_M) >= abs(A) c* P_j^m.                 (9)

Set A=P_j^-J. The input in (8) is uniformly bounded, whereas

    N_m(A_e f_M) / [M^-p N_J(f_M)]
       >= c*/(16^J A_J) M^p P_j^(m-J)
       >= C_(m,J,p) P_j^(m+2p-J).                  (10)

The last constant uses M>P_j^2/9 from (4), with p>=0. If
J<m+2p, the right side diverges. This proves (1). It is valid for any
fixed real p>=0; the original upper rule uses positive integer p.
At equality this argument gives no obstruction and supplies no upper
bound. It does NOT prove that the extra three derivatives of the existing
upper estimate can be removed.

The retained upper proof integrates radially p times, loses p torus
derivatives to the irrational divisors and uses three more derivatives
for Fourier summability. The lower family shows why the first TWO costs
cannot simply be counted as one. Taking M=P_j^2 is a convenient exact
subsequence; (4) also proves the conclusion over arbitrary M tending to
infinity. Any fixed log(M) polynomial multiplying (1) loses to a strict
positive power in (10).

This is a FINITE-JET statement. With A=P_j^-J, higher source derivatives
grow like P_j^(r-J). Therefore (5) does not supply uniform bounds at ALL
fixed r with that same source rate. It is not a counterexample to the
original all-fixed-jet class, whose theorem can legitimately request the
higher derivatives. Extra relations such as divergence, source equations,
special Fourier decay or coherent band gluing are not imposed by this
scalar comparison. No necessary lower bound for a restricted canonical
class is claimed. In particular this proof concerns A1 applied to the
residual in the exact fast-time ERROR, not a separate lower bound for the
velocity A1 Nt^-1 f, whose temporal multiplier also matters.

## 4. Join the original scale flooring and the parent's actual requests

For the original carriers

    M=c_i^rho_g epsilon^-kappa,
    (Tg S)^-1<c_i<=S^-1, S=ell^2, epsilon=2^(-h ell).

Give (5) amplitude epsilon^a P_j^-J. For J>=m, (4), (8)--(9) give

    N_J(f_M)<=16^J A_J epsilon^a,
    N_m(A_e f_M)>=C epsilon^[a+kappa(J-m)/2]
                                      S^[rho_g(J-m)/2].         (11)

The bound follows from M<=S^-rho_g epsilon^-kappa; the other flooring
inequality ensures M tends to infinity and covers the analogous J<m
case with a fixed Tg factor. If a+kappa(J-m)/2<target, no fixed
polynomial P(S+K) can bound this output by epsilon^target P(S+K).
Indeed its ratio contains exp(h log(2) ell [target-a-kappa(J-m)/2]),
which dominates every polynomial in ell. This gives the necessary
finite-jet boundary

    J >= ceil(m+2(target-a)/kappa).                  (12)

Equality is not decided: the retained S factors can be polynomially
absorbed. This is still a comparison-class theorem, not a failure of the
actual input inequality or a lower bound for the actual residual.

For the parent's normalized tail join use target=C0=2-2kappa, a=B=1-kappa
for A0 and a=E=2-3kappa for A1. For its physical transfer ROUTE use
target=(N+3/2+2h+2m)/h. That is the inherited sufficient derivative loss;
we do not prove it is necessary for every physical argument. At
kappa=1/100000, h=1/200, m=N=0 the exact comparison is:

| Request | Necessary J for this scalar finite-jet route | Retained raw J |
|---|---:|---:|
| Normalized A0 | 199998 | 200003 |
| Normalized A1 | 2 | 7 |
| Physical A0, via the retained transfer | 60200002 | 60200007 |
| Physical A1, via the retained transfer | 60000006 | 60000011 |

For a nonnegative target gain and the strict integer rule
p=floor((target-a)/kappa)+1, the raw J=m+2p+3 exceeds (12) by
between three and five derivatives. The +5 entries above arise from an
integer target gain and a STRICT positive exponent margin. They do not
turn the upper rule into a necessary minimal p theorem. These are source
indices only; older nested source/root/normal requests remain larger and
are not erased. No factorial or huge norm program is evaluated.

## 5. The precise remaining actual-input inequality

The unresolved source inputs remain, at the respective requested J,

    N_J(Qcirc gr_final) <= epsilon^B P_J(S+K_all),
    N_J(Qcirc Epostwave_z) <= epsilon^E Q_J(S+K_all),   (13)

with the original full support/endpoint and weighted-to-raw conversion
contracts, SAME profile kernels, all archived shifts and current common
domains. Here gr_final includes old means, the actual new fast velocity
gd-A1 gd, the fresh slow response and ALL six current eight-mode wave
covariances. Epostwave_z includes full forced and reconstructed pressure;
no pressure scalar rate is substituted for a pressure-field norm.

The parent supplies conditional one-way norm programs for (13), not their
independent actual evaluation or enclosure. This checkpoint proves that
simply replacing their huge finite-jet inputs by low isotropic norms is
not justified by compactness, zero Haar, exact radial moment compensation,
or the original carrier rounding. It does not establish that (13) fails.
Further progress must either validate the needed actual high-jet programs
and support identities, or prove extra quantitative structure of those
SAME sources that excludes (5) and improves their operator estimate.
Then advance the next generic signed cycle on the expanded domains.
The parent's regular order E and negative physical exponent
-3/2-3h*kappa-2m are unchanged. Higher full regular orders, a controlled
infinite diagonal, localization/forcing regularity and velocity growth
remain separate obligations before any blow-up claim.

Exact checks use Pell integer identities, full shifted-phase algebra,
real zero-Haar Fourier pairs, independently integrated finite-regularity
polynomial bump fixtures, raw derivative bounds, and the parent's exact
normalized/physical requests. The polynomial fixtures are NOT C-infinity
canonical sources; the analytical proof above uses an arbitrary fixed
C-infinity bump. Passing these checks is not a continuum pulse enclosure.

Reproduce: `python test_audit.py -v`, then `python audit.py
--verify-record evidence.json --output FRESH_PATH.json`.
