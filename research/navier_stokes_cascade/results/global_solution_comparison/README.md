# Pressure recovery, comparison and the conditional exclusion implication

This checkpoint proves a GENERAL COMPARISON LEMMA: a smooth, spatially
compact solution u on every [0,T], T<1, agrees with any smooth solution v
having the same force and initial velocity and a uniform spatial L2 bound
on that interval. No decay or growth assumption on the competitor's
pressure, and no global derivative or dissipation bound for v, is added.

The expanding-ball calculation closes with

\[
\tfrac12 E_R'(t)+\tfrac12 D_R(t)^2
 \le \|\nabla u(t)\|_\infty E_R(t)+C_T R^{-4/3},\qquad E_R(0)=0.\tag{1}
\]

This implies E_R(t)<=K_T R^(-4/3) and then v=u. Applied to the preceding
[localized field and force](../localized_force_extension/README.md), it
completes the CONDITIONAL downstream argument excluding a global smooth
bounded-energy solution: the inherited core growth contradicts smoothness
on a compact neighborhood of (0,1).

This is not an unconditional blow-up theorem. The leading profile,
background, all-order physical estimates and global auxiliary geometry
remain inherited conditional inputs. In particular, the construction of
globally compatible separated label supports has not been independently
verified. The new comparison proof does not validate these premises.
Numerical materialization of enormous constants is not logically required
for an analytic proof; validity of the continuum arguments is the issue.

Parent: `b510fd042c9f41709b5e685cd84ee865352e5f31`. Source under audit:
the proposed [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(10.15)--(10.21), Lemma 10.5. Frozen PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The proof below derives its comparison argument and retains the sharper
R^(-4/3) remainder permitted by the same estimates. Earlier evidence and
GOAL.md are unchanged.

## 1. Hypotheses and the exact difference equation

Fix T<1. Assume u,p and v,P are smooth solutions at viscosity one on
R^3 x [0,T], with the SAME smooth force and initial velocity. Suppose u
and p have one fixed compact spatial support on this interval, and

\[
\sup_{0\le t\le T}\|v(t)\|_2<\infty.\tag{2}
\]

For the inherited construction these properties are supplied conditionally
by the previous checkpoint. For the comparison lemma they are ordinary
stated hypotheses, unrelated to the special pulse construction. Smoothness
means all derivatives on each fixed compact set, not bounded derivatives
at spatial infinity. Set

\[
w=v-u,\quad\pi=P-p,\quad
g_{ij}=v_i v_j-u_i u_j=w_iw_j+w_i u_j+u_i w_j.
\]

Uniformly on [0,T], ||w||2<=W_T and sum_ij ||g_ij||1<=G_T, with finite
constants. The difference equation, with (div g)_i=sum_j partial_j g_ij,
is exactly

\[
w_t+\operatorname{div}g=\Delta w-\nabla\pi,\qquad\operatorname{div}w=0,
\tag{3}
\]
\[
w_t+(v\cdot\nabla)w+(w\cdot\nabla)u=\Delta w-\nabla\pi.
\tag{4}
\]

Neither w_t in global L2 nor grad v in global L2 is assumed. Both will
be unnecessary for the cutoff argument.

The only singular-integral theorem used below is boundedness of the
classical Riesz transforms on Lp for 1<p<infinity, applied twice at
p=3/2. An author-hosted primary reference is A. W. Knapp,
[*Advanced Real Analysis*, Chapter III, sections 5--6](https://www.math.stonybrook.edu/~aknapp/download/a2-1-realanal-clickable.pdf).
Its constants depend on p and dimension, not R or the fields. The explicit
kernel estimate, weighted Sobolev bound and interpolation used here are
derived below; this standard theorem is not replaced by numerical sampling.

## 2. Recover the pressure GRADIENT without assigning pressure decay

Let R_i have multiplier i xi_i/|xi| and define

\[
\pi_* =\sum_{i,j}R_iR_jg_{ij},\qquad
\widehat{R_iR_jf}(\xi)=-\frac{\xi_i\xi_j}{|\xi|^2}\widehat f(\xi).
\tag{5}
\]

The sign convention for a single Riesz transform is immaterial to its
double product, but the displayed double-product sign is essential. For
L1 data the Fourier transform is bounded. Since
integral_R3 (1+|xi|^2)^(-s) dxi<infinity for s>3/2, (5) defines
pi_* in L-infinity_t H^(-s); use s=2. It satisfies

\[
\Delta\pi_*=-\sum_{i,j}\partial_i\partial_jg_{ij}=\Delta\pi.
\tag{6}
\]

Equation (6) alone is insufficient: a linear harmonic pressure has a
nonzero constant gradient. We remove that ambiguity using the equation
AND the energy bound. For every a in C_c^infinity(0,T), integrate (3)
in time to get

\[
\int a\nabla\pi\,dt
 =\Delta\int aw\,dt+\int a' w\,dt-\operatorname{div}\int ag\,dt.
\tag{7}
\]

These terms belong to H^-2, L2 and H^-3 respectively. The last follows
from L1 -> H^-2 followed by one derivative. Also int a grad pi_* belongs
to H^-3. Thus H_a=int a(grad pi-grad pi_*) dt is in H^-3 and Delta H_a=0.
Its Fourier transform is a weighted L2 FUNCTION supported at {0}; a
function supported on that zero-measure set is zero. Therefore

\[
\nabla\pi=\nabla\pi_*\quad\hbox{as a space--time distribution on }(0,T).
\tag{8}
\]

This does not assume pi itself is tempered. Its time-averaged gradient is
tempered because (7) proves it. Nor does it differentiate w in global L2.

The finite-energy requirement matters. Without it, v=t e1, P=-x1 and
u=0,p=0 are smooth solutions with zero datum and zero force but different
velocities. Their difference pressure is harmonic and their kinetic energy
on a radius-R ball is (4pi/3)t^2 R^3, unbounded as R grows. They are
excluded exactly by (2), not by an added pressure condition.

## 3. Weights and a Sobolev estimate requiring only LOCAL derivatives

Choose 0<=phi<=1 smooth, equal to one on the unit ball, supported in the
radius-two ball. Write phi_R(x)=phi(x/R), chi_R=phi_R^8, R>=1, and

\[
E_R=\int\chi_R|w|^2,\quad
D_R=\left(\int\chi_R|\nabla w|^2\right)^{1/2},\quad
B_R=\|\phi_R^4 w\|_6.
\tag{9}
\]

For each fixed R the cutoff field is smooth and compact; no global
gradient assumption on w is needed. The homogeneous Sobolev bound gives

\[
B_R\le C\|\nabla(\phi_R^4w)\|_2
       \le C(D_R+R^{-1}W_T).\tag{10}
\]

One elementary proof of the first inequality is useful for its precise
scope. Coordinate-wise fundamental theorems of calculus and two applications
of Cauchy--Schwarz give, for compact smooth g on R3,
||g||_(3/2)<=product_i ||partial_i g||1^(1/3). Apply this to g=|f|^4.
Since ||partial_i |f|^4||1<=4||f||6^3||partial_i f||2, division by
||f||6^3 (unless f=0) and the arithmetic-geometric mean inequality prove
||f||6<=C||grad f||2. The argument holds for vector f as well.

## 4. Pressure commutator with an explicit radius power

For T_ij=R_iR_j and eta_R=phi_R^4 split

\[
\eta_R\pi_*=\sum_{i,j}T_{ij}(\eta_R g_{ij})
                +\sum_{i,j}[\eta_R,T_{ij}]g_{ij},
\quad[\eta,T]g=\eta Tg-T(\eta g).\tag{11}
\]

The first term has L^(3/2) norm at most C_T(B_R+1):

\[
\|\phi_R^4w_iw_j\|_{3/2}\le B_R W_T,\qquad
\|\phi_R^4w_i u_j\|_{3/2}\le W_T\|u\|_6,
\tag{12}
\]

and the analogous u_i w_j bound. The L6 norm of u is uniformly finite
on this fixed [0,T]. Riesz boundedness then supplies the stated estimate.

The off-diagonal kernel of T_ij is
(3 z_i z_j-delta_ij |z|^2)/(4pi |z|^5). Any local identity distribution
in T_ij cancels in the commutator. The boundedness and Lipschitz bound of
eta_R imply its commutator kernel is bounded in absolute value by

\[
K_R(z)=C|z|^{-3}\min(|z|/R,1).\tag{13}
\]

For p=4/3 the radial integrals are explicit:

\[
\|K_R\|_{4/3}^{4/3}
 \le C\left[R^{-4/3}\int_0^R r^{-2/3}dr
                    +\int_R^\infty r^{-2}dr\right]
 =C(3+1)R^{-1}.
\tag{14}
\]

Hence Young's convolution inequality gives

\[
\left\|\sum[\eta_R,T_{ij}]g_{ij}\right\|_{4/3}
              \le C_T R^{-3/4}.\tag{15}
\]

Both parts of min(r/R,1) matter: without near-field cancellation the
L^(4/3) integral diverges at zero; using only the Lipschitz bound makes
it diverge at infinity. An L1 stress also need not produce an L1 pressure.
For example take g_11 a radial smooth unit-mass density in the unit ball
and all other entries zero. Outside that ball pi_*=(3 omega_1^2-1)/(4pi r^3),
whose absolute integral on 1<r<R is 4/(3sqrt(3))*log R. The proof never
uses a false L1 -> L1 pressure estimate.

To justify (11) for L1 tensors, first truncate g smoothly on growing balls.
The truncations converge in L1 at each time, and in L1 of space--time by
the uniform bound and dominated convergence. Their transforms converge
in H^-2, and (15) gives convergence of their commutators in L1_t L^(4/3)_x.
For each fixed R, eta_R g is unchanged once the truncation contains its
support, so its transformed main term is unchanged too. Thus (11) holds
distributionally for g, and pi_* is locally integrable in space--time.
For fixed R, D_R is bounded in time by compact-set smoothness, which
suffices for this local-integrability conclusion without a global bound.

Test (8) against a(t) chi_R w. Since div w=0 and all test functions are
compactly supported, integration by parts gives equality of the fluxes
int pi w.grad chi_R = int pi_* w.grad chi_R for almost every time.
Time-dependent constant pressure differences therefore do not matter.

## 5. Absorb pressure and transport fluxes

Because |grad chi_R|<=C R^-1 phi_R^7, split (11) leaves the test factor
R^-1 phi_R^3 w. Interpolation gives

\[
\|\phi_R^3w\|_3\le\|\phi_R^2w\|_3
     \le B_R^{1/2}W_T^{1/2},\qquad
\|\phi_R^3w\|_4\le B_R^{3/4}W_T^{1/4}.
\tag{16}
\]

For example phi_R^3 |w|=(phi_R^4 |w|)^(3/4)|w|^(1/4), with
1/4=(3/4)/6+(1/4)/2. Combining (12), (15), (16) gives

\[
\left|\int\pi w\cdot\nabla\chi_R\right|
 \le C_T\left[R^{-1}(B_R^{3/2}+B_R^{1/2})
                   +R^{-7/4}B_R^{3/4}\right].\tag{17}
\]

The eighth power of phi is sufficient and is the smallest power for THIS
L^(4/3) commutator factorization. For chi=phi^k the commutator pairing
needs k/2-1>=3k/8, that is k>=8. This is a constraint of this proof,
not a claim that every possible comparison argument needs eight.

Pair (4) with chi_R w. Smoothness on compact sets justifies every
integration, and the exact identity is

\[
\begin{aligned}
\tfrac12E_R'+D_R^2={}&-\int\chi_R(w\cdot\nabla)u\cdot w
 +\tfrac12\int|w|^2\Delta\chi_R\\
 &+\tfrac12\int|w|^2v\cdot\nabla\chi_R
 +\int\pi w\cdot\nabla\chi_R.\tag{18}
\end{aligned}
\]

Take R larger than the support radius of u, so v=w on supp grad chi_R.
Using phi_R^7<=phi_R^6 and (16), its transport flux is bounded by
C_T R^-1 B_R^(3/2). The Laplacian flux is <=C_T R^-2. The strain term
is <=||grad u||infinity E_R. These estimates use no decay of grad v.

Insert (10) into (17) and the transport estimate. The nonconstant
dissipation powers and their Young remainders are

| Term to absorb | Dissipation power | Radius decay after Young |
|---|---:|---:|
| R^-1 D_R^(3/2) | 3/2 | R^-4 |
| R^-1 D_R^(1/2) | 1/2 | R^(-4/3) |
| R^(-7/4) D_R^(3/4) | 3/4 | R^(-14/5) |

In general R^-b D^a<=delta D^2+C_(a,delta) R^(-2b/(2-a)), for
0<a<2; maximization in D>=0 proves this directly. The offsets from
R^-1 W_T in (10) have powers R^(-5/2), R^(-3/2), R^(-5/2), and the
Laplacian term has R^-2. Thus EVERY remainder is bounded by C_T R^(-4/3)
for R>=1. Choose the finitely many delta's so their sum is <=1/2. This
proves (1). Its constants depend on T, W_T, G_T, u and the fixed cutoff;
they do not depend on R or a derivative norm of v at infinity.

## 6. Comparison on one fixed interval, then conditional exclusion

At t=0 the identical initial velocities give E_R(0)=0. Gronwall applied
to (1) yields, uniformly for 0<=t<=T,

\[
E_R(t)\le 2C_T R^{-4/3}\int_0^t
       \exp\left(2\int_s^t\|\nabla u(a)\|_\infty da\right)ds
       \le K_T R^{-4/3}.\tag{19}
\]

The exponent is finite for this fixed T<1 by smoothness and compact
support of u. On any fixed spatial ball chi_R=1 for all sufficiently large
R. Let R->infinity in (19) to see its integral of |w(t)|^2 is zero.
Continuity gives w=0 everywhere on [0,T]. No monotonicity of chi_R in R
and no a priori total dissipation bound for v is required.

Now ASSUME the earlier conditional field/force construction is valid.
Suppose a global smooth bounded-energy solution v with that same force
and zero datum existed. The general lemma applies for every T<1, so v=u
on [0,1). The earlier complete-field bound along
x_tau=(sqrt(2 Xin tau),0,0) is

\[
|u(x_\tau,1-\tau)|\ge |u_\theta(x_\tau,1-\tau)|
          \ge(e_0/2)\tau^{-A}\longrightarrow\infty.
\tag{20}
\]

These space--time points converge to (0,1). A globally smooth v is bounded
on one compact neighborhood of that point, contradicting (20). This closes
the conditional logical implication. It does not turn the inherited
profile or atlas hypotheses into verified facts and does not establish
an unforced singularity.

## 7. Reproducible diagnostics and the research gate that remains

The tests audit exact rational Holder, weight, kernel and Young exponents;
improper radial integrals; negative-Sobolev Fourier bounds; a manufactured
periodic pressure Poisson problem and instantaneous difference equation;
the complete weighted energy identity including every nonzero flux; and
an independently expanded mode-by-mode commutator, including zero output.
The periodic diagnostics test algebra only. They are NOT whole-space
simulations, numerical proof of a bound at infinity or evaluation of the
actual constructed profile. The continuum comparison proof is Sections 1--6.

Failure controls retain reversed pressure sign, omitted commutator,
premature replacement of v by w before support separation, insufficient
cutoff weight, divergent kernel variants, a non-L1 pressure from L1 data,
and an exact infinite-energy counterexample if hypothesis (2) is removed.
Earlier frozen source/evidence is preserved byte-for-byte.

```sh
python3 research/navier_stokes_cascade/results/global_solution_comparison/test_audit.py -v
python3 research/navier_stokes_cascade/results/global_solution_comparison/audit.py \
  --verify-record research/navier_stokes_cascade/results/global_solution_comparison/evidence.json \
  --output /tmp/global_solution_comparison_recheck.json
```

The next useful checkpoint is now FOUNDATION AUDIT, not another downstream
conditional transfer. Start with the global auxiliary support/atlas input:
construct and check the physical representatives on overlapping dyadic
charts, closed separation of interacting labels, bounded overlap,
coordinate-transition compatibility and derivative bounds after assembly.
Trace them to actual defining maps and cutoffs. The first missing bound
or incompatible overlap must change the stated construction or scope.
In parallel, review the leading profile, common background collar,
all-order inverses and fixed parameter inequalities. An external expert
review or independent formal certificate would add evidence these dependent
algebra checks cannot supply. No unconditional blow-up claim is made here.
