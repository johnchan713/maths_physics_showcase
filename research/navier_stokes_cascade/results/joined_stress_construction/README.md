# Joined-profile stress construction and review

Status: **`joined-stress-proposal-reviewed-derivative-budget-open`**.

This checkpoint constructs a stress-extension proposal for the actual joined
leading profile and reviews the estimates needed to accept it. The review
**does not certify the proposed frequency**. It fixes the support geometry,
derives the missing second-angular-derivative transfer equation, and gives
checks that reject two tempting but invalid shortcuts. The remaining actual
input and pressure-error bounds are recorded explicitly below.

There is no new blow-up result, complete Navier-Stokes solution, or full
admissible-stress certificate. This is the same agent's mathematical and
computational review, not independent peer review or a proof-assistant check.
The earlier leading-profile arguments are dependencies, not re-proved here.

Source: the supplied [Navier-Stokes manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
especially (4.16), (B.26), (B.34), and Appendix C. Its frozen PDF SHA-256 is
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Parent: `95560e73396c4916bb3313d021be9c400b76da41`.

## 1. Actual support locations

Use the parameters already fixed in the
[core attachment](../axis_core_attachment/README.md):

\[
X_a=4/\Lambda,\quad t_1=\kappa_0=C^{-20},\quad
t_c=\frac{t_1}{10\sqrt{\log C}},\quad
X_R=110(CP)^{10},\quad T_w=240T_d.
\]

The radius called `Xb=100` in that note is a local continuation cutoff, not
the final outer stress edge. This note calls it `Xcut` when needed.

| Purpose | Exact location |
|---|---|
| Protected analytic portion | Through \(X_a e^{t_c/4}\) |
| Left modulation endpoint | \(X_-=X_a e^{t_c/2}\) |
| End of known first collar | \(X_a e^{t_c}\) |
| Right modulation endpoint | \(X_+=X_R e^{T_d+3}\) |
| First reserved repair patch | \(X=\rho_p e^y,\ 0<y<5\) |
| Repair scale | \(\rho_p=X_R e^{T_d+T_w-23}\) |
| Actual bump supports | Union of five intervals inside \(.45<y<2.55\) |
| End of compact comparison region | \(X_{\rm rep}=\rho_p e^5\) |

The power stage starts at \(X_R e^{T_d+2}\), so the proposed right endpoint
is one unit into that stage. The gap to the first repair patch is exactly
\(T_w-26>0\) in log radius. The heat compensation occupies a later reserved
patch. The new stress repair uses precisely the first patch from the
[conditional stress note](../stress_realization_audit/README.md).

The inherited logarithmic parameter inequalities imply
\(C^{-1}<X_-<X_{\rm rep}<C^{12}\), and
\(\lambda^{-1},h^{-1}<C\). These comparisons are checked outward at 80 and
110 digits without forming C itself.

**Keep the small offsets separately.** Represent the inner endpoints as the
common base \(X_a\) and the ordered offsets \(t_c/4<t_c/2<t_c\).
Adding those offsets to a rounded \(\log X_a\) can erase their separation.
Even the manageable example \(t_c=e^{-2000}/100\) gives identical values for
`1+tc/4` and `1+tc/2` at 80 digits. The ordering is exact because \(t_c>0\);
the actual construction never substitutes the rounded, coincident endpoints.

## 2. Quantitative margins on the original profile

This paragraph extends the inequalities in the frozen attachment argument;
it does not evaluate the profile on a grid. During activation write
\(e_a=(1-\kappa_0)\sigma(y/t_1)\), \(y=\log(X/X_a)\).
At the selected left endpoint, \(y/t_1=1/(20\sqrt{\log C})\).
For this argument of the flat step,

\[
\sigma(y/t_1)\ge\tfrac12 C^{-400},\qquad
e_a\ge\tfrac14 C^{-400}\ge C^{-401}.
\]

Monotonicity preserves the last bound through activation. The preceding
argument supplies \(P_c-v_s\ge2e_a\) and

\[
R_a:=\frac{(v_s-2)_+J_c^2}{2(P_c-v_s)^2}
\le\exp(-29\log C+200B-\log8)<C^{-28}<\tfrac12.
\]

Consequently, wherever \(v_s\ge2\) in this part of the selected compact
interval, its third cone gap is at least \(C^{-401}\) and its fourth is
at least \(C^{-802}\). Both exceed the conservative floor \(C^{-1024}\).
The first collar also keeps \(v_s>2.29\): the inherited shear loss is below
.001 and the comparison error below one millionth. After activation the
interior continuation has \(v_s<1\), so it will require the loop edit.

A slightly stronger barrier gives room in \(P_c-2\) along the shape stage.
At the local cutoff X=100 the reference \(p_1>3\), with comparison error
below one millionth. Once the axial taper has ended, \(b_s=0\).
In the constant-slope and shape intervals \(S_q>1\), \(l\le.65\), so at
a proposed downward crossing of \(p_1=2.5\),

\[
D_Xp_1=XS_q/L-lp_1>X-1.625>0.
\]

Thus the older barrier at 2 can be strengthened to 2.5 there. This uses
the source and transition-error estimates of the attachment, not just the
fact that the older inequality was strict. The matching annulus already
has \(P_c>3.32\).

On the first outer slope ramp, \(U=4\eta\), \(l\in[0,.6]\), and the source
formula gives \(S_q\ge-h\): the remaining two terms are nonnegative.
Its incoming \(p_1>1.1X_R\) therefore gives

\[
p_1\ge X_R e^{-.6}
\left(1.1-\frac{.01}{.98}(e-1)\right)>.5X_R.
\]

The existing Md64 axial and intermediate certificates apply to the later
stages. In particular their axial transformed gap exceeds .417 and their
intermediate quadratic bracket exceeds .3135 at the sufficient radius floor.
These arguments supply useful margins; they do **not** supply the missing
global derivative envelope in section 6.

## 3. Why the earlier C1 acceptance record is insufficient

The loop depends on \(p_s\). The original formula (4.16) expresses \(p_s\)
using angular derivatives of cumulative moments. Differentiating the loop
with respect to eta therefore requires *second* angular derivatives of
those moments and of the matching coefficients. The frozen annulus evidence
explicitly bounds their C1 norms and establishes smoothness. It does not
instantiate this next derivative bound for the actual incoming target.

This is a missing quantitative estimate, not a counterexample to the joined
profile. For example the manufactured family

\[
g_k(\eta)=\frac{\epsilon}{2k}\sin(k\eta),\quad k\ge2,
\qquad
\|g_k\|_{C^1}\le\epsilon,\quad
\|g_k''\|_\infty=\epsilon k/2
\]

has arbitrarily large second derivative. Every member is analytic. With
\(\epsilon=10^{-16}\), the second derivative is \(5\,10^7\) for
\(k=10^{24}\), despite the same tiny C1 bound. Neither analyticity alone
nor a small C1 acceptance number can replace a quantified C2 estimate.

The exact core has higher-derivative bounds available in its weighted norm.
The next useful calculation is to transfer those through the shape formula,
the five original moment integrals, and their normalized row combinations.
This checkpoint does not assert that the transfer is impossible; it records
precisely where it has not yet been established.

## 4. Constructed second-derivative repair equation

For the *earlier joining annulus*, let

\[
Bc+Q(\eta)[c,c]=z(\eta),\qquad
L_c=B+2Q(\eta)[c,\cdot].
\]

Here z is the negative actual incoming target after axial restoration, B
is independent of eta, and Q includes the small U-square term. A dot on Q
below differentiates only its explicit eta dependence. Twice differentiating
the original equation gives the exact identity

\[
\boxed{L_c c''=z''-2Q[c',c']-4\dot Q[c,c']-\ddot Q[c,c].}
\]

Both the coefficient 4 and the last term are necessary. The explicit angular
dependence is through
\(\zeta=e^{6/5}(1+\eta^2)^2/P^2\). Its factorial-weighted C2 norm is at
most \(20e^{6/5}/256<1\). Its value, first derivative, and unweighted second
derivative are also individually less than one. The same positive bump
estimates as before therefore give
\(\|Q\|,\|\dot Q\|,\|\ddot Q\|\le100\).

Use \(\|B^{-1}\|\le1000\) and the established
\(\|c\|_\infty,\|c'\|_\infty\le r_0=6\,10^{-13}\). Then

\[
\|L_c^{-1}\|\le\frac{1000}{1-200000r_0}<1001,
\]
\[
\boxed{\|c''\|_\infty\le
\frac{1000}{1-200000r_0}(Z_2+700r_0^2)
<1001Z_2+3\,10^{-19},\quad
Z_2=\|z''\|_\infty.}
\]

`curvature_bound` implements this conditional transfer. It requires a finite,
nonnegative Z2; it does not silently take Z2 from the older C1 record. An
alternative is to prove that the actual incoming data meet the original
entry tolerance in the factorial-weighted C2 algebra. The same contraction
then gives a C2 small root, since the new zeta norm still lies below one.
Either route requires transferring the actual core data first.

`review.py` constructs a manufactured vector c(eta), forms its target with
the existing five-moment map, differentiates that target independently, and
recovers c''. Omitting the mixed term or freezing zeta produces a detected
error. Quadrature here tests a differentiation identity; it is not substituted
for exact bump integrals or an actual-profile derivative supremum.

For the needed radial derivatives of the joining bumps, the flat step has
\(\|\sigma'''\|_\infty<10^6\). To see this, put
\(g=2/x^3+2/(1-x)^3\). The logistic derivative identity is

\[
\sigma'''=\sigma(1-\sigma)
[(1-6\sigma+6\sigma^2)g^3+3(1-2\sigma)gg'+g''].
\]

On [1/4,3/4], bound the three terms by their endpoint absolute majorants
and use \(\sigma(1-\sigma)\le1/4\). On (0,1/4], use
\(\sigma\le e^{16/9-1/x^2}\), \(g\le3x^{-3}\),
\(|g'|\le7x^{-4}\), \(g''\le25x^{-5}\). The bracket is bounded by
\(32x^{-9}\), whose product with \(e^{-1/x^2}\) increases on this interval.
Reflection treats the other endpoint. Both resulting bounds are checked
outward. A joining bump of width .08 thus has second radial derivative
below \(10^6/.08^2<2\,10^8\). This explicit width loss must be retained.

## 5. A finite-frequency proposal, with its hypotheses exposed

The proposed input envelope is \(\mathcal A=C^{4096}\). This is a *target
estimate*, not a proved bound. It would have to bound the relevant field and
moment mixed jets through total order two, the loop input first derivatives,
the needed positive reciprocals and normalizations, and the inverses of the
original cone and boundary-collar margins. In particular require
\(a\ge1/\mathcal A\), \(|t_s|,|p_1|,|p_2|,v_s\le\mathcal A\), and
\(P_c(t_s)-2\ge1/\mathcal A\). Require the original third and fourth gaps
on \(v_s\ge2\), and the boundary \(v_s-2\) gaps, to exceed
\(1/\mathcal A\). Always take \(\mathcal A\ge2^{16}\).

Conditional on those input bounds, choose

\[
d_0=(16\mathcal A)^{-1},\quad
Z=16(1+2048\mathcal A^5)^2,\quad
\mu_{\max}=64\mathcal A^{3/2}e^Z,
\qquad \delta_L=e^{-\mathcal A^{20}}.
\]

The uniform-cap lemma of the stress note applies for every pressure
coordinate, including zero. Elementary comparisons give
\(\mu_{\max}<e^{\mathcal A^{16}}\),
\(T:=\sup|t|<e^{2\mathcal A^{16}}\), and
\(J:=\sup|p_2-p_1t|<e^{3\mathcal A^{16}}\).
With \(\gamma=15/(16\mathcal A)\), the chosen delta is below
\(1/2,\gamma/2,\gamma^2/[4(1+J^2)]\) and the boundary margin.
These comparisons use monomial domination on \(\mathcal A\ge2^{16}\), not
floating-point values of the exponentials. The code checks their coefficient
inequalities; their applicability still requires the stated input envelope.

There is also a useful bound for differentiating the variance root. Put
\(F(\mu,p)=\sqrt{V(\mu,p)}\) for \(\mu\ge0\), with its smooth signed
extension across zero. For \(g=\log I_0\), tilted variance gives
\(g''(z)\ge e^{-2|z|}/2\). Hence
\(V_\mu\ge d_0^2\mu e^{-4|\mu p|}\). The preceding excursion bound
\(\sqrt V\le16d_0\mu\) implies

\[
F_\mu\ge\frac{d_0}{32}e^{-4|\mu p|}>0,
\qquad F_\mu(0,p)=d_0/\sqrt2.
\]

There is no division by p or by a vanishing variance here. For
\(q=e^{z\sin\theta}/I_0(z)\), write
\((q-1)/z=\int_0^1q_z(tz)\,dt\).
The bound \(|q_{zz}|\le5e^{2|z|}\) also gives
\(|F_p|\le3d_0\mu^2e^{2|\mu p|}\).
These are tools for a later quantitative derivative estimate; they do not
bound the entire reparametrized actual loop by themselves.

Define the following candidate hierarchy symbolically:

\[
G=e^{\mathcal A^{32}},\quad H=e^{G^4},\quad
\epsilon_* =G^{-64},\quad
C_{\rm state}=D=H^{16},\quad
C_{\rm corr}=\mathcal A^{128},\quad r=\mathcal A^{-128},
\qquad \boxed{N_{\rm proposal}=1+\lfloor H^{32}\rfloor.}
\]

Under the input envelope, the loop gap estimates give a common gap floor
\(G^{-4}\), a positive a floor \(G^{-3}\), and coordinate magnitudes below
\(2\mathcal A\). Thus the old cone Lipschitz argument permits the proposed
epsilon. If the displayed error constants and tolerance were established
for the actual modulation and repair, the proposed integer would exceed

\[
\max\left\{
\frac{2(C_{\rm state}+2C_{\rm corr}\beta D)}{\epsilon_*},
16\beta^2q_*D,\frac{4\beta D}{r},1\right\},
\quad \beta=1000,\ q_*=10^5.
\]

For example, taking logs reduces the first comparison to domination of
\(64\mathcal A^{32}+128\log\mathcal A\) and fixed constants by
\(16e^{4\mathcal A^{32}}\). All quantities are finite, but **making N
enormous does not establish its unproved input or error estimates**. Neither
N nor the actual high-frequency profile has been materialized or certified.
The leading modulation N must also be fixed before the separate physical
wave parameter q of the later PDE construction.

## 6. Review decision and next calculation

The review rejects promotion of this proposal to an actual stress certificate.
Four estimates remain unclosed:

1. Transfer the actual core and shape-transition data to a uniform C2 bound
   on all five normalized incoming moment functions. Retain derivatives of
   the normalization factors. Apply the exact equation in section 4.
2. Bound the actual compact field, moment, and stress-coordinate jets and
   their positive denominators. Justify an envelope such as the proposed
   \(C^{4096}\) with a calculation, not its size.
3. Bound fixed-phase loop derivatives, their inverse phase map, and the
   modulation's effect on the original (4.16) pressure formulas. Prove
   Cstate and D; D must retain the inverse-lambda loss. Accumulate the small
   M/J difference before rounding, or bound it before subtraction.
4. Prove Ccorr and a positive-field correction neighborhood on the first
   reserved patch. Only then apply the existing finite-frequency criterion.

The first of these is the next calculation. The derivative transfer derived
here makes its required output explicit: an actual bound Z2, or the stronger
actual C2 entry bound. The retained failures concern invalid inference and
rounding, not evidence that the manuscript's construction is impossible.

Even successful completion would certify only the leading admissible stress.
The full PDE corrections, convergence estimates, smooth forcing, and actual
finite-time unbounded growth would still need their own arguments.

## Reproduce

From the repository root, with the pinned mpmath dependency installed:

```sh
python research/navier_stokes_cascade/results/joined_stress_construction/test_audit.py -v
python research/navier_stokes_cascade/results/joined_stress_construction/audit.py \
  --verify-record research/navier_stokes_cascade/results/joined_stress_construction/evidence.json \
  --output joined_stress_construction_reproduced.json
```

Passing the audit means that the geometry, conditional transfer constants,
identity controls, frozen provenance, and rejection of unsupported promotion
reproduce. It does not make any of the four open estimates true.
