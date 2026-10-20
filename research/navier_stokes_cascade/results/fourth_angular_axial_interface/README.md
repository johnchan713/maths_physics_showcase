# Fourth-angular reference axial bound and the remaining compact input inequality

Parent: `5ad730f469503d4d817e3abcd6d58390d8281152`.
Status: `same-reference-root-fourth-angular-axial-bound-and-explicit-compact-input-gate`.

The second-coefficient axis/collar estimate is now available, but its global
pressure-source bound still needs the finalized leading U0 in C_eta^4 and
the finalized first U1 in C_eta^3. This checkpoint advances the FIRST of
those inputs. It derives fourth-angular bounds for the SAME later reference
angular root, complete strength energy and axial pulse. It then joins the
axis-side estimate to an EXPLICIT compact-modulation/repair interface.
The remaining compact fourth-order inputs are not evaluated here, so a
quantified actual global U0 C4 bound is NOT declared closed.

All norms below are factorial angular norms
`|f|k=sum_(j=0)^k sup_[-1,1]|partial_eta^j f|/j!`, including endpoints.
All actual roots, supports, gauges, width C^-20, selected frequency and
schedules remain unchanged. Historical frozen evidence is never edited.
These are analytical consequences of the retained root/core/region lemmas,
not independent unconditional foundations or a continuum proof by tests.
No complete physical perturbation-pulse, covariance/mean, residual or
blow-up conclusion follows. The later REFERENCE axial pulse discussed here
is part of the leading background, not the full physical pulse construction.

## 1. Fourth angular order does not require a different angular root

Write p=1+eta^2, f=1/p and a=log(2/p). On the real eta interval the poles
at +/-i have distance at least one. The exact partial fractions and log jets
give

`f=[(eta-i)^-1-(eta+i)^-1]/(2i)`,
`a^(j)=(-1)^j (j-1)! [(eta-i)^-j+(eta+i)^-j]`, j>=1.

Thus `sup|f^(j)|/j!<=1`, `sup|a^(j)|/j!<=2/j` and |a|0<1. In particular

\[
                         |f|_4\le5,\qquad |a|_4<31/6<8.       \tag{1}
\]

Use the EXACT angular discrepancy formula in
`global_axial_angular_bound`, section 1. Its eta-independent delta0 has
size <=1/beta, beta=1-lambda>1/2, and the interpolation step has positive
mass int sigma'=1. The Banach exponential estimate therefore gives the
SAME bounds, now in C4,

`|deltaI|4 <18 exp(8)<2^21`, `|zang|4<=2^29 lambda^30`.

Here lambda=exp(-4T)<=exp(-256), T>=64. The constant matrices still have
inverse norm <3 and quadratic norm <32. Their C4-algebra contraction
gives `|cang|4<=2^32 lambda^30=:cstar`, with the SAME radius as the prior
generous C3 estimate. The pointwise accepted root and this C4 root lie in
the same uniqueness ball; `192cstar<1` identifies them. No new root,
frequency, C or smallness condition in infinitely many norms is selected.

## 2. Keep the COMPLETE fourth-angular energy in the strength equation

For interpolation, g_theta=(p/2)^theta, 0<=theta<=2. In addition to the
three historical derivatives, the exact fourth derivative is

\[
 g_\theta^{(4)}=g_\theta\left[
 {12\theta(\theta-1)\over p^2}
 +{48\theta(\theta-1)(\theta-2)\eta^2\over p^3}
 +{16\theta(\theta-1)(\theta-2)(\theta-3)\eta^4\over p^4}\right].\tag{2}
\]

Using |g|<=1, p>=1 and the individual factors gives raw derivative bounds
(1,4,12,56,408), hence |g|4<64. The interpolation energy costs at most
64Tf=64000<=1024T, since Tf=1000 and T>=64.

All subsequent UNEDITED normalized power-reference energy is a positive
constant times p^2. Its C4 norm is exactly 25/4 times its value at eta=1,
not the earlier C3 factor six. The retained value comparison gives
R_unedited<=512T, so that norm costs at most 3200T. This includes the
waiting stages and infinite power tail; no infinite interval is multiplied
by a uniform value. Together these unedited contributions cost at most
4224T. Both angular bump edits together cost at most

`25(1.2cstar+40cstar^2)<=30cstar+1000cstar^2<1`.

Consequently the COMPLETE power-reference energy entering the unchanged
strength equation obeys

\[
                                      |R|_4\le8192T.           \tag{3}
\]

This is the power reference before the heat attachment. It does not assert
that the heat-modified R by itself equals R. The earlier heat compensation
and later heat edit have exact zero net S change as functions of eta and
trial strength, so the COMPLETE strength equation is the same. Their
smoothness justifies taking these four derivatives. Full incoming pressure,
waiting/tail energy and compensation are not reset or truncated at I3/I4.
The angular root and R are independent of the axial trial strength, so
there is no circular differentiation of that strength in obtaining (3).

## 3. All four normalized derivatives of the same positive strength

Keep `ci=ai Amp+bi(eta)`, with eta-independent |ai|<=lambda^20 and
`bi=vi eta(1+eta^2)`. Its C4 norm is still 10/6 times its C1 norm, hence
`|bi|4<=2lambda^20`. The prefix polynomial
eta^2(1+eta^2)^2 has C4 norm 93, and its retained coefficient is <1/4.
Thus the SAME incoming s0 has |s0|4<512T.

Retain the full equation `k Amp^2+B Amp+D=0` from the previous checkpoint,
including both affine end bumps:

`k=K+lambda sum wi ai^2`, `B=2lambda sum wi ai bi`,
`D=-C+lambda s0-lambda exp(-26)R/2+lambda sum wi bi^2`.

The old value bounds still give k<=1, Amp<1.2 and inverse Jacobian <3.
Equation (3) gives

`|B|4<=160lambda^41`, `|D+C|4<=8192Tlambda=:epsilon<2^-20`.

For normalized jets un=Amp^(n)/n!, the exact same-root recursion is

\[
 (2k u_0+B_0)u_n=-k\sum_{j=1}^{n-1}u_j u_{n-j}
                    -\sum_{j=1}^{n}B_j u_{n-j}-D_n.            \tag{4}
\]

Using u0<=2, |Bj|,|Dj|<=epsilon, induction with |uj|<=12epsilon gives
`|un|/epsilon<=9+468(n-1)epsilon<12`, n=1,2,3,4.
This sharper finite-order padding replaces no historical bound. In
particular

`|Amp|4<=1.2+48epsilon<2`, `|ci|4<=4lambda^20`.

The raw fourth derivative in (4) includes all 8k Amp' Amp''' and
6k(Amp'')^2 terms, plus 4B' Amp''', 6B'' Amp'', 4B''' Amp', B''''Amp
and D''''. Neither odd cross terms nor linear-coefficient derivatives are
discarded; no parity of Amp is assumed.

On the original later reference axial pulse, the main R0 and the two
end bumps have pairwise disjoint supports, R0<=11, bump values<=64 and
eb exp[-(.5+lambda)y]<=1. With (1),

\[
                    \sup_{\text{later pulse }X}|U_0(X,\cdot)|_4\le5\cdot22=110.\tag{5}
\]

The other reference taper is k(y)eta, |k|<=4, with norm <=8; the remaining
post-pulse/heat stages have U=0. This closes this LATER reference-region
C4 bound relative to its retained exact roots, not all global leading data.

## 4. The compact fourth-order input cannot be replaced by its old C3 cap

Keep A=C^4096, Ahi=A^256, H and N=1+floor(H^32) unchanged. On the actual
compact modulation/first-repair rectangle the ORIGINAL U has C4 norm <=Ahi.
The modulated axial field is `U_N=U+Bcal(y,eta,Ny)/N`. Define the ACTUAL
inputs

`B4=sup_(y,phi) sum_(j=0)^4 sup_eta |partial_eta^j Bcal|/j!`,
`Z4 >= N sup_eta ||partial_eta^4 z_N||inf/4!`.

The fully transformed target z_N contains ALL five incoming moment errors,
their e normalizations, inverse-lambda row and angular derivatives of those
normalizers. Z4 is NOT an unnormalized raw pressure-moment derivative.
The all-order map theorem gives finite constants at these fixed inputs,
but the existing evaluated loop/target budget stops at C3. Neither B4 nor
Z4 is assigned H or H16 here.

Let `r=2000H^16/N` be the already proved C3 root cap. Each normalized
coefficient jet d0,d1,d2,d3 is <=r. The SAME accepted five-dimensional
root has constant B,Q, ||Q||<=1e5 and inverse pointwise Jacobian <=2000.
Its fourth normalized derivative obeys

`d4<=2000[Z4/N+1e5(2d1d3+d2^2)]<=2000Z4/N+6e8 r^2`.

The two highest unknown terms are in the Jacobian, not counted again in
the lower sum. Sum the four known jet bounds explicitly, conservatively
allowing different vector norm conventions: `|c|4<=4r+d4`.
The original unit-mass width-.1 axial repair bumps have value <=320 and
disjoint supports. The same factor e=e_star f has |e|4<=5A. Thus

\[
 K_{compact,4}\le A_{hi}+{B_4\over N}
           +1600A\left(4r+{2000Z_4\over N}+6\cdot10^8r^2\right).\tag{6}
\]

This preserves the exact repair and all its quadratic terms, not a new
small-root solve in C4. For a trial compact certificate K, (6)'s majorant
is <=K if and only if

\[
 \boxed{B_4+3{,}200{,}000 A Z_4
       \le N[K-A_{hi}-6400Ar-960{,}000{,}000{,}000Ar^2].}        \tag{7}
\]

This is an exact acceptance inequality for THIS sufficient majorant, not
a necessary condition on the true field, a new retained-scale selection
or an assertion that the actual data fail. The actual B4,Z4 and (7) remain
unverified. As an ABSTRACT norm control, n^-3 cos(n eta) sin(phi) has
C_eta^3<=8/3 uniformly but its normalized fourth eta derivative at eta=0,
phi=pi/2 is n/24. This is not an actual loop or PDE family; it only shows
why the old C3 estimate alone cannot certify (7). Qualitative smoothness
does not quantify its value at the selected enormous parameters.

## 5. Join to the full-support pressure bound, without claiming closure

The same axis-side activated envelope K(ell_axis), ell_axis=tc/2<t/20,
is holomorphic in eta on its retained outer strip. Angular Cauchy supplies

`Kaxis4=K(ell_axis) sum_(j=0)^4 delta^-j`, delta=rho/(800S).

Factorials cancel in this norm. This uses angular, NOT radial, analyticity;
the physical activation width C^-20 is unchanged. Combining the axis,
compact region and (5) gives the CONDITIONAL finite expression

`K0,4=max(Kaxis4,Kcompact4,110)`.

If K1,3 is a valid FINALIZED global first-coefficient C3 bound, use
`higher_pressure_source_moments.regular_v_bound` to get v0 C2, v1 C2
and v0 C3, the latter from K0,4. Its full k=1 signed source formula then
gives Q10,Q11 on the ORIGINAL entire [0,Xv]. Both shifted axial-viscosity
moments remain, including the weighted one. No moment ends at Xcut, I3 or
I4, and the reserved Pi1/Pi2 pressure memory remains unchanged.

Next quantify the actual B4,Z4 (including moving phase and normalized
targets), the finalized first U1 C3 and coefficient cutoff jets. Then join
the second inner estimate to the complete second repair for H2,3 and check
retained-scale acceptance. Complete physical pulse derivatives/supports,
full-wave covariance with axial velocity/pressure means, all viscosity,
cutoff and interaction residual terms, and independent inherited foundations
remain obligations. No actual global U0 C4 or full H2,3 closure is claimed.

## Reproduction

Symbolic checks cover partial fractions, log jets, every fourth interpolation
and strength term. Five manufactured positive quadratic roots independently
check all four normalized recurrences and the bounds, including eta=+/-1.
These are diagnostics, not selected enormous fields or interval certificates.
Exact rational controls check every budget, eliminated gate and pressure
interface. Historical records are byte-bound and replayed unchanged.

```
python test_audit.py -v
python audit.py --verify-record evidence.json --output fresh.json
```

The output must be fresh. Passing tests is not a continuum proof.
