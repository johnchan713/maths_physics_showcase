# Independent shared-core pressure calculus and a sharper growth partition

Parent: `0e605bb3b289ed662ea5333d0335baac8302a1ff`.
Status: `shared-core-pressure-calculus-rederived-with-sharper-gap`.

The foundation-route audit separated provider routing from independent
validation.  This checkpoint begins the latter at the first shared datum: the
axis pressure used by both `core_contraction_audit` and
`core_growth_activation`.  It independently rederives the real and complex
pressure estimates, expands the non-pressure part of the axial source, and
improves the quantitative separation used by the growth proof.

The same unedited reference measure, pressure gauge, `P,h,j,sigma`, and order
of parameter choices are retained.  No profile, root, width, coefficient,
force or historical evidence is edited.  The result is an exact calculus
lemma for the core provider; it is not an independent reproof of the whole
Banach map or activation collar, and it is not an unconditional pulse,
residual or blow-up theorem.  Exact rational tests below are failure controls,
not a substitute for the displayed continuum argument.

## 1. One datum, fixed before the core scales

For the unedited schedule, put `y=log(X/XR)`,
`f(eta)=(1+eta^2)^-1`, and
`E_id=c(y) f(eta)^theta(y)`, with `0<=theta<=1`.  The axis pressure is

\[
 \Pi_0(\eta)=-\frac12\int_{-\infty}^{\infty}E_{id}^2\,dy
 =-P^2\int_0^1 f(\eta)^{2\theta}\,d\mu(\theta).             \tag{1}
\]

The positive measure is the pushforward of `c(y)^2 dy/(2P^2)`.  The exact
inner power contributes `(5/2) delta_1`, while the already derived schedule
bound gives `mu([0,1])<5`.  Because `dX/X=dy`, the arbitrary radial scale XR
cancels before `Lambda` or `C` is selected.  Thus the two local providers use
the same datum without a pressure feedback loop:

\[
 \mu\ge(5/2)\delta_1,\qquad \mu([0,1])<5.                  \tag{2}
\]

This independently confirms the shared local input and its parameter order.
The later five-moment theorem is still needed to identify (1) with the
completed leading backward pressure; that theorem is not used to define (1).

## 2. Direct real pressure bounds and the complex contraction input

Differentiating (1) on the real interval gives

\[
 \Pi_0'(\eta)=4P^2\eta\int_0^1
   \frac{\theta}{1+\eta^2}f(\eta)^{2\theta}\,d\mu(\theta). \tag{3}
\]

Hence `eta Pi0'(eta)>=0` and, for `|eta|<=1`,

\[
 |\Pi_0|<5P^2,\qquad |\Pi_0'|<20P^2|\eta|,
 \qquad \eta\Pi_0'\ge
 \frac{10P^2\eta^2}{(1+\eta^2)^3}.                        \tag{4}
\]

The last inequality uses only the atom at `theta=1`.  It is a real-interval
bound; it does not replace the complex derivative estimate in the contraction
norm.  On the contraction rectangle, with `r0=sigma/256`,
`Re(1+eta^2)>1-r0^2`.  The common logarithm branch therefore gives

\[
 |\Pi_0|<\frac{5P^2}{(1-r_0^2)^2}<10P^2.                  \tag{5}
\]

Cauchy's estimate on the half rectangle recovers
`|Pi0'|<20P^2/r0`.  This checks the exact pressure coefficient consumed by
the all-index map; a real-axis supremum alone would not do so.

## 3. Exact pressure source and its uniform atom floor

Let

\[
 A=\tfrac12+h,\quad D=\tfrac12-h,\quad d=1-\eta^2,
 \quad U_*=4\eta+j,\quad H_*=D\eta+dU_*.
\]

The pressure contribution to the fixed axial source is

\[
 -d\Pi_0'+4A\eta\Pi_0=-\eta P^2K_p(\eta),
\]
\[
 K_p=\int_0^1\left(\frac{4d\theta}{1+\eta^2}+4A\right)
            f^{2\theta}\,d\mu.                            \tag{6}
\]

Writing `x=eta^2`, the exact contribution of the required atom is

\[
 K_{atom}(x,A)=10\frac{(1+A)+(A-1)x}{(1+x)^3}.              \tag{7}
\]

After removal of the positive denominator and factor ten, its derivative is

\[
 -4-2A+2(1-A)x\le-2-4A<0
 \quad(0\le x\le1).                                       \tag{8}
\]

Thus the minimum is the endpoint value
`K_atom(1,A)=5A/2>=5/4`; positivity of the rest of the measure gives

\[
                         K_p(\eta)\ge5/4.                  \tag{9}
\]

This proves the lower pressure part directly from the same measure used in
the complex contraction bound.  Neither a fitted atomic diagnostic nor the
later global pressure is substituted for it.

## 4. Sharpen the velocity remainder and the small-source sector

The non-pressure source is independently expanded from its definition:

\[
 R_v=-A(1-2\eta U_*)U_*-4H_*
\]
\[
 =-20\eta+(16+32A)\eta^3-(A+4)j
 +(4+16A)j\eta^2+2Aj^2\eta.                               \tag{10}
\]

For `1/2<=A<=.51`, `j<=.05`, and `|eta|<=1`, the first two terms obey

\[
 |\eta|\,|-20+(16+32A)\eta^2|\le20|\eta|.                 \tag{11}
\]

The constant term is below `5j`; using `eta^2<=|eta|`, the remaining
coefficient of `|eta|` is at most

\[
 (4+16(.51))(.05)+2(.51)(.05)^2=0.61055<1.               \tag{12}
\]

Consequently the inherited bound can be sharpened, without changing any
parameter, to

\[
                         |R_v|\le21|\eta|+5j.              \tag{13}
\]

Let `delta_*=jP^2/100` and
`Z_*=R_v-eta P^2K_p`.  If `|Z_*|<=delta_*`, (9), (13), and `P^2>=256`
give

\[
 \frac{|\eta|}{j}
 \le\frac{1/100+5/256}{5/4-21/256}
 =\frac{189}{7475}<\frac1{39}.                             \tag{14}
\]

This improves the retained `1/30` localization.  Since

\[
 \frac{H_*}{j}=1+(D+4)\frac\eta j-j^2(\eta/j)^2
                    -4j^2(\eta/j)^3,
\]

the worst absolute-term bound in (14) yields

\[
 H_*/j\ge
 1-\frac{9}{2\cdot39}-\frac1{400\cdot39^2}
       -\frac4{400\cdot39^3}>\frac{22}{25}.                \tag{15}
\]

With the unchanged `sigma=j/100`, the angular weight therefore satisfies

\[
 \chi=\frac{H_*^2}{H_*^2+\sigma^2}
 \ge\frac{88^2}{88^2+1}=\frac{7744}{7745}>.99.             \tag{16}
\]

Taking the contrapositive preserves, with more margin, the growth-provider
partition

\[
                 \chi\le.99\quad\Longrightarrow\quad
                 |Z_*|>\delta_*.                          \tag{17}
\]

Equations (14)--(17) are the new quantitative result.  They are uniform over
the stated measure family and all `P>=16`; they do not rely on sampling the
`eta~j` microregion.

## 5. What this establishes and what remains

This checkpoint independently confirms that the contraction's complex
pressure coefficient and the growth theorem's positive atom are properties
of one fixed measure and one gauge.  It also independently expands the
velocity source and strengthens the exact sector-separation margin.  No error
was found in the inherited pressure sign or parameter order; the inherited
`40|eta|+5j` estimate was conservative, not false.

The first foundation item is only partially independently validated.  Still
required are:

1. rederive every nonlinear monomial and Lipschitz constant in the infinite
   coefficient-space fixed-point map;
2. rederive the angular comparison, flat Volterra activation and complete
   collar moments, rather than inheriting their bounds;
3. continue the independent frontier through continuation, matching roots,
   compact modulation/repair, exterior and regional joins;
4. retain the later pressure/moment, weighted-edge, coefficient, source-cycle
   and complete residual obligations until their own continuum arguments are
   independently checked.

Passing the exact controls below confirms algebra, frozen inputs and scope.
It does not prove existence of the infinite fixed point or any Navier--Stokes
singularity.

```sh
python3 research/navier_stokes_cascade/results/core_pressure_gap_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_pressure_gap_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_pressure_gap_recheck/evidence.json \
  --output /tmp/core_pressure_gap_recheck_fresh.json
```
