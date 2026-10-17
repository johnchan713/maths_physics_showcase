# Prescribed continuation through shape: all fixed jets, not a copied shear

Parent: `bde6e3d9d212a706b54d8688941379400f60f79d`.
This checkpoint independently checks the defining maps AFTER the first
activated collar. It proves finite all-fixed-order bounds for the
PRE-MODULATION continuation through Xsep=110 exp(2000B), relative to the
same exact infinite pressure-family core. The reference taper, activation,
two final widths C^-20, shape length and all incoming traces are unchanged.
Historical evidence and source files are not edited.

This is an analytical regularity estimate, NOT a new cone estimate or an
unconditional actual-pulse/residual/blow-up result. Later compact modulation
changes this pre-modulation profile even inside the continuation region;
therefore its regularity is not silently attributed to the FINAL profile.
The first collar remains the actual finalized collar because later edits
start strictly to its right. Global backward pressure identification still
requires every exact fifth-row repair, not just the local density equation.

Defining sources: (4.16), (B.22), (B.26), (B.34) of the retained manuscript
[PDF](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`;
[reduced continuation](../reduced_continuation_margins/README.md) sections
1--4 and [shape/pressure](../shape_join_pressure/README.md) sections 1--2.
Their finite-order margins are not used as a substitute for all-order jets.
The manuscript's final theorem is not an input.

## 1. Which shear is prescribed after the native reference has flattened?

Let Xa=4/Lambda, y=log(X/Xa), t1=kappa0=C^-20. Write
F=E/sqrt(2X), d=1-eta^2, L=1-2h eta^2, Aphys=1/2+h,
D=1/2-h. On 0<=y<=2t1 the REFERENCE is defined by

\[
 D_y\log F_r=\psi D_y\log F_{nat},\qquad
 D_yU_r=\psi D_yU_{nat},                         \tag{1}
\]

where psi is one up to t1, tapers flatly to zero by 2t1, and is zero
afterwards. The native input is needed only on 4<=Y<=4 exp(2t1)<4.1:
C>=2 implies 2t1<1/50, and 4 exp(1/50)<=200/49<4.1.
Thereafter F_r,U_r are radially constant, but E_r=sqrt(2X) F_r is NOT.
Carry the ORIGINAL core traces of the five reference moments, and compute

\[
\begin{split}
 W_r&=1-(2D\eta M_r+dM_{r\eta})/X,\\
 Q_r&=-W_r+[(1-h)I_r-D\eta I_{r\eta}-dJ_{r\eta}
                   +2(h-D)\eta J_r]/(X\sqrt{2X}E_r),\\
 N_r&=-W_rU_r+[D(M_r-\eta M_{r\eta})+4h\eta S_r-dS_{r\eta}]/X
             +4Aphys\eta\Pi_r-d\Pi_{r\eta},\\
 p_{1,r}&=XQ_r/L,\qquad n_{s,r}=N_r/L.           \tag{2}
\end{split}
\]

These are explicit integrals and algebraic formulas for the KNOWN reference,
not a coupled equation for the later actual fields. In particular,

**after y=2t1, reference physical a_r=1-2D_y log E_r=0 and U_{r,y}=0,
whereas p1,r is positive and n_s,r need not vanish.**

The first-collar identities F_a=F_r exp(J), U_a=U_r-K in the parent
use a_r=p1,r and U_{r,y}=-X n_s,r/2, valid while reference=native.
Extending those identities past the taper would prescribe the WRONG field.
This is a domain restriction on their use, not a change to the parent's
valid first-collar statement or to the inherited continuation definitions.

Use instead the actual prescribed slopes, from their original core trace:

\[
 \log E_a=\log E_c(Xa)+\tfrac12\int_0^y(1-a(v))\,dv,
 \qquad U_a=U_c(Xa)+\int_0^y b(v)\,dv.             \tag{3}
\]

Here a=kappa p1,r and b=-kappa X n_s,r/2 through X=100;
kappa=1-(1-kappa0)sigma(y/t1) until t1 and then kappa=kappa0.
On the next width-t1 log interval multiply b by beta, flatly decreasing
one to zero. On the following width-t1 interval keep b=0 and set
a=(1-theta)kappa0 p1,r+(4/5)theta, theta increasing flatly zero to one.
Then a=4/5,b=0 up to Xi=110. There is NO inverse kappa in (3), and
no inverse p1,r is needed to define these fields. At the internal endpoints
the flat steps have all the constant-side jets, so (3) joins smoothly.
At y=0 the parent supplies the flat native/activated join.

## 2. A finite seminorm estimate, with every narrow cutoff cost retained

Let delta_c=rho/(800S) be the parent's common collar strip radius and
delta=delta_c/16. Use D_s={dist(eta,[-1,1])<s}. The following construction
uses nested domains 8delta,4delta,2delta,delta; all lie strictly within
the parent's domain. No numerical lower bound for the actual huge field
is inferred from manufactured fixtures.

For a real log interval of length T and eta in D_s define

\[
 |H|_{r;s}=\sum_{i=0}^r {\sup|D_y^iH|\over i!}.             \tag{4}
\]

These seminorms are product algebras. For a forward primitive
V=V0+int H, and a holomorphic H on the indicated outer strip,

\[
 |V|_{r;s}\le |V0|_{0;s}+(T+1)|H|_{r;s},\quad
 |e^H|_{r;s}\le e^{|H|_{r;s}},\quad
 |\partial_\eta H|_{r;s}\le s^{-1}|H|_{r;2s}.              \tag{5}
\]

The primitive bound uses V_y=H; the factorial at order i>=1 supplies
1/i, while the value costs T. The exponential bound follows by its
power series in this algebra, not an order-independent chain-rule guess.
For any width t>0, the actual flat step contributes the FINITE number

\[
 |\sigma((y-y0)/t)|_{r}\le
          \sum_{i=0}^r c_i/(i!t^i),\qquad
 c_i=\sup_{\mathbb R}|\sigma^{(i)}|<\infty.                \tag{6}
\]

For i=0 take c0=1. Smooth flat extension makes each c_i finite.
For t=C^-20, the ith term costs C^(20i); it is NOT set to one.
These inequalities are a computable seminorm recursion once finite seed
norms and step constants are supplied, not a practical numerical evaluation
at the selected enormous parameters. No uniform-in-C bound is claimed.

Apply (5)--(6), in this order:

1. Native log F, U on 0<=y<=2t1 and 8delta have every fixed radial
   seminorm finite by the parent's zero-free analytic core. Multiply their
   slopes by psi, extend the resulting integrands by zero, and integrate
   (1) to log(110/Xa). Obtain log F_r,U_r,E_r and E_r^-1 on 8delta.
   The logarithm is the continued core branch, not a new pointwise branch.
   E_r^-1 follows from exp(-log E_r), including the positive factor X.
2. Forward-integrate all five reference moments on 8delta. One Cauchy
   loss 8delta to 4delta supplies every eta derivative in (2).
   Inverses of X,sqrt(2X) and L have finite seminorms; |L|>.97 on
   these strips. Thus p1,r,n_s,r have every fixed radial seminorm on
   4delta. Their definition takes ONE extra seed/moment eta derivative,
   not an unbounded differentiation operator on the original space.
3. Multiply these KNOWN functions by the actual kappa,beta,theta and
   integrate (3). Exponentiation gives nonvanishing E_a and E_a^-1 on
   4delta. It does not require a small complex perturbation or a complex
   positivity assertion for p1,r. Fixed-width derivatives may be huge.
4. Compute the ACTUAL moments from (7) below, not reference moments.
   They are holomorphic on 4delta. Actual pressure coordinates from
   (2), with actual fields/moments, are bounded on 2delta. A further
   Cauchy estimate gives every fixed eta derivative on delta.

This finite sequence has no differential equation that feeds eta derivatives
back into its unknown. Arbitrary radial differentiation changes only the
seminorms in (5)--(6), not their strips. In particular the same delta works
for EVERY fixed radial order. Constants may depend on all frozen parameters,
the requested orders, Xa^-1 and the lengths; neither order-uniform nor
width-uniform constants are asserted.

## 3. Full forward pressure and the shape join

At every stage keep the five original incoming traces and the same
Pi_axis(eta). In log radius the ACTUAL density table is

\[
 (M,I,J,S,Cp)_y=
 (XU,\sqrt2 X^{3/2}E,\sqrt2 X^{3/2}UE,
          X(U^2-E^2/2),E^2/2),\qquad \Pi=\Pi_{axis}+Cp.    \tag{7}
\]

In particular Pi_y=E^2/2 uses the actual E. Reference memory is not
silently reused after a change to E,U. The eta-dependent axis pressure
enters N_r in (2) and hence the later U_a in (3); erasing it can change
the field itself, not merely an irrelevant pressure constant.

Set ell_i=log C+log E_a(Xi), Gi=U_a(Xi), z=log(X/Xi),
s=sigma(z/Tsh), Tsh=2000B. On 0<=z<=Tsh prescribe the same shape map

\[
 \log E=-\log C+z/10+(1-s)\ell_i+s\log f,\quad
 U=Gi,\qquad f=(1+\eta^2)^{-1}.                           \tag{8}
\]

log f has its real continued branch on 4delta, away from eta=+-i.
Equations (5)--(6) prove all fixed seminorm bounds for (8) and its
reciprocal on 4delta. At z=0 flatness of s gives log E=log E_a(Xi)+z/10
to every order: this matches the preceding a=.8, since
D_y log E=(1-a)/2=.1. The value .6 is ell=.5+D_y log E,
NOT D_y log E. At z=Tsh, every jet matches
E=C^-1 f exp(z/10), U=Gi, i.e. the unedited ideal-power continuation
before axial restoration. The equality with P f (X/XR)^.1 uses exactly
XR=110(CP)^10 and Xsep=110 exp(Tsh); no radius is moved.

Moment and pressure values at Xi and Xsep are inherited through (7),
not reset to ideal values. Smoothness of the fields does not imply five
moment matching there. The later annular join and compact repair must
establish those identities for their OWN exact roots.

Consequently, for every fixed r,k there is a finite B_rk, on the same
delta, with

\[
 \sup_{0\le X\le Xsep,\;\eta\in D_\delta}
 \sum_{H\in\{F,U,\Pi,M/X,V/X\}}|\partial_X^r\partial_\eta^kH|
       \le B_{rk}<\infty.                                \tag{9}
\]

Near the axis this is the parent's actual collar theorem. Away from it use
partial_X^r=X^-r product_{j=0}^{r-1}(D_y-j), retaining lower log derivatives.
M/X has its axis average extension; V/X is the parent's regular formula,
which takes one more eta derivative and divides only by L. Compactness
at the finite interfaces and flat matching join these estimates. This is
a relative theorem for the PRE-MODULATION leading profile, not for all
subsequent changes. The parent's exact Cartesian similarity operators
then give its leading velocity/pressure bounds uniform over 0<q<=1 on
this bounded normalized-radius region. Constants include Xsep and the
fixed width costs; no tail estimate or full pulse bound is obtained.

## 4. Meaningful checks and the exact remaining gate

Exact rational two-variable jets implement (2), (3), (7) independently of
the old raw second-jet pressure implementation. Four eta centers, three
kappa values (including zero), and jets through total order four check
the prescribed slopes, five moment flows and the full pressure's distinct
axis datum. First output stress jets are compared to the separate old
implementation, not its source strings. A manufactured reference frozen
in F,U has nonzero p1,r: extending the first-collar shortcut gives zero
shear while (3) does not. Dropping the eta-dependent pressure produces the
exact missing term 4Aphys eta c-d c_eta in N_r. These fixtures have smooth
axis-consistent reference primitives, but are NOT the actual core and do
not certify continuum cones or a vanishing residual.

Other exact controls exercise the factorial Volterra/product bounds,
flat endpoint jets, the .1 versus .6 slope distinction and the amplitude
identity in (8). Evidence freezes rational identities/counts and booleans;
floating-point comparison errors are tested with explicit tolerances,
not refrozen as historical near-zero numbers.

The next foundation gate is all-fixed-order smoothness of the SAME
finite-N compact modulation, its phase inverse and its exact five-row
coefficient repair, followed by the annular/exterior joins, complete
global moments, weighted stress edges and reserved supports. A value
state error or a C1 coefficient radius is not a bound for all higher
derivatives. With (9), that next step has finite all-order unmodulated
inputs, but its constants and compatibility still have to be proved.
The downstream conditional pulse/mean/residual/comparison chain is not
promoted to an unconditional construction or blow-up result.

```sh
python3 research/navier_stokes_cascade/results/prescribed_continuation_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/prescribed_continuation_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/prescribed_continuation_foundation/evidence.json \
  --output /tmp/prescribed_continuation_foundation_fresh.json
```
