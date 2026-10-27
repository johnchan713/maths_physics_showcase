# Finalized leading-profile certificate from the defining maps

Parent: `285b8c4251db338c866315b138aa0402884fd25a`.
Status: `finalized-leading-profile-certificate-joined-relative-to-exact-core-and-root-hypotheses`.

The parent reconstructed the fixed-order background coefficient seminorm and
real Borel selector before selecting any scale. It left one exact leaf: prove
that the collar, reserved patches and exterior support belong to one finalized
leading profile. This checkpoint closes that leaf from the frozen defining
maps, relative to their exact core fixed point and exact implicit roots. It
does not replace those hypotheses by sampled roots or numerical tolerances.

There is one correction to the inherited shorthand. A global complex-eta
strip through the heat endpoint is impossible: the exterior heat function has
zero Taylor radius there. The construction requires a holomorphic eta strip
only on the protected first collar. On later compact rectangles it requires
all fixed smooth jets, with one-sided endpoint jets at the heat boundary.

Historical evidence is unchanged. The exact controls freeze support order,
the pressure-density identity and the quantifier structure. They accompany
the continuum derivation and are not a substitute for it.

## 1. Core and activated collar

Let `t=C^-20`, `Xa=4/Lambda`, `y=log(X/Xa)=ts`, and
`a^2=Xa exp(t/4)`. The all-order core checkpoint supplies the regular axis
extensions of `F0=E0/sqrt(2X)`, `U0`, `Pi0`, `M0/X` and `V0/X`, together
with all fixed X derivatives holomorphic in eta on a strip for `X<=Xa`.

For `0<=s<=1/4`, activation uses

\[
 \sigma(s)=e^{-1/s^2}g(s),\qquad
 {1\over\sigma(s)}\int_0^s\sigma(u)u^m b(u)\,du
      =s^{m+3}{\cal V}_m[b](s),
\]

and the exact substitution `v=u^-2-s^-2` gives

\[
 {\cal V}_m[b](s)={1\over2g(s)}\int_0^\infty e^{-v}
 (1+s^2v)^{-(m+3)/2}g\!\left({s\over\sqrt{1+s^2v}}\right)
 b\!\left({s\over\sqrt{1+s^2v}}\right)\,dv.       \tag{1}
\]

Here `2<=g<=6` and `g` is independent of eta. Every fixed s derivative of
the integrand is bounded by `e^-v` times a polynomial in `v`, so dominated
differentiation applies through `s=0`. Every fixed eta derivative acts only
on `b`; Cauchy's estimate on a smaller inherited core strip therefore gives
a positive strip for each fixed mixed order. The strip may shrink with the
order. Converting s derivatives to X derivatives costs at most a fixed power
`t^-r=C^(20r)` at radial order r. This is finite after C and r are fixed and
is not asserted uniform in r.

The activation differences factor by the flat scalar
`alpha=t(1-kappa0)sigma(s)s^3`. A second use of (1), with `m=3`, gives all
five moment differences, including the pressure row,

\[
 (\Delta M,\Delta I,\Delta J,\Delta S,\Delta C_p)
   =t^2(1-\kappa_0)\sigma(s)s^6(b_M,b_I,b_J,b_S,b_C). \tag{2}
\]

Thus no `0/0` assignment occurs at the endpoint, every fixed jet joins
flatly, and the complete pressure primitive is retained. Equations (1)-(2)
extend the core certificate to the full closed collar `0<=X<=a^2`.

## 2. Why the finalized profile keeps that collar

The retained cutoffs obey the strict logarithmic order

\[
 X_{keep}=X_a e^{t/8}<X_{cut}=X_a e^{3t/16}<a^2=X_a e^{t/4}. \tag{3}
\]

The prescribed continuation starts after the first activated collar, and
every compact modulation, exact moment repair, annular root and heat edit is
supported strictly to its right. Those maps cannot alter (1)-(3). On their
own finite closed X rectangles, the continuation estimates, finite-N
modulation/repair estimates and implicit-root estimates give all fixed
smooth jets. The exterior heat endpoint supplies smooth one-sided jets, not
a false complex neighborhood. Consequently these pieces belong to one
completed leading profile while the first collar keeps its holomorphic
fixed-order seminorms.

Pressure is propagated by its forward primitive through every piece. In
particular, the exact fifth-row repairs and backward pressure identification
use the same pressure; no local density equation is used to reset its gauge.

## 3. The exact I3 patch and ordered coefficient bumps

With `rho3=XR exp(241T-12)` and `y3=log(X/rho3)`, the reserved intervals are

\[
 I_3=(0,5),\qquad I_4=(6,11).                           \tag{4}
\]

They lie before the axial pulse. On both intervals the actual leading
profile is

\[
 U_0=0,\quad E_0=e_3(\eta)e^{-(1/2+\lambda)y_3},\quad
 M_0=\mu\eta,\quad V_0=-\mu,                            \tag{5}
\]

where `0<lambda<=1/100` and `e3` is the inherited positive amplitude. Its
full pressure is

\[
 \Pi_0=K_3(\eta)-{e_3(\eta)^2\over2(1+2\lambda)}
                   e^{-(1+2\lambda)y_3}.               \tag{6}
\]

Differentiating (6) gives `(Pi0)_y=E0^2/2` exactly. The function `K3`
contains all incoming pressure, including generally nonzero heat-compensation
memory. It is not set to zero.

The three E repair bumps have centers `1/2,3/2,5/2`; the two U bumps have
centers `1,2`; each has total width `1/10`. Their closed supports strictly
interlace, end at `y3=2.55`, and lie inside I3. I4 is disjoint and remains
reserved for the complete-wave mean correction. Notice that `U0=0` on I3
does not imply `M0=V0=0` there; (5) records the nonzero accumulated flux.

## 4. Compact identities and result

After the completed axial support `Xv`, the exact total-M identity gives

\[
 U_0=M_0=V_0=0.                                        \tag{7}
\]

The global five-moment identities make the leading pressure source compact;
the exterior pressure remains the same forward pressure. The coefficient
support induction uses the five ordered I3 bumps, reconstructs positive-order
moments from the axis, and retains its separate pressure tails. Hence the
leading profile supplied here is exactly the profile used by the fixed-order
background inverse and Borel selector.

This closes the leading-profile leaf of the first current-head background
contract, conditional on the frozen exact core and root hypotheses. It does
not close the second primitive contract. The next obligation is to audit on
this same completed profile the continuum infima for the core, determinant,
normal and inward-cone gaps, with the flat-edge factorizations and one common
selector. Sampled diagnostics cannot supply those infima. Actual-pulse
covariance, complete residual improvement, force, exclusion and blow-up
claims remain open.
