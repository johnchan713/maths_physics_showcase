# Q80: support-saturated coefficient tuple on the full physical hull

Parent: `30db596f919096bc97a30be0a7cf4c81f43b04c9` (Q79).

**Roadmap correction.** Q79's request to derive the finalized coefficient-one
tuple overlooked `core_background_radial_reconstruction`, sections 1--4.
That checkpoint already composes the full-support first linear repair with
all four tuple first jets. `core_borel_cutoff_certificate` also already
proves K2<=4928 for the identified implemented Borel cutoff, and
`core_leading_background_borel_join` identifies the rule-defined infinite
schedule. These remain conditional on their exact same-field primitive
contracts. Repeating their constructions would not be a new checkpoint.
Historical records and Q79 are preserved; this correction belongs here.

The new quantitative join clips the POSITIVE-ORDER tuple estimate to its
actual common support before applying the full physical-chart comparison.
It avoids multiplying its primitive bound by the enormous exterior hull
endpoint when the primitive is already identically zero there. It does not
clip the leading field, pressure, or radial integration in the repair.

## 1. Exact cancellation precedes the support estimate

For each finalized n>=1, let U_n,E_n vanish for X>=V, with
V=Xvel=rho3 exp(3), and retain the exact first repaired moment

    integral_0^V U_n(X,eta)dX=0, for every |eta|<=1.       (1)

The established flat cutoff and repair bumps give smooth zero extensions.
Differentiating (1) is legitimate on this fixed compact interval at every
fixed available angular order. Consequently the AXIS primitive
M_n=int_0^X U_n and all its angular derivatives vanish for X>=V.
With c=D+2nh, L=1-2h eta^2 and d=1-eta^2, retain

    V_n=[2eta X U_n-2c eta M_n-d M_n,eta]/L,
    W_n=V_n/sqrt(2X), J_n=2eta M_n/[L sqrt(2X)].          (2)

Thus the entire tuple (E_n,U_n,W_n,J_n), including all its fixed jets,
vanishes there. This is stronger than local zero velocity alone, which
would not remove an incoming M_n. The five-row repair used to obtain (1)
still includes pressure source integrals through the ORIGINAL Xv and the
axis-viscosity boundary. Nothing is truncated at V when forming that repair.

## 2. A clipped native recipe, with the full chart retained

The original band chart has r<=R<=H, 1/4<=s<=4, X=R^2/(2s).
Set a=r^2/8 and b=2H^2. These remain the full chart limits. If a>=V,
the positive-order tuple bound is ZERO. Otherwise set b*=min(b,V).
Use the existing reconstruction on [a,b*], with

    Rlo=r/2, Rhi=2H, Xhi=b*.                            (3)

The upper radius Rhi is conservative even after clipping. Its preceding
primitive interval is [0,b*]. The historical global finalized F,U rows,
obtained from the local cutoff product plus the inverse-800 repair, apply
on that interval by their established zero extensions. Therefore its
through-order result supplies the tuple bound on the ENTIRE [a,b]: use
the reconstructed bound up to b*, and exact zero afterwards. Endpoint jets
agree by smoothness. No extrapolation of a local collar estimate is used.

In particular, in the reconstruction's normalized Taylor table the
zero-radial-order entries become

    |partial_eta^k M_n| <= b* k! u_0,                   (4)

instead of b k! u_0; positive radial entries still use
partial_X^j M_n=partial_X^(j-1) U_n. The first W angular derivative needs
U_n through angular order TWO; the existing first C_eta^3 row suffices.
Neither that extra derivative nor its factorial conversion is removed.
The recipe works for any supplied fixed-order rows accepted by the inherited
low-order reconstruction; it claims no new all-order coefficient enclosure.

The leading n=0 tuple must still use the FULL native image [a,b] and its
actual incoming moments. In particular its reserved plateau has M0=mu eta,
V0=-mu, so the zero-extension reasoning above would be false there.
The adapter rejects n=0 and refuses use without the exact total-moment
premise. Its Boolean guard records a required theorem, not a machine proof
that a caller's particular field satisfies (1).

## 3. Join to Q79 without changing the Borel sequence

Let M11* be the clipped coefficient-one result. For the same identified
Borel chi, the already proved Euler bound gives the sufficient first head

    C11*=max(1,9856(1+2h) M11*), Bcorr1=C11*+1/2.       (5)

The factor 9856 is 2 K2, not a radial coefficient-cutoff norm. The half-tail
still requires the retained normalized infinite schedule, as in Q79. Using
a sharper valid upper bound in the finite head does not select new scales.
No claim is made that an arbitrary historical numerical scale sequence was
chosen from this new majorant. The rule-defined same-object schedule join
and its original physical requirements remain the applicable theorem.

Insert (5) into Q79's existing full-chart composition. The leading second
jet, geometry and spectral data remain unchanged. For fixed input rows,
the inherited positive reconstruction is monotone in Xhi: replacing b by
b* cannot increase M11. The improvement can be strict when b>V, and is
exactly zero when the whole chart lies outside support. This removes a
spurious exterior-length cost, not a power of physical epsilon.

There is also an exact stronger conclusion on charts with a>=V: EVERY
positive-order velocity coefficient and the radial potential-cutoff
commutator vanish, since the latter is proportional to J_n. Local finiteness
of the same Borel sum then gives full velocity equal to leading velocity
there, even within Borel cutoff transitions. On a chart crossing X=V,
use the smooth zero extension and (5), not a discontinuous physical cutoff.
No new boundary commutator is introduced by partitioning an estimate.

## 4. The pressure exception is not a tuple component

The established support theorem has first-order pressure extending to Xv,
second-order pressure to Xv, and first-order stress to Xb. In particular on
I4 the two surviving orders remain

    Pi1=K1(eta)-mu^2/(4X), Pi2=K2(eta).                 (6)

Here K2(eta) denotes a pressure trace, not the Euler norm in (5).
Full velocity equality on I4 does NOT imply pressure equality. The radial
variation of the normalized pressure error between Xa<Xb in I4 is
chi(a1 q) q^(2h) mu^2 (Xa^-1-Xb^-1)/4. This is generally nonzero.
Moving cutoff derivatives and both traces remain in mean and residual
calculations. Neither (3) nor (5) is a pressure-support improvement.

## 5. Remaining frontier and validation

This supplies a support-aware coverage adapter for the already established
first-tuple construction and a sharper Q79 first-head input. It does not
independently prove the retained Volterra/core/profile contracts, leading
spectral gaps, actual column/root margins or higher n+6 correction budgets.
The positive-order common velocity support and exact first moment are
inherited analytical hypotheses, not inferred from finite controls.

Next join the leading spectral/cone and canonical column bounds on the
same accepted physical hull; use the original same-field selector and
retain both exceptional pressure orders in the subsequent complete means.
Uniform physical pulse estimates, radial pressure tails, full viscosity/
cutoff/interaction residual improvement and blow-up remain unproved.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`. Exact manufactured controls compare the
clipped and full recipes, check the zero-primitive mechanism and its failure
without the total moment, and preserve the reserved-pressure difference.
Tests and provenance replay are not continuum proofs.
