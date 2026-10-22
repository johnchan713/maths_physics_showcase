# Initial current Haar rows: full weighted equations and a sharper average rate

Parent: `b4505a01b491e15eec5338f9386f35b4964bf68c`.
The parent identifies a weighted CURRENT Haar-row inequality needed by the
first signed target. Here we derive its row polynomials from the complete
initial primary/fast/slow state, rather than assume an independent current
row norm. The result is conditional on the retained exact core/background,
canonical primary and fixed mean-map construction premises. It is not an
independent verification of those premises or the complete residual.

For this INITIAL prefix only, the Haar rows have sufficient weighted order
beta=2-2kappa, stronger than the retained full auxiliary-dependent row order
alpha=2-3kappa. The latter is not corrected upward: its fast graph-viscosity
term remains. No later complex state, nonlinear iteration, physical norm
contraction or blow-up follows from the sharper Haar estimate.

All fields, original roots, radius 1/6057072, width C_setup^-20, profiles,
carriers, supports, pressure gauges and accepted common tail are unchanged.
Historical evidence is preserved. No actual giant parameter is sampled.

## 1. Keep weights before aggregation and through the fixed maps

Fix an output order n and take N=n+3. Use the through-order factorial
coefficient norm in (R,Z,T,y), with all fixed R/U and common-atlas jet
conversion constants included. Work on the whole positive shell. Choose
a finite pole p large enough for ALL primary, source and operator input
jets used at this fixed order. Increasing p changes constants and pole
allowances, not the accepted first band.

The root checkpoint bounds complete primary A and epsilon^kappa B by
sqrt(zeta) times explicit polynomial S and edge losses. BEFORE unweighting,
the exact real covariance is

    W_ab=(epsilon/2)sum[A_a A_b+k^-2 B_a B_b].             (1)

Leibniz, closed support separation, k^-2<=epsilon and normalized Haar
measure give positive polynomials W_N(S) with W in epsilon*zeta*delta^-p.
Use the same coarsest common band and covering costs as the atlas checkpoint.
There is no label/deck factor and no inversion of an edge weight here.

The canonical PRINCIPAL columns match the two selected fluxes exactly.
The original squared slow partitions recover the same physical stress;
no partition is divided out or historical root substituted.
Let Sigma be the full background stress, Sigma0 its original leading part,
and retain its weighted higher-stress bound Sigma-Sigma0 in M_3. Then

    D=(bar W_rtheta-Sigma_theta,bar W_rz-Sigma_z)
       in epsilon^(2-2kappa)*zeta*delta^-p,
    D_N=(1/2)sum_{j<=N}sum_{i<=j} C_i C_(j-i)+S_N        (2)

is a sufficient polynomial after the fixed chart/atlas conversions.
C_i are the weighted normalized curl-coefficient polynomials; S_N bounds
the higher background stress. All other W entries are still only M_1.
Equation (2) uses the REAL primary's exact sine/cosine quadratures. It
does not apply after a generic complex sourced wave.

Here is why the actual initial fast/slow inputs are weighted, rather than
an inference from their unweighted norms. The initial zero-Haar source is
constructed from (1), graph multipliers and T0(gr^circ), as in the initial
scalar checkpoint. Differentiation has the already retained kappa cost;
the source therefore has order H=1-kappa and full zeta weight, with finite
poles and polynomial S norms. The graph/source formulas there are used
BEFORE edge unweighting. The compensator disappears only in this zero-Haar
pressure channel; its moving Haar moment is retained below.

The same shifted T_e chooses the primitive toward each shell edge. At
fixed integration shift its differentiated formula has no M factor. After
taking absolute values, its one-sided integrals obey the parent's kernel
bound even without a zero moment, with the ACTUAL cutoff collar retained.
Let eta_e>0 be a fixed gap on which chi_m is exactly 0 at the inner edge
and exactly 1 at the outer edge. In (2) of the parent replace its proof
collar r by r_e=min(1,L/4,eta_e); denote the resulting constant C^cut_p.
The left/right kernel proof is unchanged, and its middle weight floor is
now exp((a+b)/r_e^2). Using the default r when a cutoff switches sooner
would miss this potentially enormous fixed gap cost. No physical cutoff
or support is changed by choosing this smaller proof collar.
Terms with a differentiated cutoff are supported in the closed interior.
Thus its weighted norm constants can be constructed as

    K^w_T,e,m(p)=L_e,m*(1+C^cut_p),
    K^w_A,e,m(p)=Z_e,m*KA_e,m*Eedge(p).                   (3)

L_e,m is the finite positive sum of differentiated multiplier, cutoff,
Jacobian and fixed chart jet products in the SAME shifted-integral formula;
include both one-sided terms. Z_e,m bounds zeta^-1 on the closed support
of the A_e output (it can be padded by the through-order fixed-profile
bound). KA_e,m is the earlier unweighted compact-error constant. Every
factor in (3) is a fixed profile/hull seminorm or the explicit C^cut_p/Eedge
constant, not a new current-field estimate. This argument requires the
actual edge-directed cutoff and interior A_e support; an arbitrary bounded
unweighted operator would not have this weighted property.

The temporal inverse acts only on y, so it preserves the weight with the
same Fourier loss. For weighted RAW initial source polynomials Et,Ez and
the appropriate finite common pole, sufficient actual fast bounds are

    Ft_m=J_m Gamma_m*6S[Et_(m+4)+(1+K^w_A,1,m)Ez_(m+4)],
    Fr_m=J_m K^w_T,1,m+1 Gamma_(m+1)*6S Ez_(m+5).        (4)

Here Gamma_m=24*2^(m+4), and J_m=binom(m+5,5) safely converts the raw
maximum to the normalized through-order norm in five coefficient variables.
This is not the earlier two-slow-variable count. The input pole p must
cover these higher input jets as well. The real axial
velocity is gd-A1 gd, not gd. Its Haar average is exactly zero, as are the
other fast velocities. The weighted initial source and (3)--(4), rather
than a claim that small unweighted norms imply weighted ones, supply Ft,Fr.
The slow response profiles and their entire potential hull are in I4,
strictly inside the shell. Their explicit scalar-response estimates become
weighted estimates by the fixed positive weight floor on that hull.
The initial scalar checkpoint supplies their CURRENT scalar inputs.

Consequently, for this actual initial prefix write

    mean=(beta_f,v_f,gamma_f)+(B,u,g),
    bar(beta_f,v_f,gamma_f)=0,
    (v_f,gamma_f): epsilon^H zeta delta^-p Ft,
    beta_f: epsilon^(H+1) zeta delta^-p Fr,
    (u,g): epsilon zeta delta^-p (U,Zg),
    B: epsilon^2 zeta delta^-p X.                         (5)

The polynomials U,Zg,X are the fixed slow inverse responses already derived,
with their profile weight-floor constants. No fast/slow cross is discarded
pointwise; only its EXACT Haar average is zero in the next equations.

## 2. Recompute the full pressure and current equations

Let (b,V,G) be the same summed base, with b in epsilon*b_N and V,G
unweighted in V_N,G_N. Put

    Srr=2bB+B^2+bar(beta_f^2)+bar Wrr,
    Srz=b g+G B+B g+bar(beta_f gamma_f)+bar Wrz,
    Stt=2V u+u^2+bar(v_f^2)+bar Wtt.

The exact current Haar radial source is

    gr=epsilon partial_T B-D1 Srr-epsilon partial_Z Srz
           +Stt/R+epsilon(Delta0-R^-2)B,
    P=int gr, pm=T0(gr-rho P).                            (6)

Here D_e=partial_R+e/R and Delta0=partial_R^2+R^-1 partial_R
+epsilon^2 partial_Z^2 AFTER Haar averaging. Before averaging they are
the full graph operators. Haar kills every torus derivative, including
derivatives of the radial graph multiplier; it does not erase quadratic
correlations. The actual full mean pressure has Haar average pm in (6)
because the shifted T0 commutes with normalized Haar averaging.

The two current Haar tangential equations are exactly

    E2=epsilon partial_T u
       +D2[b u+V B+B u+bar(beta_f v_f)+D_rtheta]
       +epsilon partial_Z[G u+V g+g u+bar(gamma_f v_f)+bar Wtheta,z]
       -epsilon(Delta0-R^-2)u,
    E1=epsilon partial_T g
       +D1[b g+G B+B g+bar(beta_f gamma_f)+D_rz]
       +epsilon partial_Z[2G g+g^2+bar(gamma_f^2)+bar Wzz+pm]
       -epsilon Delta0 g.                                (7)

Every complete covariance entry and the real axial and radial velocities
remain. The fast linear time/viscosity terms vanish exactly AFTER averaging;
their full graph-viscosity row still has order 2-3kappa before averaging.
The pressure (6) is reconstructed from the same full state, not from its
leading piece. Pi1/Pi2 are retained in the same base subtraction. Actual
cutoff/A0/A1 velocities remain inside all products. Their separate additive
flat residual ledger is retained by exact subtraction, not assigned a shell
weight without proof.

## 3. Explicit sufficient row polynomials

All state/base norms below are at N=n+3; they include the finite order
needed to differentiate the radial source and then the whole pressure in
(7). Let r_N bound R^-1 in the unweighted through-order norm. Products
are submultiplicative in the factorial norm. Two weighted factors have
weight zeta^2<=zeta and pole 2p. A single weighted factor is padded to
the same pole. Define d=N+1 and

    Lvis=2d^2+r_N d+r_N^2,
    S=2b_N X+X^2+Fr^2+W_N,
    Hrad=b_N Zg+G_N X+X Zg+Fr Ft+W_N,
    T=2V_N U+U^2+Ft^2+W_N,
    Gpoly=d X+(d+r_N)S+d Hrad+r_N T+Lvis X.              (8)

Dropping only stronger positive epsilon powers bounds gr through order n+1 in
epsilon*zeta*delta^-2p with polynomial Gpoly. In particular the pressure
order is 1, not 2: 2V u/R and the other principal variances remain.
Let PiP bound the actual postslow P through order n+1 by epsilon^2 PiP, and let
Rho_N be the weighted bump norm. The exact compensated source in (6)
has zero total. The parent's two-edge primitive and chart/multiplier
coefficients give a finite positive-polynomial Kp_N(2p), hence

    Pm=Kp_N(2p)[Gpoly+Rho_N PiP],
    pm through n+1 in epsilon*zeta*delta^-2p Pm.          (9)

No independent pressure-profile norm is introduced. A monotonically padded
Kp_N overestimates the order-(n+1) operator constant; it does not require
gr jets through N. Kp_N is obtained
from C_(2p), (4) of the parent and its finite chart composition with e=0,
including the R/2 Jacobian; it is a fixed operator constant. The scalar
P is the CURRENT full defect, not a leading stress surrogate.

Set

    Rt=b_N U+V_N X+X U+Fr Ft+D_N,
    At=G_N U+V_N Zg+Zg U+Ft^2+W_N,
    Rz=b_N Zg+G_N X+X Zg+Fr Ft+D_N,
    Az=2G_N Zg+Zg^2+Ft^2+W_N+Pm,
    Epoly2=d U+(d+2r_N)Rt+d At+Lvis U,
    Epoly1=d Zg+(d+r_N)Rz+d Az+Lvis Zg.                 (10)

Termwise use of (2),(5)--(7) proves

    |E_e|_n <= epsilon^(2-2kappa) zeta delta^-2p Epoly_e(S). (11)

The limiting terms are the selected curl-square stresses in (2). The other
radial fluxes have exponents at least 2, the fast radial/tangential products
at least 3-2kappa, and the axial fluxes, time and slow viscosity at least 2.
The epsilon in front of the WHOLE pressure derivative gives order 2 from
(9). No epsilon gain is assigned to pressure itself.

Fixed pullback multiplication by R^(e+1)/2 and the nonlinear chart-jet
composition give the parent's Ehat_e,n from Epoly_e. Since 2-2kappa>
2-3kappa, (11) supplies its displayed inequality (8), CONDITIONALLY on
the lower-level canonical/background/map inputs used in sections 1--3.
There is no longer an independent current-row hypothesis for this INITIAL
prefix. The parent's scalar projection, two-edge primitive, normalized
target and original-root signed bounds may now use these polynomials;
their full moving pressure moment remains.

This does not prove every retained primary/ODE/core premise or replace
the same-domain weighted map construction by numerical evidence. The new
row polynomials can be huge. Their finite values, derivative-order-dependent
poles and fixed profile constants are not uniform in order, width, radius,
lambda or stage. The initial zero-Haar/real-primary reductions must not be
used after arbitrary later complex wave or nonzero old-Haar updates.

Next verify the actual nonzero-angular source through the signed prefix
and the first sourced inverse, with its full-path supports and complete
linear, viscous, cutoff and nonlinear terms. Then join the full pointwise
residual. Infinite realization and unconditional blow-up remain unproved.

## 4. Evidence

Independent exact current-field calculations compare (6)--(7) to the full
current conservative Haar equations before splitting fast/slow fields,
including fast squares, nonzero slow means,
pressure reconstruction, both scalar moments and complete covariance.
Negative controls detect omitted correlations, radial velocity, pressure
and inappropriate full-row promotion. Exact exponent half-lines and positive
polynomial substitutions check (8)--(11). The inherited edge and Fourier
maps are tested separately. These manufactured finite fixtures and passing
tests do not prove continuum suprema or evaluate selected gigantic fields.

```sh
python3 research/navier_stokes_cascade/results/core_initial_haar_row_budget/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_initial_haar_row_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_initial_haar_row_budget/evidence.json \
  --output /tmp/core_initial_haar_row_budget_fresh.json
```
