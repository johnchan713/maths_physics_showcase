# Q84: explicit nonlinear radial remainder envelope

Parent: `2b630e8cf4e1986e906135d87f3964fce800a514` (Q83).
Q82 isolated two original principal products from the final radial source.
Here their remaining source is expressed as one positive norm polynomial,
with an explicit primitive ceiling. It is not a new numerical enclosure of
those primitive norms. We retain the final angular-row identity, all six
wave covariance entries and the complete pressure construction.

## 1. Define the subtraction before taking norms

Use precisely Q82's FINAL mean state (beta,v,gamma), the full base (b,V,G),
and the original principal tensor P. Write

    Delta W = W_final − epsilon P/2,
    Jd = b v + V beta + beta v + Delta W_rtheta,
    Ltheta = G v + V gamma + gamma v + W_theta,z,
    Fsrc = −Etheta − epsilon partial_T v + Dz Ltheta
           − epsilon(D_r D1 v + Dz² v),
    Hrz = b gamma + G beta + beta gamma + W_rz.

With Q the Haar complement and I=c_i^-1 Nt^-1 on zero-Haar data, define

    Rem = −a1 I Q Jd − a I Q Fsrc + v²/R − tstar beta
          −(2b beta+beta²)/R − Dz Hrz
          +(Delta W_thetatheta−Delta W_rr)/R + epsilon Dz² beta,  (1)
    a=2V/R, a1=6V/R²−2V_R/R.

This is a definition using the actual subtraction Delta W, not a smallness
inferred from the two separate order-one tensors. The angular equation in
`core_angular_pressure_normal_form` gives

    A0 Q gr_final = (epsilon/2) A0[Q(Ptt−Prr)/R − a1 I Q Prt]
                   + A0 Q Rem.                                  (2)

The characteristic derivative removed in (2) is still in the local pressure
K1. The pressure primitive, moving rho Pmean, axial A1 velocity and both
exceptional reserved orders are untouched. Etheta is the already derived
FINAL angular residual; it contains no current radial pressure and supplies
no circular estimate for the radial row being bounded. Old dependencies in
its generated-source graph are retained.

## 2. One explicit finite norm polynomial

Use the factorial THROUGH-order norm (sum of absolute normalized derivatives)
in the fixed five coefficient variables, with component maxima. Products
are submultiplicative. Derivatives satisfy |partial_j f|_j0 <=
(j0+1)|f|_(j0+1). All polynomials below are positive and use stripped norm
bounds at the indicated epsilon orders. For output m put N=m+6 and d=N+1.
Let through-N data b0,V,G,x,u,z,W,Z bound respectively

    b/epsilon, V, G, beta/epsilon^(2−kappa),
    v/epsilon^(1−kappa), gamma/epsilon^(1−kappa),
    W/epsilon, (W−epsilon P/2)/epsilon^(2−2kappa).                  (3)

Let e bound Etheta/epsilon^(2−2kappa) through m+4. Let r bound R^-1,
g bound dR^d/dR, and c bound c_i, including all needed coefficient jets
through N where applicable. The radial operator is
D_r=g(partial_U+M v_r.partial_y), with |v_r|_1<2 and M<=epsilon^-kappa.
Put

    Dg=3 d (1+g), A=2Vr, A1=6Vr²+2d g Vr,
    Vis=Dg²+rDg+r²+d²,
    Jd*=b0 u+Vx+xu+Z,
    Lt*=Gu+Vz+zu+W,
    F*=e+d u+d Lt*+Vis u,
    H*=b0 z+Gx+xz+W,
    R*=d(1+4c)x+r(2b0 x+x²)+d H*+2rZ+d²x,
    Rem*=A1 T Jd*+A T F*+r u²+R*.                              (4)

T is the ORIGINAL temporal inverse bound INCLUDING its zero-Haar
projection and four-jet, coefficient, floor and profile conversions:
|I Q f|_m <= T_m |f|_(m+4). It is supplied by the inherited
`temporal_recipe(m)`, not a guessed norm-one inverse. Formula (4) may use
a padded T at the common requested pole. The outer Q costs another factor
two. For native unweighted norms, (4) yields

    |Q Rem|_m <= 2 epsilon^(2−3kappa) Rem*.                       (5)

For the original weighted contracts use the same common archived pole and
fixed-profile conversions BEFORE applying (4); include their constants in
the input bounds and T. Multipliers of different pole indices cannot be
combined merely by renaming their norms. This checkpoint introduces no new
weighted inverse theorem and asserts no universal pole-independent T.

Here is why each term has the required power. Jd has order 2−2kappa,
the non-Fsrc part of Rem has at least that order, and Fsrc has order
2−3kappa. The latter minimum comes from epsilon D_r D1 v: both radial
viscous derivative losses are retained. Lt has order 1−kappa; Hrz has
order one. The positive Dg bounds a radial graph derivative after taking
out epsilon^-kappa. Repeating the estimate at j and j+1 bounds D_r²;
coefficient derivatives of g are included by the product norm. Nt costs
at most 4d, accounting for tstar beta. In V_R inside a1 there is only the
slow radial derivative g partial_U V, since V is auxiliary independent;
charging that multiplier epsilon^-kappa would lose an unjustified power.
All epsilon factors not needed for the common minimum may be dropped only
because 0<epsilon<=1.

## 3. Requests at the radial shell index

Fsrc is required through m+4 before its temporal inverse. Its two radial
viscous derivatives request v through m+6. The covariance and mixed flux
in its axial derivative request only m+5. The first inverse needs Jd through
m+4. Other terms in (1) require at most m+2 fields. Thus (3) through m+6
and the angular residual through m+4 suffice; these are requests on the
FINAL reconstructed fields, not on an independent earlier stage.

At Q82's output n choose m=n+199996. Sufficient explicit ceilings are

    final fields/base/full W/Delta W through n+200002,
    derived final Etheta through n+200000,
    temporal inverse at output n+199996, using its four-jet recipe.  (6)

No 399990th derivative of the nonlinear final-state source is needed in
this split. Generating the data in (6) still has its own finite earlier
source/root/cutoff/mean requests; (6) is not the depth of the entire DAG.
The principal comparison uses Q83's separate leading table and its accepted
column/ODE value hypotheses. It is not included in Rem* a second time.

The exact remaining input inequality is that the stripped native norms in
(3), the final angular residual and the retained common-pole temporal
constant at (6) are bounded by positive polynomials in S, uniformly over
all labels on the same accepted domain. In particular Z must come from
the complete curl/signed-covariance SUBTRACTION estimates of Q51/Q57.
An absolute O(epsilon) bound for W or epsilon P alone does not establish
Z=O(1). The code evaluates manufactured rational bounds; it does not check
these actual field inequalities by assertion or by a finite derivative grid.

Combined with these hypotheses, Q82's shell proof and (5) supply the
nonlinear radial contribution at 3−5kappa. This completes the explicit
remainder NORM recipe relative to the inherited primitive field contracts.
Next join (6) to the finalized mean-update and covariance provider DAG on
its original common weighted domain, retaining all internal pressure and
cutoff terms. The complete residual ledger, infinite realization,
unconditional actual-pulse bounds and blow-up remain unproved.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`. Exact symbolic controls check (1), the
nonzero defect, the two viscosity losses and the finite request ceiling.
