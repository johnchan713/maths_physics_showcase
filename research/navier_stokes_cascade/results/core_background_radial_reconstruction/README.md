# Full radial-background input from the finalized coefficient and primitive

Parent: `ced9ea1e0cf26aed8c417ef957970ad029e40649`.
The parent supplies full-background and Cphi recipes but leaves the leading
velocity-tuple second jets and first-coefficient tuple first jets as inputs.
This checkpoint derives those LOW-ORDER tuple bounds from the existing
finalized F,U estimates, with the exact incoming radial primitive. It joins
the resulting constants to the parent's seven-term normal comparison.

The conclusion is an analytical implication of the retained regional,
core/Volterra, cutoff and exact moment-repair estimates. No selected giant
constant is numerically evaluated. The normalized Borel schedule, its cutoff
Euler norms, actual leading gaps, domain geometry, selector substitution,
tangent propagator and remaining supports retain their premises. Historical
frozen evidence is unchanged. No field, root, cutoff, width C^-20, band,
pressure gauge, or mean is reselected. No unconditional pulse, complete
residual improvement or blow-up follows.

## 1. Reconstruct the radial fields before taking norms

Use the original native coordinates X,eta, Rprof=sqrt(2X),
Lprof=1-2h eta^2, d=1-eta^2, D=1/2-h, A=1/2+h, lambda_n=2nh.
For each finalized coefficient (including n=0) the identities are

    M_n(X,eta)=integral_0^X U_n(Y,eta)dY,
    V_n=[2eta X U_n-2(D+lambda_n)eta M_n-d M_n,eta]/Lprof,
    E_n=Rprof F_n, W_n^r=V_n/Rprof,
    J_n^r=2eta M_n/(Lprof Rprof).                            (1)

The primitive is from the AXIS. Local vanishing of U does not reset M, V,
W or J. The full signed source integrals through the original Xv determine
the first moment repair used below; a cutoff at Xcut, I3 or I4 cannot replace
that solve. The pressure Pi1/Pi2 memories remain the inherited full ones.
This norm reconstruction does not establish a new pressure-support theorem.

On a positive annulus Rlo^2/2<=X<=Xhi, choose Rhi^2>=2Xhi. The U estimates
used for M MUST hold on all 0<=X<=Xhi, using the established zero extension
beyond Xv if needed. F bounds need only hold on the annulus. Positive-radius
E=Rprof F derivatives are not declared regular X derivatives at the axis.

For r>=1 the exact identities

    partial_X^r partial_eta^k M_n
       =partial_X^(r-1) partial_eta^k U_n,
    Lprof partial_X V_n
       =2(A-lambda_n)eta U_n+2eta X U_n,X-d U_n,eta       (2)

avoid a spurious extra radial derivative. For r=0, derivatives of M are
bounded by Xhi times the corresponding whole-interval derivative of U.
In particular a first eta jet of W contains M_etaeta. A C_eta^1 bound on U
alone does not supply it. The finalized first C_eta^3 estimate already
supplies the required C_eta^2 row; no new differentiability hypothesis is
needed on that branch of the existing proof.

## 2. An explicit positive finite-jet recipe

Fix m=0,1,2. Let f_r,u_r bound the FACTORIAL norms

    |partial_X^r F_n|_(C_eta^(m+1-r)) <= f_r,
    |partial_X^r U_n|_(C_eta^(m+1-r)) <= u_r, 0<=r<=m.     (3)

Hence a raw eta derivative of order k is bounded by k! f_r or k! u_r.
That factorial cannot be discarded. Work with nonnegative bivariate
polynomials in formal increments x,z, truncated to total degree m.
Their coefficients majorize absolute NORMALIZED Taylor coefficients.
Let

    Uhat_(r,k)=u_r/r!, Fhat_(r,k)=f_r/r!, r+k<=m,
    Mraw_(0,k)=Xhi k! u_0,
    Mraw_(r,k)=k! u_(r-1), r>=1,
    Mhat_(r,k)=Mraw_(r,k)/(r!k!),
    Metahat_(r,k)=Mraw_(r,k+1)/(r!k!).                     (4)

The extra pure-angular input in (3) supplies the last entry for r=0.
For r>=1 it uses only an already supplied lower radial row, by (2).
Set etahat=1+z, Xhat=Xhi+x, dhat=1+2z+z^2 and ell=1-2h.
An explicit reciprocal Lprof majorant is

    b_0=ell^-1,
    b_k=[4h b_(k-1)+2h b_(k-2)]/ell,
    Lhat=sum_(k=0)^m b_k z^k, b_-1=0.                    (5)

Indeed Lprof(center)>=ell, and its zero-constant variation is bounded
coefficientwise by 4h z+2h z^2. This argument is valid also at eta=+-1.
For the native radius use

    Rhat_0=Rhi,
    Rhat_j=|(1/2)_j|2^j Rlo^(1-2j)/j!, j>=1,
    Ihat_j=|(-1/2)_j|2^j Rlo^(-1-2j)/j!, j>=0,           (6)

where (a)_j is the falling factorial. These retain every inverse-radius
Jacobian term. With products truncated only after multiplication, define

    Vhat=Lhat[2etahat Xhat Uhat
              +2(D+lambda_n)etahat Mhat+dhat Metahat],
    Ehat=Rhat Fhat, What=Ihat Vhat,
    Jhat=2Ihat Lhat etahat Mhat.                           (7)

All coefficients are nonnegative rational expressions in (3) and the fixed
geometry. Multiplication bounds each raw multivariate Leibniz sum by its
normalized convolution. For each component and r+k<=m, multiply its (r,k)
coefficient by r!k! to recover a RAW bound. The maximum over these four
components and all degrees <=j defines Mtuple_(n,j), j<=m.

This is finite real jet calculus. It needs no analytic extension of the
final compact profile and introduces no physical Q power. All constants
are uniform in the physical band when their native inputs are. It does not
certify weighted stress or a pulse propagator merely from velocity bounds.

For clarity, at m=1 put c=D+lambda_n and

    H0=Xhi u_0(3+2c),
    V0=H0/ell,
    VX=[2Xhi u_1+(3+2c)u_0]/ell,
    Veta=Xhi u_0(8+4c)/ell+4h H0/ell^2.                   (8)

Sufficient radial bounds are W0=V0/Rlo,
WX=VX/Rlo+V0/Rlo^3 and Weta=Veta/Rlo. The primitive-curl bounds are

    J0=2Xhi u_0/(ell Rlo),
    JX=2u_0/(ell Rlo)+2Xhi u_0/(ell Rlo^3),
    Jeta=2Xhi u_0[2/ell+4h/ell^2]/Rlo.                    (9)

Equations (8)--(9) are also conservative at a zero-U plateau. They do not
replace the incoming primitive by its local zero velocity.

## 3. Use already available leading and first-coefficient rows

For the leading tuple at m=2, the required inputs are precisely the existing
anisotropic rows (factorial angular norms):

    F0: (F3,F_X,2,F_XX,1), U0: (KU,R1,R2).                 (10)

`global_swirl_envelope`, equation (6), supplies the first triple on its
whole required native interval. `global_axial_source_envelope`, equation
(3), supplies the second on the whole axial support, with zero extension.
These use the same Kaxis, H, Ahi, N and original regional profiles. In
particular R2= max(Kaxis,H^3 N/Xa^2) retains the mixed-frequency cost;
KU alone would not bound U_X or U_XX. The tuple second jet needs the pure
third eta derivative of U, already included in KU. No global mixed fourth
or unproved high angular row is assumed.

For the first coefficient, use the existing whole-inner C_eta^3 radial
rows J_(1,r), r=0,...,3. The first-coefficient boundary gives explicitly

    J_(1,r)=c3[r!(2Lambda)^r Mcore+p_r Mcollar+q_r],
    c3=1+2/delta+4/delta^2+8/delta^3.                       (11)

Keep the selected eta-independent coefficient cutoff k_r from
`explicit_coefficient_cutoff` and form its full product rows P_r.
`global_first_coefficient_c3` supplies the full-pressure linear repair
coefficient C13, including Q03,Q13 through the true Xv, the nonzero
axis-viscosity term, the five discrepancies, 5/e_star and lambda_patch^-1.
`global_second_coefficient_h23.finalized_first_jets` already gives

    U1,r^*=P_r+C13 B_U,r,
    F1,r^*=P_r+C13 B_EtoF,r.                               (12)

Its first bump derivatives include 121600/p and
(121600+160)/(p r_p); the F repair's inverse radius is not lost.
The two U bumps and three E bumps have the established disjoint supports.
The wrapper verifies that the whole-inner value J_(1,0) used by (12) is
the SAME local seed entering the original five-row solve. Neither a second
repair nor a new root is introduced. Set the m=1 rows in (3) to the r=0,1
entries of (12), retaining the available C_eta^3 norms as upper bounds.
Equations (4)--(7) then supply Mtuple_(1,1) explicitly from the existing
selected-parameter first-coefficient construction, rather than an independent
radial-velocity or primitive seminorm.

The code calls the historical formulas directly. It does not materialize
their selected enormous Picard exponentials; rational fixtures test the
substitution. The conditional continuum estimates of the historical
construction, not those fixtures, justify (10)--(12).

## 4. Cover the entire native image before joining Cphi

Native Rprof and normalized BAND radius R are different. Since
X=R^2/(2s), 1/4<=s<=4, Rmin<=R<=Rmax, the full chart image requires

    Rlo <= Rmin/2, Xhi >= 2Rmax^2, Rhi^2 >= 2Xhi.          (13)

All the input bounds (10)--(12) must apply on that same enlarged native
interval and on the primitive's preceding interval. The implementation
rejects a merely central-slice annulus: for example native radii [1,2]
do not cover the band-radius [1,2] chart. This makes the parent's stated
enlarged-domain premise explicit; it does not assert its hypotheses failed
or change any support. If a historical row ends before the required Xhi,
its domain must be proved sufficient or enlarged by the defining estimates,
not extrapolated from compactness or silently relabeled.

Apply (7) at n=0,m=2 for Lead_1,Lead_2 and at n=1,m=1 for M_(1,1).
With the actual Borel cutoff Euler K2>=1, the parent's increment constant
is now the explicit

    C_(1,1)=max(1,2K2(1+2h)Mtuple_(1,1)).                  (14)

Together with the SAME normalized infinite diagonal, its half-tail gives
Bcorr_1=C_(1,1)+1/2. Feeding these values into the parent's finite chart
composition supplies all seven Cphi groups without a free low-order
velocity-tuple norm. This includes the radial cutoff commutator through J,
full axial transport and whole-path unit nonzero-carrier rounding.
The Borel cutoff chi(a_n q) in (14) is NOT the coefficient cutoff kappa(X)
in (12). Its Euler K2 cannot be replaced by 729 or 2916 without a proof
about that SAME chi and derivative convention. No Borel schedule is chosen.

The unchanged normal-floor selector compatibility remains

    beta Lreq^2 >=144 r0^2 Cphi^2,
    beta=lambda_min/[4(1+U^2)^2],                          (15)

along with all original selector, rounding, domain, gap and cone premises.
Deriving a finite Cphi does not imply (15) for the original first band.
No accepted band or length parameter is increased by this checkpoint.

## 5. Checks, precise remaining inputs and scope

The exact checks compare all four reconstructed components' raw mixed jets
through order two with independent polynomial profiles and their exact
axis-to-X primitive, including eta=+-1. Plateau controls retain V=-mu when
U=0,M=mu eta and distinguish the positive-order lambda_n term.

The extra angular demand has an exact abstract obstruction. Fix smooth
unit-mass radial bumps beta_minus,beta_plus strictly before/after Xobs>0,
and put U_N=(beta_minus-beta_plus)cos(N eta)/N. If B bounds the fixed
radial factor, its C_eta^1 norm is at most B(1+1/N)<=2B. At Xobs, U_N and
all its local radial jets vanish, while M_N=cos(N eta)/N and

    partial_eta W_N(Xobs,0)
       =[N-2(D+lambda_n)/N]/sqrt(2Xobs).

This grows without bound although the final primitive is zero after the
second bump. The controls check the rational angular factors; B and the
fixed radius are harmless fixed multipliers. This is a smooth norm-input
counterfamily, not PDE coefficients or a failure of the actual selected
profile. Likewise cos(N eta)/N^2 has factorial C_eta^2 norm at most
1/2+1/N+1/N^2<1 for N>=4 but raw second derivative of magnitude one at
zero, detecting the factorial-to-raw error. None of these finite controls
samples the actual selected field or substitutes for a continuum proof.

The audit binds seven frozen evidence records and all their declared source
hashes. Reproduce with `python test_audit.py -v` and `python audit.py
--verify-record evidence.json --output /path/to/fresh.json`. Frozen outputs
cannot be overwritten. The original pressure and repair evidence is preserved.

Next verify the actual native/band coverage (13), Borel Euler K2 and the
retained normalized schedule, and insert the already derived leading gaps
and this explicit Cphi into (15) on the unchanged selector. The homogeneous
tangent propagator, actual original-root waves, complete supports and fixed
maps still require their stated checks. Then quantify the first sourced
complete covariance and fresh full-pressure/axial mean update. Uniform
physical pulse bounds, the complete viscous/cutoff/interaction residual and
blow-up remain unproved.
