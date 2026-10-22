# Current signed target: quantitative two-edge primitive and scalar routing

Parent: `4b32690c35288765c6915ef1229a327740b8e48c`.
The parent supplies explicit root/curl budgets for a prescribed weighted
signed target. Here we connect that input to the CURRENT projected mean
rows, with an explicit two-edge primitive estimate and the already bounded
scalar channels. This quantifies the historical qualitative edge-preserving
primitive lemma; it does not present that lemma or the canonical seed order
as a new result. All historical frozen evidence is unchanged.

The estimate remains conditional on weighted CURRENT Haar-row bounds and
the retained exact core/background, canonical pulse and mean-map premises.
The recent unweighted scalar/aggregate budgets do not alone prove those
weighted row bounds. Section 4 gives a precise obstruction to that shortcut.
No profile, width C_setup^-20, radius, carrier, pressure gauge, operation
order or common accepted band domain is changed.

## 1. Identify the actual signed map and its two scalar channels

Use e=2 for the angular row and e=1 for the axial row. After the initial
canonical primary/fast/slow prefix, let E_e be its CURRENT Haar tangential
mean row, with all old mean velocities and all six covariance entries.
The exact identities are

    M_2=int R^2 E_2 dR=epsilon partial_Z Jtheta,
    M_1=int R E_1 dR=epsilon partial_Z(Jz+c_rho P),
    int R^e sigma_e dR=1,
    h_e=E_e-sigma_e M_e,  Xi_e=-T_e h_e.                    (1)

The WHOLE moving product c_rho P remains under the derivative. The two
scalar moments are different; the subscript e does not denote their
epsilon order. The projection removes them from Xi, not from the pressure
or the mean system. Pi1/Pi2 cancel only against the same complete base
balance. A0/A1 and other additive flat residuals are not silently included
in this regular weighted row or asserted to have shell support.
Their actual cutoff/A0/A1 velocity corrections still remain in every
nonlinear product defining the current rows. Any additive residual split
must use exact subtraction; no flat velocity is set to zero.

On an auxiliary-independent source, the physical shifted T_e has Haar
formula

    R^e T_e h=(1-chi_m) int_left^R r^e h dr
                         -chi_m int_R^right r^e h dr.

The EXACT zero moment of h makes this the ordinary left primitive,
independent of chi_m. Its A_e error is identically zero. This simplification
is valid here because h is already Haar averaged; it is not the formula
for a general oscillatory source before averaging.

Put x=log(X/Xa), L=log(Xb/Xa)>0, R=sqrt(2sXa) exp(x/2),
s=q/Q in [1/4,4], and x in (0,L). Define

    a=t1^2, b=4, zeta(x)=exp(-a/x^2-b/(L-x)^2),
    delta=min(1,x,L-x),
    g_e(x,Z,T)=(R^(e+1)/2) h_e(R,Z,T).

Then int_0^L g_e dx=0 and
Xi_e=-R^-e F_e, F_e(x)=int_0^x g_e du=-int_x^L g_e du.
Both shell endpoints are fixed in this logarithmic chart even when s,
sigma_e and c_rho move. Coordinate changes back to native (U,Z,T) have
finite fixed-jet coefficients on the same positive shell; no new Q loss,
endpoint motion term or derivative-dependent first band is introduced.

## 2. An explicit primitive constant for the SAME full edge weight

Suppose |g|<=A zeta delta^-p, p>=0. Let r=min(1,L/4) and
q=max(0,ceil((p-3)/2)). Define

    D(c,p,r)=(1/(2c)) sum_{j=0}^q binom(q,j) j! (r^2/c)^j,
    C_p=max(exp(b/(L-r)^2) D(a,p,r),
            exp(a/(L-r)^2) D(b,p,r),
            L Eedge(p) exp((a+b)/r^2)).                    (2)

Eedge(p) is the parent's explicit constant for sqrt(zeta)delta^-p;
zeta<=sqrt(zeta), so it also bounds zeta delta^-p. The constants in
(2) are fixed and finite but can be enormous. They are not materialized
for the selected giant parameters or declared uniform in t1 or L.

For 0<x<=r, the opposite exponential is bounded relative to its value at x
by exp(b/(L-r)^2). The exact substitution
u=x/sqrt(1+x^2 v/a) gives

    int_0^x exp(-a/u^2) u^-p du
     =exp(-a/x^2) x^(3-p)/(2a)
        int_0^infinity exp(-v)(1+x^2 v/a)^((p-3)/2) dv.

The last integral is at most the finite sum in (2), since x<=r and
integral exp(-v)v^j dv=j!. This bounds |F| by
A zeta delta^-p C_p (and actually gains x^3 on that collar).
The right collar uses the BACKWARD primitive and c=b. On the middle
[r,L-r], zeta>=exp(-(a+b)/r^2) and |F|<=A L Eedge(p), giving the
third term of (2). Consequently, on the entire SAME shell,

    |F(x)|<=A C_p zeta(x) delta(x)^-p.                     (3)

Zero total moment is essential for the right collar. Neither small total
moment nor small unweighted norm can replace it at a flat endpoint.

Use factorial-normalized total jets J_n in (x,Z,T), and nondecreasing
p_n with J_n(g)<=epsilon^alpha zeta delta^-p_n H_n(S).
For purely slow derivatives integrate as in (3). Every jet with i>=1
x-derivatives is a jet of g of order n-1, with normalized factor 1/i<=1.
With H_-1=0 this gives

    J_n(F)<=epsilon^alpha zeta delta^-p_n
                    [C_(p_n) H_n+H_(n-1)].               (4)

All slow derivatives of the exact total moment vanish as well, so the
backward formula applies to them. Multiplication by R^-e adds a finite
positive-polynomial coefficient convolution. The native chart conversion
is treated by the composition coefficients below, not by pretending that
it is linear. Neither operation needs a higher incoming source jet.
Smooth zero extension follows because
every finite inverse-distance power is dominated by the SAME zeta.

## 3. Substitute the actual scalar bounds, then feed the original root map

Let Ehat_e,n(S) be a sufficient weighted jet bound for the pullback
(R^(e+1)/2) E_e, with edge exponent p_n. The projector bumps lie on the
reserved closed interior hull in X; hence the pullbacks
s_e=(R^(e+1)/2)sigma_e have finite weighted jet bounds B_e,n(S) with
no endpoint pole. They remain moving functions, not frozen numerical bumps.
If m_e,n bounds epsilon^-alpha J_n(M_e), products in (1) give

    H_e,n=Ehat_e,n+sum_{j=0}^n B_e,j m_e,n-j.               (5)

For this INITIAL prefix use alpha=2-3kappa. The root/aggregate/initial-scalar
chain supplies, at every fixed slow order n,

    m_2,n <= (n+1) Theta_theta,n+1,
    m_1,n <= (n+1)[Theta_z,n+1+C_rho,n+1 Theta_P,n+1],      (6)

because the actual moments have epsilon^3 and 3>alpha. The extra positive
epsilon power is safely dropped, not assigned to J itself. C_rho,n is the
normalized through-order slow bound for c_rho. Equations (5)--(6) use the
CURRENT postslow defects. They do not reset the old fast/slow velocities,
delete old/new correlations or substitute a projected Xi estimate for a
discarded scalar moment.

Let K_e,n be the output coefficient polynomial from (4) and multiplication
by R^-e in the logarithmic chart. Dividing Xi by its full zeta to form the
normalized signed input xi_e costs the differentiated reciprocal weight.
For finite
reciprocal-weight coefficients L_j, with
J_j(zeta^-1)<=zeta^-1 L_j delta^-3j, set

    Vx_e,n=sum_{j=0}^n K_e,j L_n-j,
    r_n=max_{0<=j<=n}(p_j+3(n-j)),
    d_N=max(0,max_{0<=n<=N}(r_n-3n)).                      (7)

For the smooth native-to-logarithmic chart map Psi, let C_j(S), j>=1,
bound the sum over its components of the order-j normalized derivative
norms. With C(z)=sum_{j>=1} C_j z^j define

    V_e,n=sum_{k=0}^n Vx_e,k [z^n] C(z)^k,

where [z^0]C^0=1. Finite multivariate Taylor composition proves this
bound; its higher chart derivatives cannot be omitted. Since r_k<=r_n
for k<=n, no larger edge exponent is needed. Thus, for n<=N,
J_n(xi_e)<=V_e,n delta^(-3n-d_N), exactly the input required by the
parent's ORIGINAL-root signed budget. The pole allowance d_N may depend
on e as well; use V_1,n+V_2,n and the larger pole allowance for its
two-component target. This finite-order allowance does not change the
common band domain.
The primitive itself adds no source derivative. A signed curl at order n
reads (6) through n+1: this scalar-moment channel alone can require the
underlying initial primary data through n+13. The separate weighted
Haar-row input Ehat_e,n+1 may require further jets. No total derivative
count for the complete current residual is asserted.

The remaining concrete inequality is

    sup_on_chart [epsilon^-alpha zeta^-1 delta^p_n
        J_n((R^(e+1)/2) E_e)] <= Ehat_e,n(S)               (8)

for the actual CURRENT full regular mean rows on every accepted band,
plus the retained exact moment and support identities. Historical
qualitative class claims give a conditional route to (8); the recent
unweighted coefficient/scalar constants alone do not supply it. This
checkpoint derives (2)--(7), not an independent proof of (8).

## 4. A precise obstruction to using only unweighted bounds

Choose a nonnegative smooth phi supported in (0,1/4), with integral m>0.
For eta>0 small, put

    g_eta(x)=exp(-a/(4eta^2))
        [phi((x-eta)/eta)-phi((x-7eta/4)/eta)].

Each source has EXACT zero total moment, smooth compact support in (0,L),
and every fixed unweighted derivative norm tends to zero as eta->0.
At x=3eta/2, however,

    F_eta=eta*m*exp(-a/(4eta^2)),
    F_eta/zeta(3eta/2)
      >=eta*m*exp(7a/(36eta^2)).                           (9)

Hence this quotient defeats EVERY fixed inverse-distance power bound
with an eta-independent constant. For a desired power d choose any integer
k>(d+1)/2; exp(y)>=y^k/k! makes the divergence elementary.
The smooth bounded R^-e factor does not alter it. Such sources can be
converted to radial sources by h=2R^(-e-1)g, preserving their exact
weighted moment and compact support.

This is a counterexample to an operator inference, not an actual canonical
field or a disproof of its conditional weighted estimates. It shows why
small unweighted current norms and exact moment cancellation, even with
all fixed unweighted derivatives controlled, cannot discharge (8).

Next prove (8) from the complete CURRENT mean equations and weighted
fast/slow/pressure operators, then verify the actual sourced nonzero rows
and whole-path supports. Viscosity, moving pressure, axial means, A0/A1,
cutoffs and every nonlinear interaction remain in the complete residual
audit. Unconditional pulse bounds, infinite realization and blow-up remain
unproved.

## 5. Evidence

Exact finite Gaussian-tail polynomials independently verify the transformed
primitive and its ODE; rational envelopes check (2). Current-field moment
and projection controls replay all six covariance and whole moving-pressure
identities. Polynomial projection fixtures detect a missing total moment;
scaled two-bump controls check the exponent in (9), not a canonical PDE
counterexample. Finite-order substitutions test (5)--(7) and the scalar
channel jet count. These diagnostics are separate from the analytical
inequalities above and do not certify continuum suprema.

```sh
python3 research/navier_stokes_cascade/results/core_signed_target_edge_budget/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_signed_target_edge_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_signed_target_edge_budget/evidence.json \
  --output /tmp/core_signed_target_edge_budget_fresh.json
```
