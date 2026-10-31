# Q52: two successive viscous terms in the final angular row

Q51 raises the initial signed target to order two on its canonical
componentwise route. It does not raise the earlier fast difference. Here
we identify the *field*, including its sign, left at the limiting order
in the final angular equation. Write O_a for the existing fixed-order
coefficient bounds at epsilon^a times a positive polynomial in S, with
the appropriate shell weight and a sufficient finite pole on a common
hull. All equalities modulo O_a need the lower field estimates at the
finite higher orders actually differentiated. We use the same frozen
common chart and the original mean-map intertwiners in both updates.

Let W denote the original complete physical angular covariance component
W_rtheta in normalized velocity units (order one). It is not its Haar
average, the signed target, or the later signed covariance. Put
j=d R^(d-1), L=v_r.grad_y, I=c_i^-1 Nt^-1 on zero-Haar functions, and
Q=I_d-Haar (I_d is the identity). M and c_i are frozen; j is independent
of y. Our conditional angular estimate is

    Etheta_final = epsilon^2 M^5 j^5 c_i^-2 L^5 Nt^-2 Q W
                   + O_(3-4kappa).                                 (1)

The leading term is O_(3-5kappa). This gives one kappa of separation;
it does not assert that this upper bound is attained.

## Derivation on the unchanged first-real prefix

Before either initial mean, the projected angular source is
Q[D2(W-Sigma_theta)+Dz W_theta,z]. Sigma is auxiliary independent.
D2=partial_R+MjL+2/R, Dz=epsilon partial_Z. The initial fast angular
velocity is consequently

    u0 = -Mj I L QW + O_1,             u0 in O_(1-kappa).

The inverse has a polynomial cost and finitely many extra jets but no
epsilon cost. For either fast update at tangential order h, write its
angular velocity as u. The exact row difference after removing c_i Nt u is

    Rtheta = -epsilon M^2 j^2 L^2 u + O_(h+1-kappa).                (2)

Indeed the vector Laplacian minus its two-graph term contains at most
one graph derivative, including the derivative of j, the cylindrical
connections and the slow axial viscosity. Slow time costs h+1. Radial
transport contains the base radial velocity or the increment's radial
companion (extra epsilon); old mean products have at least the extra
old tangential order 1-kappa as well as one radial factor. Axial flux
is differentiated by Dz. This is the complete angular subtraction in
[the fast-row difference](../core_fast_row_gain_join/README.md);
there is no angular pressure gradient in this mean equation.

Use (2) first at h=1-kappa. The old slow difference is O_(2-2kappa).
The complete pointwise signed covariance response is O_(2-kappa) from
Q51 (full tensor at two and at most one graph loss). The Q49 initial
A1 time remainder is axial only. Q preserves the leading expression.
The new angular fast source and its velocity are therefore

    F1 = +epsilon M^3 j^3 I L^3 QW + O_(2-2kappa),
    u1 = -epsilon M^3 j^3 I^2 L^3 QW + O_(2-2kappa).

Now h=2-3kappa. In (2) applied to u1 the discarded source error
costs two graph derivatives and gains one epsilon, reaching 3-4kappa.
The other groups in (2) also reach at least 3-4kappa. The complete
final angular row equals the post-signed Haar row plus this difference
and the new slow difference, with the old means PLUS the new fast field
in the latter. For clarity, at target order two the three real-pair error
orders are 3-2kappa, 3 and 4-2kappa; the post-signed Haar row is at least
3-2kappa (the moment bump is order three, other complete covariance
terms enter with Dz or the compensated Haar pressure derivative).
The postfast scalars have order min(2,2,(1-kappa)+h,2h)=2, retaining
both old correlations and the single signed scalar change. Thus the
slow difference is at least 3-2kappa even charging two graph losses.
These are the same complete maps and moments as in
[the corrected initial target ledger](../core_initial_target_rate_gain/README.md),
with the Q51 target rate. Substitution gives (1). Both inverse signs and
both negative vector-viscosity signs matter. No commutation of partial_R
past j has been used in the leading operator.

The axial row also contains A1 Dz Q Phi of the **post-signed, pre-new-fast**
state. The radial pressure tail A0 Q gr uses the **final** state and still
has the separate compensator -rho P. Neither occurs as an additive term in
(1). The axial velocity gd-A1 gd, full reconstructed pressure, local
primitive terms and moving pressure moments still enter the underlying
fields and scalar estimates. Equation (1) is for the regular mean angular
row only; it does not change the physical additive-tail ledger.

## A test for improving the angular rate

No common transverse zero is needed to test the *oscillation* of W.
On the unit auxiliary torus set G=L^5 Nt^-2 QW. On a dyadic shell
K<=|k|_1<2K write l=v_r.k, t=v_t.k. The inverse symbol modulus is
|t|^2/((2 pi)^3 |l|^5). Q50's quadratic norm argument separates the
l values, with zero adjoined, by 1/(12K). Summing on either side gives
sum |l|^-10 <=4(12K)^10, while |t|<=2K. Thus the square-sum norm
of the inverse symbols is <=8*12^5 K^7. Splitting by the larger
coordinate of k, Parseval bounds the square-sum of the shell coefficients
of G by 2*(2/K)^8 J8(G), J8=max_i sup|partial_yi^8 G|. Cauchy--Schwarz
and sum 2^-h=2 prove

    ||QW||_infty <= C J8(G),     C=2*8*12^5*2*2^8.

Use (1) through eight auxiliary jets. Since j>=j0>0 and
M=c_i^rho epsilon^-kappa, recovering G costs at most
j0^-5 (Tg S)^3 epsilon^(-2+5kappa): c_i^(2-5rho)<=(Tg S)^3.
Thus an angular-row bound at 3-5kappa+eta requires the original
auxiliary oscillation of W at order 1+min(eta,kappa), with the
inherited polynomial and weight/pole allowances. The Q51 improvement
of the signed target cannot establish this inequality on its own.
If W has a zero on the common torus at the same slow point, its Haar
mean is then bounded by the same oscillation. Separate holes in
different label supports do not by themselves prove a hole in their sum.

Next compare this condition to the actual primary stress and the
simultaneous support of the original summed wave, or bound the
oscillation directly. No lower bound along physical bands, failure of
the weaker required order 3-6kappa, uniform actual-pulse enclosure,
improvement of the complete residual, or blow-up is proved.
