# Quantitative initial scalar prefix from the complete canonical primary

Parent: `fa58d2832bd4d58e2f17fabb04ca4db63ddc781a`.
The canonical initial mean checkpoint already proved a conditional S_2
defect order qualitatively. This checkpoint supplies its explicit norm
budget from the complete real primary coefficients, via the newly quantified
fast and slow maps. It also separates the pressure's zero-Haar fast input
from its indispensable moving Haar moment. The S_2 order itself is not new.

This is the INITIAL primary/fast/slow scalar prefix, not the later
signed/sourced canonical finite-cycle seed. The reviewed primary ODE/root,
aggregate coefficient, core/background and edge hypotheses remain. No
first-band constants, full residual, infinite realization or blow-up are
certified. Historical evidence, profiles, tilt, width C_setup^-20, carriers
and the accepted common tail are unchanged.

## 1. Complete real covariance supplies a quantitative input envelope

Use the actual real canonical complete primary

    w=sqrt(epsilon)[A cos(k Phi)-k^-1 B sin(k Phi)],
    B=curl_*[n cross A/|n|^2],
    W_ab=(epsilon/2)sum_labels[A_a A_b+k^-2 B_a B_b].     (1)

All cutoffs and cylindrical curl connections stay in A,B. The mixed
quadrature averages vanish for this REAL primary only. Distinct-label
crosses vanish by the retained closed support separation. No later complex
wave is assigned that cancellation. Discrete carriers are frozen under
coefficient differentiation and k^-2<=epsilon on the accepted domain.

Let N_n be the fast checkpoint's raw coefficient norm in (U=R^d,Z,T,y).
Assume positive polynomial aggregate bounds a_n(S),b_n(S) for the label
vectors of every derivative through n of A and epsilon^kappa B, respectively.
Precisely the pointwise label-l2 derivative norms must be bounded uniformly
in the band; individual per-label bounds without an overlap/sum bound do
not suffice. Cauchy--Schwarz and Leibniz then give all six entries

    N_n(W_ab)<=epsilon Omega_n(S),
    Omega_n=2^(n-1)[a_n^2+b_n^2],                       (2)

since epsilon^(1-2kappa)<=1. The primary foundations provide the qualitative
coefficient route under their retained hypotheses; (2) records exactly
which aggregate norms are being used. It does not derive them from support
alone or evaluate the selected enormous coefficients.

## 2. The initial fast input keeps full pressure but projects a slow term

Start with zero mean velocities. The SAME full mean equation gives

    gr=-D1 Wrr-epsilon partial_Z Wrz+Wtt/R,
    P=int bar(gr), pm=T0(gr-rho P),
    Etheta=D2(Wrt-Sigma_theta)+epsilon partial_Z Wtz,
    Ez=D1(Wrz-Sigma_z)+epsilon partial_Z(Wzz+pm).          (3)

The background stress Sigma and compensator rho are independent of the
auxiliary torus. The complete same-base equation, including Pi1/Pi2, is
subtracted against itself; its pressure is not approximated by a leading
term. Base flat errors remain separate. The cutoff in T0 is slow, so
Haar averaging commutes with its translated radial integrals. Therefore

    pm^circ=T0(gr^circ),
    ftheta=(D2 Wrt+epsilon partial_Z Wtz)^circ,
    fz=(D1 Wrz+epsilon partial_Z pm+epsilon partial_Z Wzz)^circ. (4)

This projection removes rho P and Sigma ONLY from the zero-Haar fast
channel. Their full means, pressure shape and scalar moments remain.
In particular partial_Z pm^circ is differentiated as a whole moving
operator output; no derivative is commuted through its cutoff.

For an explicit graph bound, put j_n=N_n(d R^(d-1)), r_n=N_n(R^-1),

    D_e,n=2^n(3j_n+e r_n).

The exact floor M=c^rho_g epsilon^-kappa with c<=S^-1<=1 implies
M<=epsilon^-kappa. Since |vr|_1<2,

    N_n(D_e F)<=epsilon^-kappa D_e,n N_(n+1)(F).

With H=1-kappa define

    G_n=D_1,n Omega_(n+1)+Omega_(n+1)+2^n r_n Omega_n,
    Et_n=2[D_(2,n)+1]Omega_(n+1),
    Ez_n=2[D_(1,n)+1]Omega_(n+1)+2 KT_(0,n+1) G_(n+1).   (5)

Thus N_n(gr)<=epsilon^H G_n and
N_n(ftheta),N_n(fz)<=epsilon^H Et_n,epsilon^H Ez_n.
The factor two bounds the zero-Haar projection. The last term of Ez_n
uses KT_0 at order n+1 and includes every moving-cutoff derivative.
It comes from epsilon N_(n+1)(pm^circ), with an extra epsilon dropped
using epsilon<=1. Consequently Et_n uses Omega through n+1, Ez_n
through n+2. There is no independent initial fast-residual norm assumption
beyond the explicitly displayed primary and operator inputs.

## 3. Initial and postfast scalar polynomials

Let M_j=(d0-a0)max(a0^j,d0^j) on a positive hull containing the WHOLE
primary and fast potential/cutoff support. This hull need not be I4.
Put J_n=(n+1)(n+2)/2 to convert raw slow jets to normalized slow sums.
Integrating the complete covariance identities yields

    I_P,n=J_n[2M_-1 Omega_n+M_0 Omega_(n+1)],
    I_theta,n=J_n M_2 Omega_n,
    I_z,n=J_n[2M_1 Omega_n+(M_2/2)Omega_(n+1)].           (6)

The initial P,Jtheta,Jz are bounded by epsilon times these polynomials.
The factors of epsilon on the Drz derivatives were safely dropped, not
the Drz covariance itself.

Substitute Et,Ez from (5) into the fast checkpoint's actual Ft_n,Fr_n,
including g=gd-A1 gd. Denote its zero-old-mean pressure/axial fast budget
at output n by F_P,n,F_z,n (using Ft,Fr at n+1). Then

    Pi_P,n=I_P,n+F_P,n,
    Pi_z,n=I_z,n+F_z,n,
    Pi_theta,n=I_theta,n+M_2 Ft_n^2.                     (7)

Indeed all fast linear Haar terms vanish, while its quadratic changes
have order 2H=2-2kappa>=1. Old mean velocities were zero BEFORE this
fast update; this does not suppress its new velocity squares or A1 terms.
The three CURRENT postfast scalars satisfy epsilon Pi at every fixed jet.
Their old fast velocities have exactly zero Haar mean, including B.

## 4. Fresh slow solve gives explicit S_2 and the initial sharper Haar rates

Use the SAME five-map at the postfast scalars (7), with its explicit
fixed-profile inverse constants. At field order N=n+2 let

    U_N=H_N(Kp_N Pi_P,N+Kz_N Pi_z,N),
    Z_N=H_N Kg_N Pi_theta,N,
    X_N=(N+1)H_(N+1)Kpsi_(N+1)Pi_theta,N+1.             (8)

The actual slow increments obey |u|_N,|g|_N<=epsilon(U_N,Z_N) and
|B|_N<=epsilon^2 X_N. Let |b|_N<=epsilon b_N on the reserved hull,
where G=0; b is not set to zero. Because all OLD fast Haar velocities
are zero, the coupled slow polynomials specialize to

    s=2b_N X_N+X_N^2, h=b_N Z_N+X_N Z_N,
    t=U_N^2, z=Z_N^2, t_n=(n+1)+(n+1)(n+2),
    Theta_P,n=M^s_0 t_n X_N+M^s_-1(s+t)+(n+1)M^s_0 h,
    Theta_theta,n=M^s_2 U_N Z_N,
    Theta_z,n=M^s_1 z+[M^s_2 t_n X_N
                             +M^s_1(s+t)+(n+1)M^s_2 h]/2. (9)

Here M^s uses the entire reserved slow hull, separately from (6).
Directly recomputed current linear rows cancel the incoming scalars;
the remaining actual defects obey

    |Pnew,Jtheta_new,Jz_new|_n<=epsilon^2
                              (Theta_P,n,Theta_theta,n,Theta_z,n). (10)

This uses C=1, not the later-cycle C>=2-2kappa theorem. It follows
independently from the initial zero-Haar condition and 2C=C+1=2:
tangential squares have order 2; slow time, radial b-products and axial
base fluxes have order at least 3. Compact radial vector viscosity has
zero plain/weighted integrals but is retained pointwise; axial viscosity
also stays. Old fast times slow terms vanish only AFTER Haar averaging.

The initial sharper Haar majorants for future finite stages are now
a_initial,n=U_n+Z_n and r_initial,n=X_n: tangential orders 1, radial
order 2. The full fast+slow velocities retain the weaker full rates
1-kappa and 2-kappa, with bounds Ft_n+U_n+Z_n and Fr_n+X_n.
No infinite sum or complete mean-residual estimate follows from (10).

The whole axial residual moment remains

    |epsilon partial_Z(Jz_new+c_rho Pnew)|_n
      <=epsilon^3(n+1)[Theta_z,n+1+C_(n+1)Theta_P,n+1].  (11)

Equations (5)--(9) use Omega through n+10 for scalar output n and
through n+11 for (11). Specifically fast tangential/radial output r
uses Omega through r+6/r+7; postfast P,Jz through r+8 and Jtheta
through r+6; the slow u/radial B responses then use n+10/n+9.
Omega controls both A and the ALREADY CURLED B. Reducing B to its
underlying amplitude/normal data can cost another jet; n+11 is not
a claim about uncurlled amplitude data alone.

Next verify the remaining actual aggregate primary/root/edge hypotheses
and extend these initial scalar budgets to the signed/sourced prefix and
its complete pointwise residual. This checkpoint does not establish that
later finite-cycle seed, a complete residual improvement or blow-up.

## 5. Evidence

Exact Haar/moment diagnostics distinguish pressure's zero-Haar channel
from the whole moving scalar product. Independent current-source
differences retain fast squares and the zero-Haar initial simplification,
and detect its failure for a later nonzero old Haar state. Positive
polynomial evaluation checks the full norm chain and exponent ledgers.
All fixtures have finite endpoint regularity and are not selected actual
canonical fields. The conditional analytical estimates above are not
deduced from passing tests.

```sh
python3 research/navier_stokes_cascade/results/core_initial_scalar_seed_budget/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_initial_scalar_seed_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_initial_scalar_seed_budget/evidence.json \
  --output /tmp/core_initial_scalar_seed_budget_fresh.json
```
