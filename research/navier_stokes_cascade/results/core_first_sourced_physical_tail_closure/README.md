# Conditional first-sourced compact physical-tail closure

Parent `e028fb3687f788bdf92c456d6dcb1becb61c3f7e` separates the compact
tail closed for the first-real state from the improved regular residual of
the changed first-sourced state.  This checkpoint carries the dependency
graph through that first sourced wave, its complete complex covariance,
full pressure and moment changes, and the actual fresh fast/slow means.

For every fixed physical jet `m` and requested power `q^N`, the changed-state
`A0/A1` inputs are finite acyclic outputs of the same named all-fixed-order
primitive contracts.  Thus the existing raw compact estimates conditionally
make the **post-sourced** compact vector flat on the unchanged accepted
domain.  This closes the stage-synchronization gap as a conditional
implication; it does not certify the primitive continuum contracts or improve
the regular residual.

Historical evidence, fields, roots, widths, cutoffs, gauges and accepted
bands are unchanged.  The exact controls distinguish finite dependency
bookkeeping from numerical evaluation of the enormous requested jets.

## 1. The changed-state graph is acyclic

Let `R` be the first-real state and let `v` be the first sourced complete
cut curl with its full forced pressure.  The construction order is

1. derive the regular nonzero source of `R` at `B=1-kappa`;
2. solve the constrained sourced equation with zero entrance data and cut
   the potential and pressure before taking the complete curl;
3. form every entry of

       Delta W=Bil(w_R,v)+Cov(v);                           (1)

4. reconstruct `Delta P`, `Delta p_m`, `Delta Jtheta` and `Delta Jz`,
   retaining the moving `c_rho Delta P` product and all old pressure memory;
5. form the postwave `Qcirc` tangential rows;
6. apply the exact fresh fast map, including the actual axial velocity
   `g_d-A1 g_d` and its `A1` response from that postwave row;
7. recompute pressure and all three scalar defects after the fast change;
8. apply the reserved five-map slow response to those current defects;
9. form the final full radial source `gr_S`, its unit-mass compensated
   pressure and only then `A0 Qcirc gr_S`.

The first sourced wave depends only on `R`.  The postwave `A1` response is
an input to the fresh fast field, but that response depends on the pre-fast
postwave row.  The final `A0` output is not used by any earlier node.
Therefore the graph contains no compact-output cycle.

At the final first-sourced state the exact compact residual is

    Tcompact_S=(-A0 Qcirc gr_S,0,+A1 Qcirc Ez_postwave).    (2)

The axial source here is the same pre-fast postwave row appearing in the
fresh-fast identity.  It is not the stale first-real `Phi_R`, and the radial
source is not `gr_R`.  Equations (1)--(2) keep the complete changed
covariance, axial velocity and local mean pressure.

## 2. Finite backward request ledger

Put

    B=1-kappa,   D=3/2-2kappa,
    Lambda=3/2+2h,   0<h<1/100.                            (3)

The changed final radial source has order `B`; the postwave axial row has
order `D`.  For physical jet `m` and requested power `q^N`, the retained raw
compact rule chooses

    p0=max(1,floor(((N+Lambda+2m)/h-B)/kappa)+1),
    p1=max(1,floor(((N+Lambda+2m)/h-D)/kappa)+1),          (4)
    J0=m+2p0+3,   J1=m+2p1+3.                             (5)

These are sufficient strict-rule counts, not necessary lower bounds for
every possible proof.  The complete post-sourced call graph gives the
following safe backward recurrence.

- `A1 Qcirc Ez_postwave` reads its postwave source through `J1`.  Its nested
  original compact recipe reads the original first-real means through
  `J1+200004`.
- `A0 Qcirc gr_S` reads `gr_S` through `J0`, hence final fields through
  `J0+3`.
- A final field requested through order `r` reaches the fresh slow row at
  `r+1`, its postfast state at `r+4`, the fresh fast radial row at `r+9`,
  the nested old `A1` raw source at `r+200010`, and its original first-real
  mean input at `r+200013`.
- The sourced-wave/covariance branch from the same final-field request reaches
  the original first-real state no later than `r+28`; it is dominated by the
  retained nested `A1` branch.

Thus one sufficient common first-real input ceiling is

    Jstar=max(J1+200004, J0+3+200013, J0+3+28)
         =max(J1+200004, J0+200016).                       (6)

For the normalized final-row target `C=3/2-kappa`, the existing counts are
`J0=100005`, `J1=7`, and (6) gives `Jstar=300021`.  For the illustrative
physical values `h=1/200`, `m=N=0`, `kappa=1/100000`,

| branch | integrations | raw source | first-real ceiling |
|---|---:|---:|---:|
| radial `A0` | 30,100,002 | 60,200,007 | 60,400,023 |
| axial `A1` | 30,050,003 | 60,100,009 | 60,300,013 |

The radial branch dominates.  All values are finite symbolic requests.
Nothing here evaluates a derivative, factorial, profile seminorm or selected
giant pulse at those orders.

## 3. Conditional high-jet implication and physical flatness

Each edge in the graph is one of the already retained finite operations:
Leibniz product, chart composition, root/reciprocal recurrence, constrained
ODE variation, complete curl, harmonic convolution, Haar integration,
compact primitive, pressure/moment map, or fixed fast/slow inverse.  For any
one finite order these maps preserve finite positive-polynomial constants,
provided their named normal/root/value gaps, same-domain weights, full-path
supports and smooth endpoint extensions hold.

Consequently the primitive all-fixed-order contracts imply, at the orders in
(5)--(6),

    N_J0(Qcirc gr_S) <= epsilon^B P_J0(S+K_all),
    N_J1(Qcirc Ez_postwave) <= epsilon^D Q_J1(S+K_all).    (7)

The original separate `A0` and `A1` kernel profiles, width, gap, edge floors,
common pole embeddings and whole-line extension are retained.  Applying
(4) and the same-domain polynomial absorption yields

    |Dphysical^m Tcompact_S| <= C_(m,N) q^N               (8)

for each fixed `m,N`.  Constants may depend on that fixed request and the
fixed profiles.  Equation (8) is not uniform in derivative order, stage,
width, radius, `h` or an infinite diagonal.

The original Gaussian commutators and assigned base/exterior remainder stay
as separate conditional flat summands.  They are not sent through the
sourced inverse.  Spatial/slow cutoffs remain inside the complete regular
operator with viscosity, pressure and every interaction.

## 4. The complete residual is still blocked by the regular order

The synchronized first-sourced split is now conditionally

    Rphysical_S=Rregular_S+Tcompact_S
                    +TGaussian_S+Tbase_exterior.           (9)

Its regular nonzero order is still `D`.  Before polynomial absorption its
physical exponent is

    hD-Lambda-2m=-3/2-h(1/2+2kappa)-2m<0.                (10)

After reserving `eta>0`, a later complete regular order `a` can imply
`q^N` only if

    h a > N+Lambda+2m+eta.                                (11)

The present `D` fails (11), even for `m=N=0`.  This is an insufficiency of
the available upper bound, not a residual lower bound.

The next obligation is the later common-family sourced induction with a
stage-synchronized compact graph at every retained stage, followed by a
controlled physical diagonal and the remaining localization/forcing product
estimates.  Independent validation of the primitive background, cone,
normal, pulse, support and base/exterior contracts also remains open.  No
unconditional pulse family, complete residual decay, infinite realization or
blow-up is asserted.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
