# Q92: explicit complete regular linear-source coefficient

Parent: `b54fe4824b9f1303b8e6d7e8653f1b7dee20f6a7` (Q91).
Q91 composes the mean-wave part using corrected full-pressure means. This
checkpoint constructs the remaining REGULAR LINEAR coefficient from the
actual full operator, conditional on its lower whole-domain wave and fixed
coefficient norms. It does not evaluate or certify those input inequalities.
The unchanged field, pressure, width, endpoint and domain contracts remain;
all historical evidence is preserved.

## 1. Exact operator and its cancellation

Use regular_source_foundation (1)--(2), on a fixed label and nonzero integer
harmonic mu, omega=k mu. Write the complete regular amplitude as a=t+r,
where t is the principal homogeneous amplitude with its prescribed regular
spatial cutoffs and r is the COMPLETE potential-curl remainder. The regular
split must preserve the principal equation of t; cutoffs constant along the
principal pulse flow commute with that equation. Any other cutoff defect
must be retained explicitly, not assumed to satisfy it. The full regular
pressure is pi. The homogeneous principal equation cancels exactly

    c Nt t+K t+epsilon omega^2 |n|^2 t+i omega n pi=0.

Consequently the regular remainder before the separately assigned temporal
commutator is exactly

    E(t,pi)+[c Nt+K+epsilon omega^2 |n|^2]r+E(r,0).       (1)

E is the FULL differential expression of regular_source_foundation (2),
including vector viscosity and cylindrical connections. Complete r is not
transverse and is not discarded by the principal cancellation. The normal
in this section is the phase normal; source derivative order is denoted j.
All actual spatial cutoff derivatives belong in t,r,pi and (1), with their
original bounds. This algebra does not supply a missing cutoff norm.

## 2. Common weighted raw-jet inputs

At output order j put D=j+2 and use factorial through-order raw coefficient
norms on the entire original enlarged lifted rectangle in the SAME R chart.
Assume, with a common pole p covering every derivative through D,

    ||t||_D <= epsilon^(1/2) sqrt(zeta) delta^-p P T,
    ||r||_D <= epsilon^(1-kappa) sqrt(zeta) delta^-p P R,
    ||pi||_D <= epsilon sqrt(zeta) delta^-p P Pi.          (2)

These mean raw derivatives divided by the displayed weight, not jets after
dividing by P and then forgetting its derivatives. T,R,Pi are positive
scale polynomials on the whole common domain. Set B=1-kappa. The signed
homogeneous addition fits these coarse orders under the retained target
hypotheses. It must use the identical original root and current target.
No independent square root or final sourced state can be substituted.

Use unweighted coefficient norms through D (pad upward when needed): b
bounds the radial background divided by epsilon; G the axial background;
g>=1 the normalized total radial-graph derivative coefficient; r0 the
reciprocal radius; N the l1 phase normal; e the phase defect divided by
epsilon; bR the radial background derivative divided by epsilon; Bz the
full axial background derivative divided by epsilon; c the pulse-frequency
coefficient; K the principal matrix. Every variable here denotes a norm
coefficient, not the physical field itself. The explicit n' and inverse-normal
costs needed to derive Pi remain upstream requirements.

At fixed j, one ordinary coefficient derivative costs at most D times a
through-D norm. The repeated radial derivative must differentiate its
multiplier. Writing Dr=sum v_i partial_i, with ||v||<=epsilon^-kappa g,
Leibniz bounds

    ||Dr^2 a||_j <= 2 D^2 epsilon^(-2kappa) g^2 ||a||_D. (3)

This requires v jets through j+1, contained in the chosen D. One term in
(3) differentiates a twice; the other differentiates v then a. Replacing
Dr^2 with a frozen-multiplier second derivative is not an exact operator.
Dz=epsilon partial_Z supplies its own gain; Nt's direction has l1 norm
less than two on the original torus. Discrete k, mu and labels are frozen.

## 3. Explicit coefficient polynomial

Let hmu=|mu|, d=j+2. The following positive polynomials bound the velocity
part of E at gain 1/2-kappa above its input order:

| Retained operator group | Coefficient |
|---|---|
| slow time, radial and axial transport | d(1+b g+G) |
| phase defect | 2 hmu e |
| background derivatives and connection | bR+Bz+b r0 |
| amplitude vector viscosity | 2 d^2 g^2+d^2+r0 d g+r0^2 |
| mixed phase/amplitude viscosity | 4 hmu N d(g+1) |
| phase divergence | 2 hmu N[d(g+1)+r0] |
| angular vector-Laplacian connection | 4 hmu N r0 |

Call their sum E0. The corresponding gains are respectively 1-kappa,
1/2,1,1-2kappa,1/2-kappa,1/2-kappa,1/2. They all dominate 1/2-kappa
for 0<kappa<1/4. We use k<=2 epsilon^-1/2, retaining mu everywhere.
All product bounds use the mixed multiindex Leibniz formula in factorial
through norms; no scalar-only differentiation convention is substituted.

The pressure gradient has coefficient d(g+1) Pi at order B. The principal
operator on r has coefficient

    H0=2 c d+K+4 hmu^2 N^2.                             (4)

It preserves order B; in particular Nt is a coefficient derivative, not an
extra radial graph loss. E(r,0) is better by 1/2-kappa, so it can be
embedded at B on 0<epsilon<=1. If Psi bounds the existing temporal cutoff
multiplication in the REGULAR grammar through j, (1)--(4) yield

    Clin_j=Psi{E0 T+d(g+1)Pi+(H0+E0)R},
    ||f_lin||_(j,weighted) <= epsilon^B Clin_j.           (5)

The weighted norm removes sqrt(zeta) delta^-p P. Fixed coefficient norms
in the table are unweighted, so products do not erase the wave factor or
add a new edge pole. Poles of t,r,pi must be padded BEFORE using (5).
An actual first-real Fourier sum uses its original finite coefficient mass
and the largest required fixed harmonic norm, never output-mode division.

This removes an independent Clin norm callback relative to explicit input
inequalities (2) and the fixed table. It is not a proof of (2) from the
existence of the norm program. The code returns each operator group as
well as its sum, so pressure, repeated graph differentiation and principal
viscosity on r can be inspected separately.

## 4. Regular versus cut wave slots and the frontier

Q27's regular and actually cut wave templates remain distinct. Clin uses
the regular T,R,Pi and Psi. Q91's nonlinear U uses the COMPLETE actually
cut velocity, including cutoff curl. One cannot infer the former from the
latter by dividing by a vanishing cutoff. Temporal commutators remain
separate additive tails under their original Gaussian separation theorem.
Spatial/transverse cutoffs, Borel jets and later localization are not put
in that tail merely because a cutoff is smooth. Their proper coefficients
and any noncommuting defects must remain in the exact regular split.

The conservative source-input request is j+2, with root/target jets through
j+3 for a complete-curl coefficient. Deriving Pi at j+2 requires moving-normal
jets through j+3 and inverse-normal/matrix jets at the corresponding order.
These are local requests before expanding the complete current target and
mean dependency DAG. No uniformity in derivative order or construction
stage is asserted; no later sourced compact source is transplanted here.

Next supply (2), actual cut U, and the fixed table on the SAME current
factory's whole lifted R domain from the original-root regular/cut templates,
then bound the wave-wave term and add (5) to Q91. Q20--Q21 and Q90 convert
that actual full-source estimate to the inverse's normalized path, entrance
and normal-pressure requirements. No unconditional pulse bounds, complete
residual improvement, infinite realization or blow-up follows here.

Exact diagnostics independently differentiate a variable radial graph,
check mixed factorial derivative estimates, and detect pressure/curl
omissions. The inherited source-foundation suite checks the full vector
operator on 672 manufactured exact linear-jet components. These tests are
not continuum certificates for (2). Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json`.
