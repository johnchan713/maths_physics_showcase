# Corrected diagonal field: closed-time atlas and punctured terminal jets

Parent `b74a2cce993f585e3c4d77cac26d3912ce7b5992` constructs one
conditional locally finite diagonal field from the synchronized corrected
finite states.  This checkpoint proves the corresponding closed-time
implication for that **same** field: its complete potential, separately
divergence-safe direct swirl, velocity and fixed-gauge pressure have
compatible smooth one-sided jets at `t=1` on punctured local compacts.

The result remains conditional on the retained all-fixed-order background,
root, inverse, compact-map, support, physical-transfer and closed-chart
estimates.  The exact controls check endpoint arithmetic and failure modes;
they do not evaluate the actual terminal field or prove those continuum
contracts.

## 1. Use the same complete representatives

Keep `tau=1-t`, `A=1/2+h`, `D=1/2-h` and

    q=tau+z^2 q^(2h),       X=r^2/(2q),       eta=z/q^D,
    L=1-2h eta^2 >= 1-2h > 49/50.                         (1)

The corrected diagonal has

    B_j=3/2-2 kappa+j(1/2-3 kappa),       kappa=1/100000. (2)

Write its actual seed in the operator-selected gauge as

    u^[0]=curl(A_seed)+B_seed.                              (3)

`A_seed` contains the summed background potential and every finite seed
wave and mean potential.  `B_seed` contains the base and finite direct mean
swirls.  All of these are actual representatives supplied by the retained
construction, not potentials inferred from a velocity estimate.

With the exact cutoff scales and increments from the parent, define

    A=A_seed+sum_(j>=1) chi(a_j q) A_j,
    B=B_seed+sum_(j>=1) chi(a_j q) b_j e_theta,
    p=p^[0]+sum_(j>=1) chi(a_j q)(p^[j]-p^[j-1]),
    u=curl(A)+B.                                            (4)

Equation (4) is precisely the parent's field.  Every cutoff-curl
commutator stays inside `curl(A)`.  Direct swirl remains separate and
axisymmetric.  The full seed pressure, both exceptional same-base pressure
tails, all changing covariance pressure and axial/moving-moment memory stay
in one fixed gauge.  No pressure reset, new label, new diagonal schedule or
endpoint forward inverse is introduced.

The seed differs from older current-state seeds by only finitely many
actual operations.  Under the inherited all-fixed-order endpoint contracts,
those finite additions have the same closed supports, smooth zero
extensions, one-sided slow-parameter jets and uniform closed-chart bounds.

## 2. Closed enlarged atlas and the exact endpoint condition

Use the full accepted chart range

    1/4 <= q/Q <= 4.                                       (5)

On `c<=q<=qmax<1`, the active dyadic bands lie in the finite closed window

    max{ell0,ceil(log2(1/(4qmax)))}
       <= ell <= floor(log2(4/c)).                          (6)

For `c=2^-24,qmax=2^-20`, this is `ell=18,...,26`; the narrower
`1/2<=q/Q<=2` formula would omit both accepted boundary bands.  Every
closed slow box has finitely many mesh centers, including the already
selected terminal-time indices `a_T=-1,0,1`.

Put

    q_e=min(qstar/16,2^(-ell0-2)),
    tau<=q_e/16,       |z|<=(q_e/16)^D.                     (7)

Then `q<=3q_e/8`.  Since `Q<=4q`,

    4Q<=16q<=6q_e<=3qstar/8,
    Q<=3*2^(-ell0)/8.                                     (8)

Thus every enlarged chart stays strictly inside the inherited physical
domain and above the admitted first-band threshold.  This restriction is
only an observation cylinder; it does not retune the field or labels.

The diagonal scales satisfy `a_j>2^j/qstar`.  If `q>=c` and stage `j` is
active, then `a_j c<1`, hence

    j < log2(qstar/c).                                     (9)

Only finitely many stages and bands survive on each such compact.  The
arbitrarily growing stage constants are therefore harmless at the endpoint.
The decisive inherited inequality is

    sup_(E+ x [1-delta_t,1))
      |D_x^alpha partial_t^s(A,B,p)| < infinity             (10)

for every fixed `alpha,s`, including the next time derivative.  The
closed-pole physical pulse bounds, one-sided compact maps and finite source
grammar supply (10) conditionally.  If any of those estimates fails to be
uniform through `|eta|=1` or `T=0`, equation (10) is the precise blocker.
Interior smoothness, bounded velocity or a flat residual at the origin do
not imply it.

## 3. Cover every punctured terminal compact

Let a terminal compact in (7) satisfy `|x|>=d>0`.  Choose

    b=min(q_e/4,d^2/(16Xext),1/4),
    zeta=min(d/4,(b/16)^D),       tau<=b/16.                (11)

On `|z|>=zeta/2`, equation (1) gives a positive lower bound for `q`, so
(10) applies.  On `|z|<zeta`, one has `r>d/2`, `q<=3b/8`, and

    X>=d^2/(3b)>=16Xext/3>Xext.                             (12)

This entire equatorial neighborhood is the unchanged heat exterior.  The
two regions overlap and cover positive-radius equatorial points and axis
points with `z!=0`.  The singular origin is deliberately excluded.

The exterior pressure conclusion uses the retained global pressure moments
and support hull, not velocity support alone.

## 4. Full exterior heat and pressure jets

In the heat exterior the fixed-gauge representatives are `A=0`, direct
swirl `K e_theta`, and

    K=C r^(-1-2h) H(4tau/r^2),
    H(Z)=Gamma(1+h)^(-1) integral_0^infinity
         exp(-v)v^h(1+Zv)^(-h) dv,
    p=-integral_r^infinity K(rho,tau)^2 d rho/rho.          (13)

For `A_j=(h)_j(1+h)_j` and
`S_j=sum_k binomial(j,k)A_k A_(j-k)`, dominated differentiation gives

    partial_t^j K(r,1-)=C 4^j A_j r^(-1-2h-2j),
    partial_t^j p(r,1-)=-C^2 4^j S_j
       r^(-2-4h-2j)/(2+4h+2j).                             (14)

The complete Leibniz convolution `S_j` is essential; already at `j=2`,

    S_2-2A_2=2A_1^2>0.                                    (15)

The exact radial recursion

    c_(n+1,k)=-(a+n+2k)c_(n,k)+2c_(n,k-1)                 (16)

retains both the changing homogeneity and the derivative of
`Z=4tau/r^2`.  It yields all fixed mixed radial/time bounds on
`r>=rmin>0`.  Cartesian frame derivatives add only fixed inverse-radius
costs and axial derivatives vanish.

For nonzero exterior constant `C`, the time Taylor coefficient ratio is

    (4/r^2)(h+j)(1+h+j)/(j+1)>4j/r^2.                      (17)

Moreover `A_j>=h(j-1)!j!` and `S_j>=2A_j`.  Hence the exterior swirl and
full pressure have zero-radius terminal time Taylor series.  They still
have smooth one-sided jets.  A later extension must use smooth jet gluing;
analytic continuation or zero extension of nonzero pressure jets is not
licensed.

## 5. One-sided traces and remaining boundary

Using (10) at the next time order gives

    ||D_x^alpha partial_t^s G(t2)-D_x^alpha partial_t^s G(t1)||
       <=M_(alpha,s+1)|t2-t1|,       G=(A,B,p).             (18)

Every fixed jet therefore has a uniform one-sided terminal limit on an
enlarged punctured compact.  Spatial coordinate integrals identify all
spatial derivatives, and time integration identifies successive time
jets.  For every fixed finite `k`,

    ||D_x^alpha[G(1-tau)-sum_(j=0)^k
       (-tau)^j T_j^G/j!]||
       <=M_(alpha,k+1) tau^(k+1)/(k+1)!.                   (19)

Curl commutes with these finite limits.  A velocity `m`-jet needs
potential order `m+2` and direct-swirl order `m+1`; pressure needs `m+1`.
There is no uniform-in-order bound and no interchange with an uncut
infinite series.

This closes a **conditional punctured endpoint implication for the same
corrected diagonal field**.  It supplies no regular velocity or pressure
at the singular origin and no future solution.  Next join axisymmetric
spatial/time localization and zero initialization to this field, retaining
all curls, pressure, viscosity and nonlinear terms; then construct the
smooth compact force, prove whole-space energy/dissipation, and audit the
global comparison.  The primitive analytical contracts remain before any
unconditional pulse, forcing, complete-residual or blow-up claim.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
