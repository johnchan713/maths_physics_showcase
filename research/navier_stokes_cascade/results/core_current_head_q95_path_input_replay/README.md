# Q95: explicit whole-path, Gaussian quotient and normal-forcing input replay

Parent: `26c44bc5c8c48948e711742e6cbcc1e0572467aa` (Q94).
Q94 supplies the three-group first-real regular source coefficient F_j at
B=1-kappa and pole p_j=2p_(j+2). This checkpoint composes THAT stateful
source with Q20/Q21 and Q90, keeping finite derivative quantifiers, poles,
raw-versus-factorial norms and unprojected pressure forcing explicit.
The result remains conditional on the original whole-domain native provider
inequalities. No frozen evidence, field, cutoff, pressure memory, width,
selector or accepted q-domain is changed. No sourced response is built yet.

## 1. Fix the requested ceiling before embedding the source jets

Fix R first and compute Q94 F_j for EVERY 0<=j<=R on the same exact state,
label and whole enlarged lifted rectangle. Q94's factorial through-j norm
gives the raw derivative bound j! F_j, not F_j alone. Choose

    pF=max_(0<=j<=R) 2p_(j+2).                         (1)

Because delta<=1, every source jet embeds into this one common pole before
any convolution. No monotonicity of independently supplied coefficient
polynomials is assumed. A fixed-jet estimate cannot be treated as an
all-order norm or as an arbitrary later-stage source bound.

Q21's original affine lift holds slow R,Z,T fixed, keeps the transverse
coordinate, and satisfies v(H_w)=w. The column absolute sums are <=2 for
auxiliary columns, <=2/S_* for the w column, and one for slow columns.
All higher map derivatives vanish. Thus a mixed raw path jet with l slow,
k auxiliary and r pulse-time derivatives, j=l+k+r, obeys

    |Dslow^l Daux^k Dw^r f_R(x,H_w)|
      <= 2^(k+r) S_*^-r epsilon^B sqrt(zeta) delta^-pF
                                               P_i(w,x) j! F_j(S). (2)

The accepted S_*>=1 gives the uniform raw majorant

    A_j=2^j j! F_j.                                    (3)

The original Gaussian P_i is evaluated at the SAME mapped point; there
is no ratio between a Gaussian at w and one at a physical trace. The
original lift, labels, length and endpoints are frozen when differentiating.
Neither a new smaller-q floor nor iteration of Q26's chart factor 16 occurs.
Physical graph derivatives and moving physical corner charts are separate.

## 2. Quotient jets with explicit log-pole bookkeeping

Suppose ell_j bounds every raw order-j mixed jet of log P_i on the WHOLE
path, with finite pole e_j>=0. Q20's unweighted log-Gaussian contract is the
special case e_j=0. If such unweighted bounds have not been supplied, poles
cannot silently be discarded. Constants and polynomials ell_j are actual
lower-provider inequalities, not measurements or consequences of a norm tag.
The original fixed-label reference-rate formula supplies them when its
whole-path rate jets are proved; length/endpoint derivatives must be added
if a different chart moves these frozen data.

For either exp(+log P) or exp(-log P), raw Bell majorants satisfy

    H0=1, h0=0,
    Hj=sum_(i=1)^j binom(j-1,i-1) ell_i H_(j-i),
    hj=max_(i=1)^j [e_i+h_(j-i)].                       (4)

Here hj is a sufficient pole for Hj. Different summands may be padded to
that pole since delta<=1. The mixed multiindex Leibniz coefficients at each
total split sum to the displayed binomial coefficient. The zeroth log value
is UNUSED; neither |log P| nor 1/min P enters (4).

Writing g=f_R/P_i, apply the product rule AFTER the raw path bounds (2):

    Gj=sum_(i=0)^j binom(j,i) Hi A_(j-i),
    |D^I g| <= epsilon^B sqrt(zeta) delta^-q G_j,
    q=pF+max_(j<=R) hj, |I|=j<=R.                      (5)

This is a quantitative same-domain input to the retained inverse, with no
additional SOURCE derivative order. If all log jets are unweighted, q=pF;
otherwise the explicit Bell pole is added. Omitting derivatives of P in
this calculation would replace Gj incorrectly by Aj.

The native R chart has five coefficient coordinates. The redundant path
chart has three slow coordinates, two auxiliary coordinates and w: SIX.
The raw bound (5) is dimension independent. If a factorial through-order
norm in that six-variable chart is wanted, a sufficient coefficient is

    Gthrough_R=sum_(j<=R) (6^j/j!) G_j.                 (6)

It is not the native five-variable count. A separately reduced transverse
chart may have five independent coordinates, but its jet norm requires its
own coordinate conversion before replacing (6). This distinction changes
only finite constants; it does not produce an epsilon reserve.

## 3. Unprojected normal forcing and the pressure input

Keep b=n/|n|^2 from the SAME Q93 full native phase and inverse-normal floor.
Its native raw jets are bounded by the Q93 binomial convolution of the
normal and reciprocal-square jets. Convert their raw maxima to a native
factorial sum with 5^j/j!, then apply the same atlas and R-chart map used
by Q93. A further j! and 2^j safely yield raw path-jet coefficients B_j.
The normal/inverse-normal contracts here are unweighted, as in Q93; if a
different provider has poles, they must also be added to (5).

The exact full forced pressure remains

    pi/P_i=i[(K^T n-n').(t/P_i)+n.g]/(k m |n|^2).

Equivalently its SOURCE-NORMAL summand is

    pi_f/P_i=i b.g/(k m),
    Cj=sum_(i=0)^j binom(j,i) Bi G_(j-i).                (7)

The program computes Cj from the actual same-state normal recipe, not a
new pressure callback. For each fixed nonzero input harmonic m, 1/k<=sqrt
(epsilon), so this summand has sufficient order B+1/2=3/2-kappa and the
common pole q. This is the forcing contribution, not a bound for the entire
response pressure: the moving-normal response term needs the inverse's
own t/P_i estimates and corresponding phase jets. Normal rows need phase
jets through R; the Q94 wave/pressure branch already requests through R+3.

A bound for the tangent projection of f is insufficient: a purely normal
source can have t=0 but pi_f!=0. No projection is applied to g before (7).
Pi1/Pi2 background memory, cumulative compact pressure, axial gd-A1 gd and
the differentiated c_rho P moment remain in the identical upstream means
and target generating Q94's F_j. They are not reset by Gaussian division.
Original temporal/base/exterior and compact mean tails stay outside g.

## 4. Entrance and support are separate from magnitude

A finite norm (5) does not prove any zero trace. Use the original first-real
cutoff grammar: every regular source term has a cut wave/pressure factor,
or is the prescribed temporally cut regular linear remainder, with all
required zero entrance jets on the whole lifted collar. Q90 then supplies
entrance flatness of the full unprojected f and g, and of BOTH the zero-entry
velocity response and its pressure. This is conditional on those actual
cutoff/support and smooth-multiplier premises, not a boolean certificate
in the program. Pulse entrance w=0 is distinct from global one-sided
physical time tau=0. No physical corner extension follows here.

The raw outgoing solution need not vanish. The existing construction cuts
complete potential AND pressure and retains their cutoff residual. Whole
slow/transverse path support, not a single entrance value, is what uniqueness
uses to preserve support. Old additive residual tails are not fed into the
inverse, and first-real compact tails cannot be paired with a changed state.

## 5. Result, scope and next obligation

PathInputs accepts only Q94 CurrentSource, computes all source coefficients
through fixed R, performs (1)--(6) and derives (7). It returns the original
source pole, extra Bell pole and common inverse-input pole separately. The
Q94 source-factory reserve R+2, target reserve R+16 and profile reserve R+29
remain conservative local requests before expanding the complete primitive
DAG. Log-reference-rate jets through R are an additional original provider
requirement, not new derivatives of the source itself.

Established: an explicit CONDITIONAL replay of Q94's first-real source into
normalized whole-path jets and source-normal pressure forcing. Native pulse,
root, phase, atlas, support, log-Gaussian and fixed-background inequalities
are not certified by these polynomials. Next apply the actual same-domain
zero-entry constrained inverse with these source coefficients and its retained
propagator/gap bounds, derive complete potential/curl/pressure response jets,
and recompute the changed covariance and all mean/compact sources. No
unconditional pulse bounds, complete residual improvement, infinite
realization or blow-up is established.

Exact manufactured checks compare mixed Gaussian quotient derivatives,
exercise nonzero log-pole accumulation and the six-variable count, and
verify the typed current-source composition and normal-pressure retention.
They are not selected-field continuum certificates. Run `python test_audit.py
-v` and `python audit.py --verify-record evidence.json --output FRESH.json`.
