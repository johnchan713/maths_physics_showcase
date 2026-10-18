# Regular-source envelopes and weighted edges; correction to the tail selector

Parent: `ffd0b394d7f0e5df08c9d171937ae3fb1d362d8e`.
Status: `regular-source-envelope-edge-estimates-with-explicit-freezing-correction`.
On 18 October 2026 the parent's checkpoint CI step passed; the full workflow
was still running without a failed step. Historical evidence is unchanged.

This audits the functional inputs to the sourced inverse: derivatives of
f/P, complete nonlinear wave factors, signed reciprocal roots and oriented
mean primitives at both flat shell edges. It gives explicit all-fixed-order
estimates, not an evaluation of the actual giant fields. Full moving/source
normal pressure, axial velocity, all covariance entries and Pi1/Pi2 memory
remain. The result is relative to the reviewed background/frame and current
target-jet/moment estimates; no unconditional pulse/residual/blow-up theorem
or uniform-in-stage bound is claimed.

## 1. A missing displayed acceptance condition must be restored first

The parent defined E using pulse length and whole-pulse smallness, with
ell_old explicitly NOT presumed analytically accepted. That does not by
itself supply target-direction FREEZING. The earlier primary-covariance
foundation, Section 3, independently requires

    tau_freeze <= Cfreeze*S^(-3), S=ell^2,
    ell^6 >= 8*Cfreeze/gamma0, gamma=gamma0/2.                    (1)

Indeed an unfrozen cone |W0|<=(1-gamma0)V0, V0>0, and coordinate changes
|Delta V|,|Delta W|<=tau_freeze*V0 give

    V>=(1-tau_freeze)V0,
    |W|/V <=(1-gamma0+tau_freeze)/(1-tau_freeze).

The sharp half-margin threshold is tau_freeze<=gamma0/(4-gamma0).
The padded gamma0/8 suffices. The constant is obtained from the SAME smooth
limiting stress DIRECTION and closed-box derivative bounds; one must not
divide by a target that vanishes at a shell edge. An unbounded T0/zeta
does not prevent this relative-direction estimate.

No displayed domination of Cfreeze by the parent's Lreq was proved. This
is an omitted SUFFICIENT condition in that selector, not a demonstration
that the actual construction fails: its actual fixed constants have not
been numerically evaluated. The parent theorem remains conditional on its
retained cone inputs; its claim to remove all value acceptance needs (1).

The corrected INITIAL choice in that construction template is

    Lplus=max(Lreq(gamma0/2), Lcarrier),
    Eplus=ceil(max(16,ell_old,512/h^2,sqrt(Lplus/(2r0)),
                  (8*Cfreeze/gamma0)^(1/6))),
    qplus=min(qstar/64,2^(-Eplus-2)).                              (2)

Here Lcarrier=2u/(gstar*nstar) may be retained if the additional signed
axial-carrier floor from leading_frame_foundation is used; it is not a
replacement for the full normal floor and can otherwise be zero. All the
existing normal/Riccati/cone/root bounds in Lreq must use the certified
frozen margin gamma0/2 (or a smaller one). A supplied Lreq at an unrelated
larger margin is not that budget. Every threshold is a finite VALUE
condition; none depends on a requested jet or nonlinear stage.

For ell>=Eplus, (1), Ls>=2r0 ell^2>=Lplus, and the parent's eighth-series
smallness bound hold. Enlarged charts meeting q<qplus have ell>Eplus and
4Q<=16q<qstar/4. Thus (2) repairs the stated acceptance contract. It is a
choice BEFORE realization/localization/Borel scales; it does not identify
the resulting field/force with either historical one. Any later domain
restriction or scale replacement still needs the parent's distinct-object
propagation. The historical first band and force are not recertified.

## 2. Normalize the source without hiding division by a tiny envelope

On a fixed accepted label P>0 on the entire closed pulse, P<=1. The
Gaussian/phase estimates supply all fixed jets of l=log P with polynomial
S bounds. For multi-indices let K_I>=sup|D^I l|, I!=0. Define nonnegative
Bell bounds B_0=1 and, choosing any e with I_e>0,

    B_I=sum_(J<=I-e) binom(I-e,J) K_(J+e) B_(I-e-J).              (3)

Differentiating P^(-1)=exp(-l), inductively,

    |P D^I(P^(-1))| <= B_I.

If |D^J f|<=F_J*w*P on the full path, w=sqrt(zeta) times the retained
inverse-edge weights, quotient differentiation yields

    |D^I(f/P)| <= w * sum_(J<=I) binom(I,J) F_J B_(I-J).           (4)

Choose common maximum edge exponents through the finite order under
consideration. Each B_I is a finite polynomial in finitely many K_J:
(4) adds polynomial S costs, not an exponential inverse-entrance cost.
It applies to vector components and then to -Bleft*Proj(f)/P, using the
current frame/normal floors. It neither differentiates frozen carriers nor
assumes a bound for arbitrary unweighted forcing.

Every generated regular source vanishes smoothly on the entrance collars:
linear regular expressions contain the smooth temporal cutoff; nonlinear
ones use COMPLETE cut velocities. Since P(0)>0 for each band, all mixed
entrance jets of f/P vanish. This supplies the exact entrance compatibility
needed for high-harmonic mixed pulse/parameter inverse bounds. Small P(0)
is not set to zero. Smooth support is different from Gaussian smallness.

The full pressure remains
pi_m=i(n.K*t_m-nprime.t_m+n.f_m)/(k*m*|n|^2). Pure normal f can have
t_m=0 but pi_m!=0, so normalization does not give pressure the extra two
harmonic gains of the potential. The complete Pi1/Pi2 background remains
in K and in the separate full mean/base equations.

## 3. Check the actual finite regular-source grammar

Ordinary coefficient derivatives of a complete wave preserve
epsilon^alpha sqrt(zeta) P times finite polynomial S/edge losses. In the
complete cylindrical linear remainder, every term is a coefficient jet of
the current t, complete curl remainder or pressure, multiplied by bounded
background/phase coefficients. Formula (4) therefore checks its full-path
normalized input. Radial advection, slow time, amplitude pressure gradient,
vector viscosity, mixed viscosity, radial multiplier derivatives and angular
connections are retained; their epsilon ledger is the reviewed regular-
source one, not a principal-amplitude substitution.

For two same-label wave factors, zeta*P^2<=sqrt(zeta)*P. Leibniz supplies
all fixed jets. The complete divergence identity, not principal
transversality, removes k from phase transport and leaves n/m; that finite
harmonic cost is retained. Zero outputs go to the complete mean equations,
never to an inverse with denominator m+n. Different labels use their
separated CLOSED enlarged supports, not a same-phase argument.

Mean maps may fill the auxiliary torus and need not retain P. Their actual
coefficient bounds times a wave still retain that wave's P, sqrt(zeta),
closed support and all entrance collars. A mean is not forwarded alone.
The physical mean-wave phase derivative retains its genuine half-power
loss. The actual compact A0/A1 velocity and all cutoff curl terms stay in
these products. Covariance averaging is bounded on each fixed coefficient
jet and removes P only by P<=1; it does not produce P^(-1).

Thus the finite regular expression grammar has the inverse's full-path,
edge and entrance inputs whenever its current mean/signed target jet and
moment bounds hold. There is no need for a separate arbitrary-source
inverse assumption. This is not an automatic verification of every later
target merely because an operation has the right name.

## 4. Quantitative primitives at both flat edges, with the moment retained

Put x=ya, Y-x=yb, a=t1^2>0, b=4, t1=C_setup^-20 UNCHANGED, and

    zeta=exp(-a/x^2-b/(Y-x)^2),
    g(x)=x^(-p)(Y-x)^(-q)*zeta, p,q>=0.                         (5)

For x<=dL, with

    dL=min(Y/2,1,sqrt(a/[2(p+1)]),(a*Y^3/[32(b+1)])^(1/3)),

direct differentiation gives
(log g)'=2a/x^3-2b/(Y-x)^3-p/x+q/(Y-x)>=a/x^3.
For 0<s<x this implies g(s)<=g(x) exp[-a(x-s)/x^3], hence

    integral_0^x g(s)ds <= (x^3/a) g(x).                         (6)

At the right edge interchange (x,a,p) with (Y-x,b,q). If a source f has
zero TOTAL radial moment and |f|<=C*g, its primitive has BOTH forms
F(x)=int_0^x f=-int_x^Y f. Therefore it has the original full zeta weight
at BOTH edges, with the extra three powers in (6). All fixed differentiated
primitives use finitely many differentiated input bounds; bounded positive
radial chart factors are included in C. The middle closed interval has a
positive weight floor, so the estimate extends over the whole shell with a
finite constant. The shifted compact mean primitive integrates toward each
edge and uses this same bound; auxiliary shifts do not change zeta.

The compact shifted T_e uses SEPARATE edge-directed pieces and an interior
A_e remainder. It does not require a pointwise zero radial moment merely
to be compact. In the zero auxiliary mode this is F=int_0^x f-chi(x)M,
M=int_0^Y f, with derivative f-chi'(x)M. Even for M!=0 this has zero
endpoint values with a nonzero interior A_e remainder. Its zero-Haar moment
is needed to make that remainder flat/solve the temporal mean. This is
distinct from the zero TOTAL needed for an ordinary pressure primitive
to have the same left and right formulas. Neither condition is dropped.

Without the zero total, a positive smooth f=zeta has F(Y)>0 and fails
every flat right-edge bound. This is a scalar support obstruction, not a
Navier--Stokes counterexample. Pressure/moment corrections and zero-Haar
conditions are therefore mathematical inputs, not optional bookkeeping.
Reserved mean bumps lie on a fixed compact interior hull: zeta has a
positive lower bound there, so their finite jet costs preserve this class.
No bound uniform in t1^-1, derivative order or lambda^-1 is inferred.

For any fixed integer p>=0 and c>0 one explicit edge absorption is

    sup_(x>0) x^(-p) exp(-c/x^2)
       <=max(1,ceil(p/2)! / c^ceil(p/2)).                        (7)

For x<=1 use the indicated exponential-series term; for x>=1 use x^-p<=1.
Derivatives add finitely many inverse powers. This proves smooth zero
extension with the actual half-exponential constants a/2 and b/2=2,
without changing the pulse width. Finite tests do not prove all orders.

## 5. Signed targets lose one half weight and recover the other half

On the supported INTERIOR the current canonical squared root satisfies
y_sigma>=c*zeta/s0 and every fixed jet is at most C_I*(zeta/s0) times
finite S/edge powers; s0=r0^2/sqrt(Ls) is the physical Haar scale. For
each term in a reciprocal-square-root derivative,

    D^I(y^(-1/2)) = sum C*y^(-1/2-n)*product_(nu=1)^n D^(Jnu)y,

the zeta power is -1/2-n+n=-1/2, not increasingly negative. A CURRENT
signed target with |D^I Xi|<=C_I epsilon^alpha zeta times edge powers,
and d=Hc^-1(Xi/epsilon), delta a=d/(2ac), therefore obeys

    |D^I delta a| <= C_I epsilon^(alpha-1) s0^(-1/2)
                         sqrt(zeta)*edge_powers.                (8)

Jets of the Haar factors are absorbed into finite polynomial S constants.
Both 1/epsilon and the physical r0 cost remain. No new positive root,
cutoff division or smallness of Xi is used. At the boundary the reciprocal
root alone diverges; its product in (8) has smooth zero extension by (7).
T0/zeta and ac/sqrt(zeta) need NOT be bounded: the actual outer pole
yb^-3 (and root pole yb^-3/2) are admitted by the finite edge powers.

## 6. Exact remaining gate and reproduction

These estimates supply the functional input checks for regular sources
generated by the reviewed finite operators, relative to their CURRENT
target bounds. The remaining concrete gate is to verify for the actual
recomputed signed/full-covariance mean targets, at each finite cycle,

    |D^I Xi| <= C_I epsilon^alpha zeta ya^(-pI) yb^(-qI),
    int R^2 bar(Etheta)=epsilon*partial_Z Jtheta,
    int R bar(Ez)=epsilon*partial_Z(Jz+c_rho*Pmean),               (9)

with current axial velocity, all six covariance entries, complete pressure,
the whole moving-product derivative, both velocity moments, zero-Haar data
and the reserved hull. The canonical finite-cycle checkpoint gives an
expression-level argument; (9) must not be certified solely by protocol
stamps or by a historical target from another family. Pi1/Pi2 may cancel
only against the SAME complete base. Flat base/cutoff residuals remain
additive: physical flatness does not imply a normalized source bound or a
flat forward response.

The new checks use exact mixed two-variable jets for exponential quotient
and reciprocal-root calculus, exact root-finding and cone-freezing controls,
edge derivative/primitive identities and moment obstruction controls. They
are not actual fields or continuum supremum certificates. Run
`python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
Fresh output is required. Historical source bindings and relevant inverse,
source, root and mean evidence are reproduced unchanged before publication.
