# Q40: derive the current row norm before the next fast inverse

Parent: `37bd02fc95254cacb2d59fbf2bf5d3866e5abcad` (Q39).
Parent focused CI 41140974189 passed; broad CI 41140974132 was still
running at inspection on 30 October 2026. Q39's independent new-fast-field
norm slot is now supplied by an explicit current-row/compact-inverse
composition. Its initial A1-tail, compact/time-inverse and original
background contracts remain conditional. No frozen evidence is altered.

## 1. The current source is a sum of retained contributions

Let F0,S0 be the original initial fast and slow means. The real signed
wave changes no mean velocities. The inherited exact initial mean identity,
followed by the signed covariance update, gives

    Qcirc E_current = Qcirc [Rf(0,F0)+Rs(F0,S0)
                       +(0,A1_initial Fz_initial)+DeltaE_wave]. (1)

Here Qcirc=I-Haar has norm at most two, and Rf/Rs omit only the tangential
c Nt term cancelled by the original fast solve (or zero for a slow field).
The initial Haar row disappears under Qcirc. Its axial A1 remainder does
NOT generally disappear: it is kept explicitly in (1). Neither the initial
fast velocity nor its A1 velocity component is discarded.

The implementation obtains F0 and S0 from the SAME Q35 initial budget,
using weighted slow fields on the entire original potential hull. It calls
the complete historical difference functional twice, first from zero to
F0, then from F0 to F0+S0. These functionals retain graph/cylindrical
products, slow time, vector viscosity, radial pressure source and its
projected primitive. They are not just quadratic velocity differences.

DeltaE_wave is computed from Q37's COMPLETE covariance change, all six
entries, including principal crosses and signed squares. Its radial source
and rho-compensated pressure are recomputed with Q38's new-pole pressure
map and whole-hull moment bound. No selected-flux cancellation is used to
replace the full covariance norm by its smaller error norm.

At row order n, all field/covariance/fixed norms are requested through
N=n+3. The output positive row bounds are

    Rtheta_n=2(Rf_theta,n+Rs_theta,n+DeltaE_theta,n),
    Rz_n=2(Rf_z,n+Rs_z,n+DeltaE_z,n+A1initial_n).              (2)

There is no independent current row callback in this construction.

## 2. One common row pole and separate initial estimates

For a finite maximum row index M, set Nmax=M+3, let p0 be the initial
budget's common pole and d the Q36 signed-target allowance, and take

    q=max(6p0,6Nmax+9,6Nmax+d+6), v=2q.                     (3)

The retained old full-row contracts fit the conservative 12p0 allowance;
the complete covariance at every requested index fits v. Lower-pole
terms embed without extra factors because the retained delta is capped
by one. Old fixed input recipes are still requested at their own p0;
pressure and moment constants are freshly evaluated at v. They are not
inferred by monotonicity from an old pressure constant.

The initial A1-tail norm must declare order at least n and pole at most v.
It must bound the SAME initial remainder in the weighted norm. This check
is only domain bookkeeping: it does not certify the original raw cutoff
kernel's high derivatives or the epsilon rate of a user-supplied number.
The full source/flat-extension and background contracts remain essential.

With B=1-k, alpha=2-3k and D=alpha-k, the original fast-row difference has
order alpha and the initial slow difference is no worse. The signed full
covariance row response has order D. Thus (2) has order D provided the
retained initial A1 remainder has order at least D. The original compact
integration-by-parts recipe must supply that inequality; no zero tail is
assumed. These are first-real rates, not later generic/complex rates.

## 3. Apply the complete fast inverse with the correct derivative orders

For fast output j, use current rows at j+4 for tangential/axial fields and
j+5 for the radial field. Let I bound c_time^-1/S on the unchanged graph,
KT be the compact primitive bound and KA the axial compact correction.
The historical inverse yields sufficient coefficients

    Cj=binom(j+5,5) 24 2^(j+4) (I*S) (j+4)!,
    u_new <= Cj Rtheta_(j+4),
    z_new <= Cj (1+KA_j) Rz_(j+4),
    x_new <= binom(j+5,5) 24 2^(j+5) (I*S) (j+5)!
                                      KT_(j+1) Rz_(j+5).    (4)

The polynomial S is retained explicitly. Omitting it would misread the
scaled time-floor input. KT is requested at j+1, not j; KA at j. The
factor 1+KA retains the actual axial gd-A1 gd velocity. Bounds (4) have
conditional orders D,D,D+1 and common pole v. The compact/time provider
must declare coverage at both orders and at v. Q38's radial pressure map
alone is not a proof of these oscillatory compact-map estimates.

`FastFields` binds (4) to the derived `CurrentRows` object and returns
exactly Q39's field/order/pole interface. It can therefore replace Q39's
independent new-fast norm callback. For Q39 slow output through m, reserve
new fast through m+2, current rows through m+7, Q37 pressure ledger through
m+8 (covering its covariance through m+10), Q36 signed target through
m+11, and its profile/primary covariance through m+24. These are finite
composition ceilings before canonical/background expansion and before the
much larger compact-tail graph. Low-order fixtures do not cover that graph.

## 4. Validation and remaining gate

The tests evaluate the current-row decomposition using the actual Q36/Q37
norm chain with manufactured lower inputs. They detect removal of the
initial A1 remainder, preserve both mean differences and their viscosity
and pressure terms, and check the original root and covariance pole.
Independent exact fast-factor controls verify the S factor, KT(j+1),
KA(j), j+4/j+5 calls and distinct radial/tangential constants. The high-order
fast composition is an analytical finite recipe, not a numerically
materialized continuum field or full compact-tail certificate.

The remaining inequalities are the SAME initial A1 weighted-tail estimate,
original compact KT/KA and time-floor estimates at (3), and the actual
background/current input contracts throughout the nested request graph.
The code does not infer those from domain declarations. The postslow
scalar reset and complete final row differences still need composition,
followed by the unaveraged wave/mean source with whole pressure memory.
Uniform actual pulse bounds, complete viscosity/cutoff residual improvement
and blow-up remain unproved.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json` here.
