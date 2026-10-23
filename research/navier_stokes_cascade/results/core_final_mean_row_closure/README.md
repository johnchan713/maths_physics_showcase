# Final regular mean rows: fresh cancellation and the pointwise A0 estimate

Parent: `cbb0f3f352bbf8d83671c4adfdfdea96bbe0ba23`.
Its localized canonical stress check passed in hosted CI when inspected;
the overall workflow was still running without a failed step. The parent
proved complete regular tangential mean rows at D=3/2-2kappa at the
POSTWAVE stage. This checkpoint follows the actual fresh fast/slow means
and bounds the final pointwise radial compact remainder. Relative to the
original inputs, all THREE final regular angular-zero rows have order
C=3/2-kappa. The exact additive physical tail ledger is still separate.

This is a conditional analytical implication on the same accepted bands,
width C_setup^-20, canonical roots and fixed profiles. It independently
evaluates neither the actual enormous fields nor their high-jet bounds.
Historical sources and frozen records are unchanged. No new background,
cone, operator, smallness or band-threshold assumption is introduced.

## 1. Fresh means change the row: use their cancellation, not their size alone

Use B=1-kappa, C=B+1/2, D=C-kappa and 0<kappa<1/4. Let Ftheta,Fz
be the Q1=I-Haar parts of the postwave tangential rows. The preceding
fast-row checkpoint derives F in M_D and the SAME actual fresh fast
increment (X,U,Z) with tangential size D, radial size D+1. Its exact
identities, including the compact axial velocity, are

    c Nt U=-Ftheta,
    c Nt Z=-Fz+A1 Fz.                                  (1)

Here Z is gamma_d-A1 gamma_d, not just gamma_d. These increments and
their radial companion have zero auxiliary Haar average and preserve
both velocity moments. The fresh slow increment has tangential size C,
radial size C+1 and is auxiliary-independent on its entire prescribed
potential hull. The full current G=0 and same-V-power-law reserved-hull
contract used to derive that slow estimate remains an original premise.

Subtract the COMPLETE mean equations at each operation with the same
wave and base. `core_fast_row_gain_join.difference` gives the entire
remainder excluding only c Nt of its tangential increment. In particular
the pressure difference is exactly

    dpm=T0(dgr-rho dP),  dP=int bar dgr dR,              (2)

with unchanged moving rho, both full base pressure tails, all covariance,
radial/axial flux, connection and viscosity terms. The common base tails
cancel only in the subtraction; they remain in the full residual ledger.
Actual old means are coarse B/B+1, not set to zero. Applying the same
positive polynomial to a fast increment at D and slow increment at C gives

    Rfast in M_(D+1-2kappa),
    Rslow in M_(C+1-2kappa).                            (3)

The smallest displayed order is from epsilon times two graph derivatives
in viscosity. A pressure difference itself is only of order D or C;
its axial derivative contributes epsilon. A radial transport uses either
b=epsilon bN or a radial mean. Old/new mean products have an additional
positive B factor. Slow c Nt vanishes exactly. Both exponents in (3)
are greater than C, since 1-3kappa>0. No graph coefficient is omitted.

Combining (1) with exact subtraction, AFTER both fresh maps,

    Q1 Efinal=(0,A1 Fz)+Q1 Rfast+Q1 Rslow.              (4)

Do not carry the old F a second time in (4) or identify the postwave
and final rows. Pressure in each row is reconstructed from that row's
own full state. The code's sharp fast size is its derived FastRowGain
recipe, not a coarse B norm relabeled D.

## 2. The retained compact estimate gives the fresh A1 error

For output raw coefficient order n, integrate the SAME whole-line
compact operator p times by parts. Put J=n+2p+3. For e=0 or 1 its
original estimate is

    N_n(Ae Q1 f) <=
      2^(n+J) L Bn AJ [8*3^p*2^(n+p+3)]
                       M^-p N_J(Q1 f),                 (5)
    Bn=N_n(k_e chi_U),  AJ=N_J(f_e),
    M^-p < c_i^(-rho_g p) epsilon^(p kappa).

The original irrational radial divisor is bounded below by 1/(3|k|_1).
Compact smooth zero extension kills every boundary term. Fourier lattice
summability and frozen-M coefficient differentiation are exactly those of
`core_fast_inverse_budget` (8). Derivatives of the cutoff/profile factors
are included. A generic Ae operator norm does not bound their product.
The original positive profile recipes do. For through-order output norm,
Q1 costs at most two, 6^n safely converts raw output jets in five
variables, and raw N_J is bounded by J! times the corresponding through
norm. The implementation records these factors, powers and factorials
without allocating their giant selected integer values.

For the fresh A1 Fz in (4), p=2 suffices:

    A1 Fz in M_(D+2kappa)=M_(3/2).

Since 0<kappa<1/4, (3) has still higher order. Thus the final Q1
tangential rows have order 3/2>C. The unchanged localized principal
stress identity and complete final Haar bound of the parents give

    bar Efinal_theta, bar Efinal_z in M_C,
    Efinal_theta, Efinal_z in M_C.                       (6)

There is no new large integration count for THIS fresh A1 application:
J=n+7. Its postwave source already contains the PREVIOUS A1 recipe,
which must remain nested. At kappa=1/100000, output zero reads its old
source through 200008. The fast radial input for (3) needs old source
through 200009 and original initial means through 200012. These finite
dependencies, already present in the next-source ledger, are not truncated.

## 3. The exact zero mode of the radial A0 remainder

At the FINAL state let gr be its full graph radial source, including all
six complete covariance entries and actual mean products. Define

    Pnew=int bar gr dR, h=gr-rho Pnew, pm=T0 h,
    int rho dR=1, Dr T0 h=h-A0 h.

The base-subtracted regular radial row is EXACTLY

    Er=-rho Pnew-A0 h.                                 (7)

The parent derived Pnew in M_(B+C) from the actual fresh slow reset,
retaining the whole moving axial-pressure product. Its Haar statement
does not alone bound A0 pointwise. The additional exact fact is

    int bar h dR=0,
    A0 bar h=0,
    A0 h=A0 Q1 h=A0 Q1 gr.                             (8)

For e=0, f_0(U)dU=dR, so the zero Fourier coefficient in the whole-line
integral is precisely the moment in (8), not a different U moment.
The moving compensator is auxiliary-independent. Its full derivatives
are retained; the cancellation holds for each slow parameter before
differentiation. It is an integrated cancellation, not bar h=0 pointwise.
If rho is not unit mass, an M-independent term proportional to
(1-int rho)Pnew remains. It cannot be called oscillatory or Gaussian-flat.

The complete gr has the retained coarse order B: a graph derivative of
covariance has order 1-kappa=B, and V times the coarse tangential mean
also has order B. Radial mean time/transport and viscosity terms are
stronger. No pressure inverse or temporal inverse is needed to bound gr.
Apply (5) with e=0 to Q1 gr, choosing

    p0=floor(1/(2kappa))+1,
    B+p0 kappa>C, J0=n+2p0+3.                         (9)

The strict exponent margin preserves the same band domain, with finite
positive polynomial constants. Combining (7)--(9) with B+C>C gives

    Er in M_C.                                         (10)

This is a POINTWISE auxiliary bound for the regular radial mean row.
It does not promote pm, bar pm or the actual mean velocities to order C.
For kappa=1/100000, p0=50001 and J0=n+100005. Output zero therefore
requires the CURRENT final gr through raw order 100005, obtained from
current fields through at least 100008. The underlying original mean
constructor entails further finite shifts. These actual inputs remain
unverified, as do the already larger nested A1 dependencies above.

## 4. Weights, domains and the resulting limited closure

Use a common original finite pole p covering EVERY finite jet required
in this composition, including (9) and the nested old A1 inputs. Do not
reuse a low-order pole when adding high jets or infer band-uniform
seminorm bounds from smoothness alone. The established reweighting puts
mean inputs at 6p and their products/radial source at 12p. Operator
domains at 6p,12p and the estimate argument S+K remain unchanged.

For (5), first unweight the source with Eedge(12p,original width). The
output is on the ORIGINAL interior chi_U hull, with its actual gap. The
weight-floor inverse is bounded by the finite recipe
3^ceil((width^2+4)/gap^2), and c_i^-rho_g p is bounded by the original
time-floor polynomial to power ceil(rho_g p). The raw e=0 kernel
profiles must be those of A0; A1's operator constant cannot substitute
for them. Existing lower-pole embeddings retain max(1,L)^ceil(12p).
The deferred `current_radial_source_polynomial(J0)` applies that embedding
to the actual final fields before the complete graph-source products;
it does not relabel an unembedded coarse row norm.
All are finite profile-dependent constants, not physical scale changes.
No derivative-dependent first band or materialization requirement is added.

Equations (6),(10) close all three regular angular-zero rows at C for
this same FINAL first-sourced prefix. The retained multiharmonic nonzero
regular source remains at D, so the combined regular ledger has order D.
This is NOT a bound for the entire physical residual: the independent
base residual and exact additive Gaussian, spatial and scale cutoff tails,
wave pressures and all cross interactions must still be joined with their
physical losses. No additive tail is passed to the next sourced inverse.

Independent finite controls check moving unit-mass compensation, the
nonzero pointwise compact A0 mode, exact compact pressure endpoints,
nonunit-mass obstruction, the full pressure-difference and viscosity rate
budgets, and the original moderate Fourier constant. They use finite C5
polynomial data and rational formal endpoint phases; they do not certify
the original C-infinity profiles or their high jets. Moderate composition
controls retain the exact large-order RECIPES. Some even moderate-order
coefficient integers are large, so references to their exact existing
positive-polynomial formulas are retained rather than serialized decimals.

Next join the exact additive tail ledger and every physical scale loss
to this regular closure. Independently discharge actual background/cone,
root/profile, high-jet, support and endpoint contracts before claiming
uniform actual pulses, complete physical residual improvement, an infinite
field or blow-up. No such claim follows from these tests.

```sh
python3 research/navier_stokes_cascade/results/core_final_mean_row_closure/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_final_mean_row_closure/audit.py \
 --verify-record research/navier_stokes_cascade/results/core_final_mean_row_closure/evidence.json \
 --output /tmp/core_final_mean_row_closure_fresh.json
```
