# Q133: final leading angular assembly and the order of spatial suprema

Parent: `efea38379f0d75edfe3b66f4a565e7196da0987f` (Q132).
Fresh fetch on 3 November 2026 agreed with this head. Q132 research run
41719704419 and field run 41719704534 were in progress with no failed job
steps reported. The first public API request returned 403; a retry and the
connected GitHub job reader supplied fresh results. Current notes, Q132
frozen evidence, Q127/Q126 and Q125's norm convention were read first.
Historical evidence, actual roots, supports, pressure gauge, N, physical
scales and C_setup^-20 width remain unchanged. The setup result stays local
leading core/first activation collar.

Q132 closes the same leading repair's conditional angular jets. This
checkpoint assembles a FINAL compact U0 bound, then joins it to the same
axis, taper and exterior regions. It also makes two DIFFERENT spatial
norm conventions explicit. The historical maximum rule is valid for
Q125's intended norm; its stronger reinterpretation is not. No historical
numeric result or physical counterexample is asserted to be false here.
Actual profile, cover, target and acceptance gates remain premises.

## 1. Final compact U0 from the SAME modulation and repair

Let C be the original whole leading compact region through the last repair
bump, including unmodulated subregions and gaps. The PREMOD cap u and eta
radius used in Q131 must bound the SAME native PREMOD U uniformly on C,
not only on the active modulation subinterval. This whole-cover binding is
an actual remaining input; a local bound does not establish it.

Use Q131's same axial modulation recipe: Pmod=u+MB/N and Bmod. Its
moving-phase primitive is extended by zero only using the original flat
collars, so the same angular bound holds on unmodulated parts as well.
Q132 gives the same repair amplitude bound Prep Brep^m(m!)^3 before the
fixed bump factor. Q63's conservative exact B_0=162*32600 bounds every
original unit-mass bump value. Supports are disjoint and independent of
eta, including flat endpoints, so no extra factor of five is introduced.
Thus for the FINAL field after BOTH operations,

    P=Pmod+B_0*Prep, B=max(Bmod,Brep)>=1,
    sup_C,eta |partial_eta^m U0| <= P B^m(m!)^3. (1)

The triangle bound is valid even without exploiting the original spatial
separation of modulation and repair. No new cutoff, root, support or
frequency is selected. The actual Q132 target must still include ALL five
baseline defects and incoming modulation differences, including pressure.
Neither ending the modulation nor U0=0 on a reserved patch resets them.

In the angular factorial SUM norm the resulting regional upper bound is

    C_m=P sum_(j=0)^m B^j(j!)^2.                (2)

The conversion divides each raw derivative by j!, leaving (j!)^2. Using
the raw m-th derivative bound directly as this sum is not the conversion.
An integer binary ceiling is

    b_C=ceilLog max(1,P)+ceilLog(m+1)
                     +m ceilLog B+2m ceilLog m, (3)

with the last term zero at m=0. This follows by bounding every term of
(2) by P B^m(m!)^2 and m!<=m^m. High order 100005 requires no factorial
or large power materialization. This is a conditional final compact
estimate, not a certificate that its whole-region premises have been
verified for the selected enormous profile.

## 2. The intended angular algebra norm and its valid regional maximum

Q124 defines |f(X,.)|_m=sum_(j<=m) sup_eta |partial_eta^j f(X,eta)|/j!.
Q125 explicitly asks for sup_[0,Xv]|U0|_(m+2): the spatial supremum is
OUTSIDE that angular sum. Write

    W_m=sup_X sum_(j=0)^m sup_eta |partial_eta^j U0(X,eta)|/j!.

For a cover by the same axis, compact, reference taper and exterior pulse
regions, ordinary locality of sup_X gives

    W_m <= max(Axis_m,C_m,Taper_m,Exterior_m).   (4)

This is the intended Q126/Q127 maximum rule. It is valid for W_m. With
Q127's actual leading axis Cauchy cap Kaxis and original delta, and Q126's
same exterior root contracts,

    Axis_m=Kaxis sum_(j=0)^m delta^-j,
    Taper_m<=8, Exterior_m=40(2^(m+1)-1).        (5)

All later leading axial fields vanish. Equation (2) now supplies the
separately finalized compact input missing in Q127, CONDITIONAL on the
original cover and all Q128--Q132 premises. A PREMOD bound alone remains
insufficient, and the helper rejects it. The global function also checks
order and raw-versus-binary conventions. It does not validate an identity
label or numeric cap merely because it was passed to a calculator.

For Q125's discrepancy order m, the interface needs W_(m+2), not W_m.
The first_input helper checks that order and returns the same angular
algebra norm, with no raw-factorial relabeling. The actual F0 norm on its
whole support, core/collar inputs, full pressure identities and support/
scale contracts are STILL needed before that first-discrepancy estimate
or any physical selector can be certified.

## 3. A stronger norm: a precise obstruction to swapping sup and sum

A stronger quantity, useful to distinguish during raw derivative work, is

    S_m=sum_(j=0)^m sup_X,eta |partial_eta^j U0|/j!.

Always W_m<=S_m, but the same upper bound for W_m cannot simply be
reinterpreted as a bound for S_m. For S_m the general regional gluing is

    S_m <= sum_regions S_region,m,
    NOT in general max_regions S_region,m.     (6)

This does not invalidate the historical intended W_m inference in (4).
The source check clarified this distinction before publication; there is
no claim that the actual Navier-Stokes construction realizes the following
manufactured counterexample.

Take the smooth field U(X,eta)=100[1-X+(X/4)sin(4eta)] on [0,1]x[-1,1],
split at X=1/2. Since sin attains 1 at eta=pi/8 and cos attains 1 at 0,
for fixed X the value and first angular derivative suprema are 100-75X
and 100X. Hence W_1=125; the two regional W_1 values are 112.5 and 125.
Their maximum is correctly 125. In contrast,

    S_left,1=100+50=150,
    S_right,1=62.5+100=162.5,
    S_whole,1=100+100=200 >162.5.                (7)

Different derivative orders attain their spatial suprema in different
regions. Multiplying a radial interpolation by smooth monotone flat
weights yields the same example with flat endpoints and a zero extension;
the obstruction is a norm-order issue, not a jump or nonsmooth gluing.
At m=1 the inherited taper/exterior padding 8 and 120 would not repair
this stronger maximum inference: max(150,162.5,8,120)<200. The old helper's
arithmetic test of a maximum does not prove an S_m statement. Nor does
this generic counterexample refute its intended W_m statement.

When the stronger convention IS required, the same Cauchy and Gevrey
caps control each regional S_m as well. A safe upper bound is

    S_m <= Axis_m+C_m+8+Exterior_m.             (8)

A sharper, moderate-order alternative is the orderwise maximum

    sum_(j=0)^m max(Kaxis delta^-j,
             P B^j(j!)^2, taper_j, 40*2^j),     (9)

where taper_0=taper_1=4 and higher entries vanish. Equation (9) is at most
(8). For binary budgets (8) costs at most TWO extra bits beyond the
maximum of the four regional ceilings. That extra reserve is needed only
for this stronger convention; Q125's intended W_m uses the maximum.
Neither corrected norm routing nor a finite numeric bound revalidates an
old retained-scale inequality automatically. Every requested scale gate
must use the norm/order convention it actually requires.

## 4. Checks, scope and remaining obligation

Exact smooth controls verify (7), both sup/sum orders and the inherited
maximum's proper W_m interpretation. Independent finite jets of a
modulated-plus-repaired field check (1), and positive factorial sums check
(2),(3). Same-stage/order/convention controls reject PREMOD substitution
and an m-versus-m+2 source input. Moderate-order sum and orderwise bounds
are compared; high binary budgets avoid factorial expansion. These are
manufactured diagnostics, not sampling of the selected profile. The audit
replays Q132 and checks the frozen source/evidence hash closure; fresh
output must match byte for byte. Passing tests are not a continuum proof.

Next prove the actual whole-compact PREMOD caps, zero extensions, same
accepted target/root/positive-field gates and complete regional cover.
Then apply the intended W_(m+2) input with the actual F0/core/collar data,
and check the low-order state/cone margins and retained scales. Even S_m
controls only pure angular derivatives: Q121 requires RAW MIXED X,eta
profile bounds, with U through m+1 over the whole axis primitive interval.
No radial derivative estimate follows by relabeling (2),(4) or (8).
Global coefficient/reserved-mean supports, uniform actual physical pulse
bounds, full wave-covariance mean corrections and complete viscous/cutoff/
interaction residual remain open. No unconditional pulse, full residual
improvement or blow-up follows. All historical frozen records are preserved.
