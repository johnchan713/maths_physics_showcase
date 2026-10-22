# The original length selector already closes normal nonvanishing

Parent: `5ba2fad301fac1ac7a68885c0d79d70caa7cdb76`.
This checkpoint closes the NORMAL NONVANISHING implication for the original
rule-defined accepted primary tail, relative to its retained valid primitive
background/reference bounds. It proves that the existing Riccati length
dominates the existing normal threshold. No separate Ctan or Cphi
compatibility inequality is needed for that conclusion. Nothing is selected
again. Unknown adherence of an arbitrary older numerical band remains a
different issue. The full propagator, covariance and residual are not closed.

## 1. Use actual pulse length, keeping complete axial transport

Use the SAME reference and original phase as `frame_column_foundation`,
section 1. Write b_-=bmin<=Bs<=bmax=b_+, and let R,R0 be in the original
[Rmin,Rmax], |R-R0|<=Dbox S^-3, 1<=L<=S. The unrounded carriers obey

    ptilde/R0=Bs[Ktheta-sigma*u*gtheta/(L|g|^2)],
    pz=Bs[Kz-sigma*u*gz/(L|g|^2)], |p-ptilde|<=1/k.

The original tangent normal is (p/R,pz-epsilon*v*H_Z),
H=pF+pzG, 0<=v<=L. Subtract Bs K without estimating the whole radial
normal. In particular the angular identity is

    p/R-Bs Ktheta = (p-ptilde)/R
       +Bs Ktheta*(R0/R-1)
       -sigma*Bs*u*R0*gtheta/(L*R*|g|^2).                     (1)

This separates box displacement from tilt; a factor u/gmin need not be
charged to the box term. With the SAME valid Mz bounding BOTH |F_Z| and
|G_Z| of the COMPLETE background, define the original constants

    Pmax=Rmax*b_+(1+u/gmin)+1, Zmax=b_+(1+u/gmin),
    Zc=(Pmax+Zmax)Mz,
    Dt=b_+ Dbox/Rmin+b_+ u Rmax/(Rmin gmin)
       +1/Rmin+b_+ u/gmin+Zc.                                (2)

The extra +1 in Pmax retains nearest-NONZERO rounding. Equations (1)–(2)
give a component-sum bound, hence also a Euclidean bound,

    |n_tan-Bs K| <= b_+ Dbox/(Rmin S^3)
       +b_+ u Rmax/(Rmin gmin L)+1/(Rmin k)
       +b_+ u/(gmin L)+epsilon*L*Zc.                         (3)

Under the ORIGINAL all-path slow gate S^2(epsilon+epsilon^2+k^-1)<=1,
each term is bounded by its own coefficient divided by L, since 1<=L<=S.
Thus |n_tan-Bs K|<=Dt/L. No pressure reset, leading-only H_Z, omitted
axial velocity, cutoff curl omission, or radius-free replacement L=S is
used. The geometric proof still supplies the actual L=2r0/ci.

It follows that L>=Da:=2Dt/b_- gives |n_tan|>=b_-/2 and

    |n|^-2 <= 4/b_-^2.                                      (4)

This is the ORIGINAL reference-Bs floor. The linear quantity called beta
in `frame_budget.py` is b_-/2; it must not be confused with the different
squared phase quantity beta_phase=lambda_min/[4(1+U^2)^2]. A coarse old
b_- may have b_-^2<beta_phase. Then (4) does NOT prove the stronger
4/beta_phase upper bound. Downstream estimates demanding that particular
stronger constant must still justify it or explicitly use (4) in their
own budgets. Nonvanishing itself needs only (4).

## 2. The retained Riccati length dominates the normal threshold

Use the original positive formulas in `frame_column_foundation/frame_budget.py`:

    sr=3u/2, S0=1+sr, Da=2Dt/b_-, Ds=(Dr+sr*Dt)/(b_-/2),
    DH=D11+D12+D21,
    CE=kM*DH+S0*S1+lambda_max*Ds,
    mu=lambda_min/(1+sr), LR=8CE/mu, kM>=1.                  (5)

All these are functions of primitive input BOUNDS, before any ODE solve.
The existing generator comparison is not invoked to prove its own normal
denominator. We first establish the denominator using the positive budget
formulas, and only then may those generator/Riccati estimates be applied.
This avoids an apparent circular bootstrap.

The actual leading growth formula gives

    lambda0^2=2F0*(-Ntheta)(|g0|-2F0*(-Ntheta))
       <=2F0|g0|.

Therefore VALID bounds on the nonempty reference family satisfy

    lambda_min^2<=2 Fmax gmax,
    lambda_min<=Fmax+gmax.                                  (6)

These are consequences of the retained reference formula, not an added
spectral assumption. Arbitrary positive diagnostic inputs need not obey
(6); the new verifier checks it explicitly.

The SAME budget formulas give, using only nonnegative terms,

    DH >=Da[4Fmax+(1+S0)gmax]
       >=2Da(Fmax+gmax) >=2Da*lambda_min,
    CE>=DH, CE>=lambda_max*Ds>=lambda_min*Ds.

Consequently

    LR>=16(1+sr)Da>=Da,
    LR>=8(1+sr)Ds>=Ds.                                     (7)

The original normal threshold is Ln=max(1,Da,Ds). Thus

    max(25,LR)>=Ln.                                         (8)

This is a uniform analytic inequality at the actual finite parameter
values, however large the constants. No numerical giant field, new
cutoff norm or new independent floor assumption is needed to establish it.

## 3. Join the SAME original selector without reselecting anything

`core_primary_tail_recheck`, section 1, defines Lreq to include the original
Riccati and normal thresholds, along with 25 and its cone/comparison
lengths. Its retained max selector guarantees, on EVERY accepted band,

    2r0 ell^2>=Lreq>=max(25,LR),
    Cslow ell^4(epsilon+epsilon^2+k^-1)<=tau, Cslow>=1,tau<=1.

The original geometry gives 2r0 S<=L<2r0*Tg*S<S. Hence 1<=L<=S,
the slow gate used in (3) holds, and (8) proves (4) on the rule-defined
accepted tail. The normal threshold was even explicitly included in Lreq;
(7) shows it is already dominated by its Riccati part for valid inputs.

This closes that original denominator gate relative to the SAME primitive
bounds. The later 36 Cphi^2<=beta_phase*S^2 and 4 Ctan^2<=beta_phase*S^2
tests are sufficient alternatives with their own quantitative constants.
Their separate compatibility margins are not necessary conditions for the
original positive normal denominator. Failing one of those stronger recipes
does not negate (4)–(8).

`verify_existing` checks a SUPPLIED old band and old length input, using
the original band inequalities. It never invokes select_tail and never
returns a replacement band or length. An independently fixed historical
band or length with unknown semantic provenance still has the exact
remaining checks: its original band inequalities, valid primitive bounds,
and Lreq_existing>=max(25,8CE/mu). It is not retroactively declared accepted.
Lreq>=Ln is a weaker sufficient option when only this floor is needed.

The selector's fixed dictionary input `epsilon` is used by its original
tilt/deficit certificate; it is not substituted for the band's Q^h. The
band slow gate is checked through the original all-half-line series bound.

No larger later bound for Mz, Db or CE is silently substituted into an old
Lreq. The original Mz can be justified from its actual suprema or from a
valid full-background first-jet estimate; the preceding q^h diagonal
bound supplies one such finite estimate. A newly conservative numerical
upper enclosure need not be dominated by an independently fixed old length.

Removed independent gate: an extra Ctan/Cphi selector compatibility test
for original normal NONVANISHING on the original rule-defined accepted
tail. Added analytical assumptions: none. This route uses the already
retained frame/Riccati primitive budgets, including Db, whose actual
validity is not newly proved here. The independent low-order Ctan route
still applies without those full frame assumptions. Core/coefficient
regularity, actual primitive bounds, original physical chart/support and
q-domain gates remain prerequisites of their respective accepted domains.

Next obligation: propagate the original normal bound through the sourced
inverse with its complete pressure and actual denominator constant, then
close the generated complete covariance/mean and residual inputs. Keep
any unknown historical-band membership separate from this rule-defined
tail; do not reopen its normal nonvanishing as an independent premise.

## 4. Reproduction

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output fresh.json`. Checks cover the original exact
budgets, both tilt terms, full axial transport, algebraic physical
consistency, the strict factors in (7), unchanged-band success/failure,
the distinction between b_- and beta_phase, and frozen scope/source
bindings. Manufactured rational controls are not actual profile enclosures,
continuum proofs, unconditional pulses, complete residual bounds or blow-up.
