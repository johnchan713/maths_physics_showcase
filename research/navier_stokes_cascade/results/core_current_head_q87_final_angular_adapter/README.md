# Q87: the final angular norm from exact fast cancellation

Parent: `b5b1844eede2160ae1f4588db32b5982bf68f058` (Q86).
This supplies Q84's remaining e coefficient slot on the original first-real
prefix, relative to the corrected factory and its physical rate/domain
contracts. It does not compute a continuum supremum or a complete residual.

Q86 correctly rejects using an absolute final-field norm or the projected
INCOMING row as the final angular norm. Its suggested use of Q52's finer
leading-viscosity decomposition is unnecessary for this particular, coarser
slot. The old H=2−3kappa bound concerns the earlier fast difference that
GENERATES the new source; it is not the rate of the FINAL new-fast row
difference. We preserve Q86 and document the narrower route here.

## 1. Derive the Haar term on the same pressure state

Use exactly Q86.FinalMeans and its corrected Q47 Postslow object. Write
O for the full original fast+slow means, F for the complete new fast field,
and S for the refreshed slow field. The wave covariance is fixed during
these mean updates. The post-signed angular row Etheta^s changes from the
initial row only through the complete signed covariance change Delta W_s.
Its Haar average obeys exactly

    Haar Etheta^s = Haar Etheta_initial
          +(partial_R+2/R) Haar Delta W_s,rt
          +epsilon partial_Z Haar Delta W_s,tz.               (1)

There is no current angular pressure gradient. Haar removes the auxiliary
graph derivatives only AFTER the products are formed. It does not remove
fast-fast correlations in the initial row, which are already included in
the actual InitialHaar.rows coefficient program.

The initial Haar coefficient is retrieved from the SAME signed target's
source.haar.rows(j), at its proved C0=2−2kappa rate. Q85's same-root signed
coefficient (excluding the original curl square) bounds Delta W_s at rate
TWO on Q51's target premise. Thus both derivative terms in (1) are at
least C0. The radial derivative costs a slow derivative of Delta W_s, not
a graph epsilon loss. At output j a sufficient positive signed Haar
coefficient is

    ((j+1)(g+1)+2r) signed_(j+1),                            (2)

where r bounds 1/R and g=dU/dR. The epsilon in the axial derivative was
safely dropped for this coarse rate. The signed coefficient already keeps
both orientations, signed squares, full curl and actual cutoffs.

Unweight the initial Haar row at its own pole 2p and the signed change at
its OWN returned covariance pole, using the inherited fixed shell factors.
Then add (1). This is not promotion between unrelated weighted norms. The
same whole positive hull and accepted fixed-q domain must apply to both.
No original or later pressure trace is reset to obtain (1).

## 2. Cancel the new fast term before estimating the two differences

In the regular mean system the complete first-real angular identity is

    Etheta_final = Haar Etheta^s + Rtheta(O,F)
                                      + Rtheta(O+F,S).        (3)

Indeed c_i Nt Ftheta cancels minus the zero-Haar incoming row exactly;
c_i Nt Stheta=0. F includes its axial gd−A1 gd velocity in every product.
The refreshed slow field is computed from the current postfast scalars,
including the signed covariance change ONCE and full old/new fast products.
This is the same pressure state used to construct (1), not a sourced complex
state or a previous final row.

For old components (x,u,z), increment (X,U,Z), and full base (b,V,G),
the existing angular difference polynomial at output j is

    rt=bU+VX+xU+Xu+XU,
    at=GU+VZ+zU+Zu+ZU,
    Rtheta* = dU+(g d+2r)rt+d at+Vis U,
    d=j+4, Vis=d²(g²+1)+r g d+r².

In the implementation g in (4) is padded to 1+2g, using |v_r|_1<2, so it dominates the
slow derivative as well as the graph part even if the supplied graph
multiplier bound is less than one.                           (4)

It is the ANGULAR component of the full difference in
core_fast_row_gain_join, with its two-graph radial viscosity and axial
viscosity retained. No unverified full-torus Haar-pressure norm enters (4).
The complete pressure remains in the fields' upstream construction.

For B=1−kappa, new-fast order D=2−4kappa and refreshed-slow order
alpha=2−3kappa, the worst difference orders from (4) are respectively
D+1−2kappa=3−6kappa and alpha+1−2kappa=3−5kappa. Other time,
transport and full old/new product terms are no worse. Old O+F stays at B,
including its radial companion B+1. For 0<kappa<1/4, both displayed rates
are at least C0=2−2kappa. Hence (1)--(4) give a positive polynomial bound
for Etheta_final/epsilon^C0 at every fixed requested j.

This does NOT raise the entire angular residual to a finer rate than Q52.
The Haar estimate in (1) is deliberately coarse. Equation (3), not an
absolute estimate of c_i Nt v_final, is the reason C0 is available.

## 3. Executable coefficient adapter and derivative requests

The adapter takes only Q86's corrected object. It retrieves the initial
Haar row, Q85 signed covariance, full old means, actual complete fast fields,
refreshed slow fields and fixed multipliers from that object's constructors.
There is no independent outgoing e or anonymous fine-row callback. It checks
that the pressure/covariance object identity remains unchanged after creation.
Q86's shape/domain/factory reserves and Q85's rate-2 premise are retained.

For angular output j, (4) uses mean and fixed data through j+3; (2) uses
signed covariance through j+1; InitialHaar.rows(j) retains its own previous
mean/covariance/pressure requests. We keep factory reserve j+3 rather than
assert a new minimal global primitive ceiling. The returned coefficient is
unweighted in the factories' R chart. It is NOT yet a native-U coefficient.
For Q84 the adapter's native_coefficient applies the finite coordinate-jet
proof of Q43 on the SAME hull. For 0<1/d<3 the inverse-power coordinate
R=U^(1/d) has normalized positive-degree jets bounded by
Rmax binom(k+2,2) Umin^(-k), since |binom(1/d,k)|<= (3)_k/k!.
Add four identity-coordinate jets at degree one to form A(t). A sufficient
normalized-R-input to raw-U-output factor is
j! sum_(k<=j)[t^k](1-A(t))^-1. Multiply additionally by 243 to bound the
factorial-through U norm. Only finite positive sums occur. The supplied
Umin must enclose the native image of this very same positive hull; positivity
alone does not certify that premise. No derivative or epsilon power is lost.
This makes explicit a chart conversion that must also accompany Q85/Q86's
R-chart coefficient slots before any subsequent native-U composition; a
norm-factor conversion by 243 alone is not a radial chart conversion. A harmless factor 243 additionally pads any inherited raw-to-
factorial through-order conversion in the five coefficient variables, as
in Q85/Q86; already-converted terms may be overestimated by this padding.

To supply Q84 at radial output n set j=n+200000. The angular constructor
then needs factory reserve n+200003, ONE more than Q86's final-field slot
alone. The same builder's conservative profile/original-covariance reserve
is n+200030. This is not the total depth of the archived nested A0/A1 DAG.
No derivative-dependent first band, new physical cutoff or new root is chosen.

Together Q85, Q86 and this adapter give the covariance, defect, mean and
angular coefficient slots in Q84 relative to their original physical
contracts. Base, temporal and common weighted-domain constants still need
the same-hull supplied bounds. The rational Q84 control evaluator does not
certify EP-valued actual norm inputs merely by evaluating a fixture.

Next compose these slots and the original base/temporal/operator records
on ONE archived pole/hull, and check the regular and additive residual
ledgers on that same prefix. In particular do not pair first-real compact
tails with a later sourced regular residual. Unconditional physical pulse
bounds, complete residual improvement, infinite realization and blow-up
remain unproved.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`. Manufactured leaves test object identity,
coefficient composition and requests; independent symbolic controls test fast
cancellation and the angular functional. They are not continuum enclosures.
