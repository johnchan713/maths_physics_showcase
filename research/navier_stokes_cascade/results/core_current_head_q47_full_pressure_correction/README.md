# Q47 correction: the initial tail and current rows require full compact pressure

Parent `9e4097774930b205204170f4e14394197877f3fc` (Q46).
While tracing the original axial A1 remainder, we found that Q40 assigns
Q38.recipe(...).Kp to ALL THREE non-Haar row changes. Q38 proves that
bound for the ordinary two-edge primitive on radially zero-integral
inputs. Its Haar pressure calculation remains valid. Its proof does not
cover the operator needed in these full rows. Q40's executable positive
polynomial and passing tests do not resolve this domain mismatch.
Historical frozen files are preserved. Use this corrected factory when
feeding the current fast field to Q39/Q41.

## Exact obstruction

Write h=gr-rho P, P=integral Haar(gr). Only integral Haar(h)=0 is
guaranteed. For a smooth interior bump phi of positive radial mass,
h(R,y)=phi(R) cos(y1) has zero Haar but its ordinary left primitive
equals a nonzero multiple of cos(y1) past the bump. The edge weight
zeta delta^(-v) tends to zero there at the shell endpoint, for any
fixed finite v. Thus a weighted ordinary left primitive estimate on
this class with a finite constant is false. This is a domain witness,
not the selected physical source, and does not assert that a particular
numerical Q38 constant is smaller than the norm of the actual T0.

The pressure is pm=T0 h with the ORIGINAL translated auxiliary argument.
The right branch of T0 on the outer collar and the full A0 remainder are
essential. In particular D0 pm=h-A0 h. Physical angular averaging of
the wave tensor does not kill dependence on the auxiliary angle.
The earlier core_full_compact_pressure_response already distinguishes
these operators; the error here was in the subsequent Q40 composition.

## Correct row estimate at the same state and pole

For n let N=n+3 and v be Q40's common row pole. Apply Q46.build with
e=0 at (N,v) on the same profile shell and width. Its KT bounds the
full T0 in the normalized R norm, by the fixed-shift raw integral
estimate, both chart transfers and the two plateau branches. For each of
the old-fast difference, the slow difference based on the OLD fast state,
and the complete signed covariance change, keep its original Gr and P
bounds and put

    pm_n <= KT0_(N,v) (Gr + rho_N P),    P <= Kmom_v Gr.

This is an upper bound on the compensated pressure, including its moment;
the plus sign is from estimating the subtraction. The axial row still
includes (n+4)pm_n, and the mean differences retain slow time, viscosity,
and all cross terms. Apply Qcirc with norm at most two only after adding
the original initial A1 term. Q40's existing conservative factor of two
on that term is harmless. Q38's Haar operator remains in Q37/Q39 scalar
calculations. No full pressure is replaced by its latest increment.

The factory instantiates TWO channels from the same original cutoff
shape and witnesses: T0 for full pressure, T1/A1 for the axial fast
inverse. Their multipliers are different even on the same hull. For
Q41 output m the largest row is m+10, so the conservative pressure
request N is m+13. Q46's combined operator interface requests the shape
through m+14. The old m+7 cutoff count covered only the new axial fast map;
it is insufficient for this combined factory. This is a sufficient count,
not a lower bound for a KT-only implementation. The corrected row changes
coefficients, not the conditional epsilon powers at fixed geometry/order.

## What remains for the original A1

The tail in Q40 is A1 Qcirc Ez_prefast of the ORIGINAL primary, before
either initial mean or the signed wave. It is not the later post-signed
A1 Dz Qcirc Phi. The initial-cancellation checkpoint gives the retained
strict budget p=99999, J=n+200001 at kappa=1/100000. With source order
B=1-kappa this gives B+p kappa=2-2kappa, sufficient for Q40. Its needed
inequality is a weighted bound of that order for the ACTUAL shifted tail.
The raw kernel estimate uses N_J(f1), N_n(k1 chi_U), and N_J(Qcirc
Ez_prefast); obtaining the last requires full shifted T0 pressure through
J+1 and original covariance through J+3. For n=m+10 this reaches original
covariance m+200014. Neither low-order Q46 cutoff checks nor Q41's
m+27 profile ledger supplies these higher jets. Pressure-compensator terms
can drop out of the projected zero-mean prefast source only by the exact
projection identity, not by applying the Haar primitive to its other modes.
Next bind this separate original-source budget, with its actual fixed
profiles and high pole, before asserting tail closure. Full residual
improvement, uniform actual pulses and blow-up remain unproved.
