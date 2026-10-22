# The retained frame estimates close the tangent propagator input

Parent: `935dae292935babd245a1839126ca246339cb04a`.
The weighted energy argument already exists in `sourced_inverse_foundation`.
This checkpoint substitutes the ORIGINAL explicit CE,Cd and normal-length
budgets, proves the conversion to the physical Euclidean tangent norm,
and feeds the resulting constant and original length bound to the first
sourced inverse. Thus its tangent-propagator and length upper bounds are
no longer independent inputs, relative to the retained valid primitive
background/frame bounds. No new accepted-band inequality is needed.

The canonical-frame entrance mismatch is NOT an extra obstruction for the
zero-entry forced inverse. It still prevents silently identifying the two
homogeneous primary pulses or their roots. No such identification is made.
Actual primitive validity, full source/mean/covariance closure and the
complete physical residual remain open. Historical evidence is unchanged.

## 1. Same normal, same reference envelope, whole closed interval

Use the full original normal n and An, including full F,G and H_Z and
the moving-normal pressure term, on 0<=v<=L. The accepted-tail result in
`core_normal_length_closure` gives the original normal thresholds and
1<=L<=S=ell^2. Retain the ORIGINAL primitive bounds, not newly enlarged
ones silently inserted into the old length selector. In the inherited
actual-normal frame Ba=U Ma from `frame_column_foundation`,

    Ba_left(An Ba-Ba') = diag(lambda,-lambda)+E,
    ||E||_2 <= CE/L,
    d=epsilon*k^2*|n|^2,  |d-dref|<=Cd/L,
    lambda=lambda0/sqrt(1+sref^2)>0,
    dref=lambda0*(1+sref^2)/(1+u^2)^(3/2).

These are equations (8),(9) of that checkpoint, with its full connection
and operator-norm constants. The SAME reference envelope satisfies

    P'/P=lambda-dref,   P(L/2)=1,   P>0.

The normalized amplitude is z/P, not z divided by a numerical lower
bound for P. No inverse Gaussian bound is needed, including near endpoints.
Continuous slow and transverse parameters stay fixed along each common
lifted path; its v coordinate ranges over the same [0,L]. Discrete labels,
carriers, c0, epsilon and k are frozen for this path. Derivatives and
cutoff/support properties retain their existing separate recipes.

## 2. Two-time energy estimate, all nonzero harmonics

For the homogeneous constrained equation t'=An t-m^2 d t, t=Ba z,
put y=z/P and delta=m^2-1>=0. Then

    y' = [diag(0,-2lambda)+E-(d-dref)I-delta*d I] y.

For real or complex y the first diagonal matrix is negative semidefinite.
Taking the real Hermitian pairing gives

    (1/2)(|y|^2)' <= [(CE+Cd)/L-delta*d] |y|^2.

Set H=CE+Cd. On EVERY subinterval 0<=w<=v<=L, integration gives the
coordinate propagator estimate

    ||V_frame,m(v,w)||_2
       <= exp(H*(v-w)/L) * P(v)/P(w)
          * exp[-delta*integral_w^v d]
       <= exp(H) * P(v)/P(w) * exp[-delta*d_* (v-w)].   (1)

Here d_* = bmin^2/4, because epsilon*k^2>=1 and the original normal
floor is |n|>=bmin/2. There is no m^2 Cd loss: split off the EXTRA
ACTUAL nonnegative damping, instead of comparing all m^2 d with dref.
Equivalently V_m=exp[-delta integral d] V_1 by scalar commutation and
uniqueness; E and the diagonal part need not commute with each other.
This introduces neither a higher-harmonic threshold nor a new first band.

The finite constant exp(H) is allowed in a polynomial-in-S estimate:
H depends on the fixed primitive bounds, not on S or label. This is not
the unusable exp(C L) obtained by ignoring the dissipative reference
part. If the supplied primitive bounds vary with S, uniformity must first
be established; this checkpoint does not hide that dependence.

## 3. Explicit physical Euclidean conversion

Write a=n_tan/|n_tan|, b=(a_z,-a_theta), sa=n_r/|n_tan|,
wa=sqrt(1+sa^2), and c=c0. In the inherited frame

    Ba = [er-sa*a,b] [[1,1],[c*wa,-c*wa]].

The two columns of U=[er-sa*a,b] are orthogonal with squared norms
wa^2 and 1. Consequently

    Ba^T Ba = wa^2 [[1+c^2,1-c^2],[1-c^2,1+c^2]],

whose eigenvalues are 2wa^2 and 2c^2 wa^2. On the actual tangent plane,
the specified Ba_left is the inverse coordinate map, so for tangent t

    |Ba_left(w)t| <= |t|/[sqrt(2)*wa(w)*min(1,|c|)].

This does NOT replace the specified ambient left inverse by a
Moore--Penrose inverse or assert equality of their ambient operator norms.
Together with |sa|<=S0 from the original normal/frame thresholds,
cmin<=|c|<=cmax and wa(w)>=1, the endpoint conversion costs at most

    Kframe = (1+S0)*max(1,cmax)/min(1,cmin).

Thus the physical tangent propagator satisfies

    ||V_T,m(v,w)||_2 <= T0 * P(v)/P(w)
                       * exp[-delta*d_* (v-w)],
    T0=Kframe*exp(H).                                  (2)

The bound is for the restricted map n(w)^perp -> n(v)^perp, including
arbitrary tangent initial vectors, not merely the growing primary column.
Omitting the two endpoint frame factors would not prove the requested
physical Euclidean bound.

For an exact positive rational polynomial budget put N=ceil(H)>=0 and

    T(S)=Kframe*3^N,    Length(S)=S.                     (3)

The series for e gives e<3 (n!>=2^(n-1), strictly for some n>=3), hence
exp(H)<=3^N. Equations (2),(3) close the T and length slots in the
unchanged Gaussian-adapted `ForcedBudget`. `forced_from_frame` derives
both slots and uses the original 4/bmin^2 denominator via the preceding
bridge. There is no caller-supplied tangent certificate to assume.
The exact finite integer recipe can be retained without materializing
an astronomically large selected constant. The executable materializer
is exercised only on supplied manufactured finite inputs in this audit.

CE, Cd and Kframe themselves depend on the original primitive floors,
including bmin. Thus this composed bound need not be polynomial in
q=4/bmin^2 when those primitives change. The preceding checkpoint's
q-degree statement explicitly holds all OTHER inputs, including T, fixed.
The present conclusion is uniform polynomial dependence on S for the
fixed original primitive bounds; these are different assertions.

For |m|>=2, (2) also gives a tangent zero-order residence bound
Kframe*exp(H)*min(L,1/[(m^2-1)d_*]); |m|=1 retains length L.
This residence gain is NOT substituted for the length in the adapted
ambient derivative inverse: its analytical normal extension evolves at
P'/P, and derivatives still charge all damping coefficient jets.

## 4. Canonical entrance and the frame mismatch

Let Br=Ba T be the manuscript-reference mixing frame. The exact identities
already proved in `frame_column_foundation` give

    G_r = T^-1 G_a T - T^-1 T',  f_r=T^-1 f_a,
    V_r(v,w)=T(v)^-1 V_a(v,w) T(w).

Multiplying by the endpoint frames shows identical physical propagators.
For zero entrance data, T(0)^-1*0=0 exactly; transforming the entire
forcing and the derivative term preserves the SAME physical solution and
pressure. We may therefore ESTIMATE the canonical zero-entry physical
inverse using Ba, without replacing any actual pulse, source, cutoff or
root. No CEref smallness gate or reset of a nonzero primary entrance is
needed for (2). The path entrance is fixed in the lifted coordinates;
coordinate descriptions with moving entrance must retain their endpoint
chain terms before applying the inherited mixed-jet recurrence.

By contrast, a homogeneous primary initialized by P(0)eplus in each
different frame generally gives different physical data. That previously
recorded mismatch, and its covariance/root consequences, remain. The
present upper propagator estimate proves no identity of those primaries.
Likewise the ODE pressure retains n.f even if Pi*f=0. Global Pi1/Pi2
memory, moving moments and fresh axial velocity means are unchanged;
this local propagator estimate is not a replacement for them.

## 5. Closure ledger and validation

Removed independent inputs: physical Gaussian-ratio tangent bound and
whole-path length upper bound for this original accepted-tail inverse.
Removed overstated obstruction: canonical entrance-frame mismatch for
the ZERO-entry forced inverse. Added analytical assumptions: **none**.
Retained: validity of the original full background/frame primitive bounds
and same accepted-tail/domain; actual coefficient/source jet and whole-path
support bounds; full covariance, original roots and fresh mean identities.
An unknown separately fixed historical numerical band is not certified.

Exact rational controls prove the Gram/singular-square formulas, verify
positive-semidefinite energy gaps for ALL vectors in manufactured cases,
and check both endpoints of the frame change and noncommuting scalar
damping. Negative controls detect omission of damping deviation and frame
conversion. The original normal/denominator and full-pressure callback
joins are exercised. These finite controls corroborate the analytical
proof; they do not enclose an actual giant fluid profile.

Next obligation: pass the generated COMPLETE covariance through the fresh
axial velocity and full-memory pressure means and justify the resulting
full-path source/support estimate on one domain. Resolve primary-frame/root
identification wherever that calculation uses it. All physical scales,
viscosity, cutoff and interaction residual terms remain obligations before
unconditional pulse bounds, residual improvement or blow-up.

```sh
python3 research/navier_stokes_cascade/results/core_tangent_propagator_closure/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_tangent_propagator_closure/audit.py \
 --verify-record research/navier_stokes_cascade/results/core_tangent_propagator_closure/evidence.json \
 --output /tmp/core_tangent_propagator_closure_fresh.json
```
