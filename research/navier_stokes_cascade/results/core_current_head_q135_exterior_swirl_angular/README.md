# Q135: the whole exterior leading swirl angular envelope

Q134 closed real heat endpoint jets and the same three-row compensator.
This checkpoint assembles those estimates with the original exterior
schedule, covering ALL X>=Xstar, rather than only the first reserved patch.
The result is a conditional leading E0/F0 angular bound. The whole-profile
axis and finalized compact F0 input, mixed radial bounds and physical pulse
estimates remain obligations. Historical frozen evidence is unchanged.

## Same reference and same roots

Retain the factorization in [global_swirl_envelope](../global_swirl_envelope/README.md),
section 3: E=S gamma R_H R_A T_o M_H. Here S is positive, independent of
eta, nonincreasing and bounded by its entrance cap S0 (the retained native
choice is S0=2 Ahi). Every stage endpoint is independent of eta. The
factors use their original real radial supports, equal to one outside the
corresponding edits; gamma alone changes from f=(1+eta²)^-1 to 1/2.
The scalar schedule includes the entire axial profile pulse and subsequent
release, hold and waiting stages. Q77 prescribes E=e_b f exp(alpha y)
on the axial stage; its angular strength and end corrections change U,
not that prescribed E. No estimate for a physical high-frequency pulse
is inferred from this slow-profile factorization.

Use 0<lambda<=2^-256, h=lambda², Xstar>=100. All accepted-root and
reference-profile premises of Q126/Q134 are retained. In particular the
later angular root uses the COMPLETE PRE-HEAT future mass and energy.
The heat compensator cancels the complete terminal differences, including
the marginal renormalized-I tail. Its compensation-only partial pressure
and energy are not zero. Pi(X)=Pi(Xstar)+integral E²/(2x)dx retains the
actual edited field and all five incoming traces. This assembly neither
reselects either root nor discards any of that pressure memory.

## Five real-angular factor caps

All raw caps below have form P B^m(m!)², at every finite order m>=0,
uniformly over their entire real radial stage. They also hold for the
identity extensions of relative factors outside their original supports.

| Factor | P | B | Derivation |
|---|---:|---:|---|
| gamma | 81 | 2 | Radius-1/2 disks for f and log(2f) only |
| R_H | 1+15360000 g | 512 | Q134 same heat root, g=h^4/Xstar, three bumps of value <=64 |
| R_A | 1+128 r | 2 | Q126 same angular root disk cap r=2^19 lambda^30, two bumps of value <=64 |
| T_o | (1-h/1000)^-1 | 0 | Original eta-independent terminal taper |
| M_H | 1+4h/Xstar | 32 | Q134 real Gamma bound and original 0<=chi<=1 |

For gamma, before interpolation |f|<=2 on each radius-1/2 disk.
During interpolation gamma=exp[-log 2+theta log(2f)] with theta in [0,1]
and |log(2f)|<=4, hence |gamma|<=exp(4)/2<81. Afterwards gamma=1/2.
Cauchy gives 81*2^m*m!, which is bounded by the stated Gevrey-2 cap.
This uses no complex extension of a radial cutoff or of Gamma.

For R_H, Q134 bounds the value of each accepted coefficient by 6400g
and positive derivatives by 80000g*512^m*(m!)². Their common padded
prefactor is 80000g. Three original unnormalized width-.3 sigma-prime
bumps give 3*64*80000g=15360000g; adding one handles order zero.
For R_A, Q126's SAME complex angular root has modulus <=r on the
radius-1/2 disks. Its two original bumps give the factor 128r; Cauchy
again supplies a Gevrey-1 bound that can be padded to Gevrey-2.

The terminal taper is exactly
[1-rho sigma((3-y)/2)]/(1-rho), rho=h/1000. Its value lies between one
and (1-rho)^-1, and its positive eta derivatives vanish. The real Gamma
multiplier is 1+chi(H-1). Where it differs from one, X>=XK>=Xstar,
so Q134 gives the last row including one-sided pole limits. No complex
Gamma disk, negative argument or analytic slow-angle heat root is used.

## Product and factorial norm

Leibniz gives the real Gevrey-2 product rule
(P,A)*(Q,B) -> (PQ,A+B): the normalized factorial weight is
1/binom(m,j)<=1 and the remaining sum is bounded by (A+B)^m.
Thus the four angular bases sum to 2+512+2+32=548. With

    P = 81 S0 (1+15360000g)(1+128r)(1+4h/Xstar)/(1-h/1000),

we obtain, uniformly over the WHOLE exterior,

    |d_eta^m E0| <= P 548^m(m!)²,
    |d_eta^m F0| <= P 548^m(m!)²,  X>=Xstar.             (1)

The second inequality uses F0=E0/sqrt(2X) and sqrt(2X)>=1; this division
is used only away from the axis. It supplies no radial derivative gain.
The possibly overlapping supports were bounded by a product, so overlap
does not require a new disjointness assumption.

For Q125's intended factorial angular algebra norm the bound is

    sup_(X>=Xstar) sum_(j=0)^m sup_eta |d_eta^j F0|/j!
       <= P sum_(j=0)^m 548^j j!.                       (2)

The same cap controls the stronger sum of global spatial suprema, since
(1) is uniform order by order. A binary budget is
ceilLog(max(1,P))+ceilLog(m+1)+10m+m ceilLog(m), with the last term zero
at m=0. This is a finite bound at arbitrary order, not an analyticity claim.

## Verification and boundary of the result

Exact checks cover factor constants, scalar/analytic padding, products,
factorial normalization, large-order binary budgets, domain/identity guards
and an independent byte-identical evidence replay through Q134. These
formula checks are not proofs of the retained continuum premises.

The next inequality remains the WHOLE-profile F0 norm through m on the
same cover as Q133's finalized U0 norm through m+2. Supply the actual
regular-axis and finalized compact F0 bounds, then join their pointwise-X
norms with (2) by a regional maximum. Do not divide by sqrt(2X) at X=0
or use pre-modulation bounds for the final compact profile. Q121 further
requires mixed X,eta F0 derivatives through m and U0 through m+1 on the
whole physical hull; pure angular estimates cannot replace them.
Actual supports, low-order gates and retained scales still need validation.
Uniform physical pulse bounds, full covariance means (including axial
velocity and pressure), complete residual improvement and blow-up remain
unproved. Run test_audit.py and audit.py --verify-record evidence.json
--output a fresh path.
