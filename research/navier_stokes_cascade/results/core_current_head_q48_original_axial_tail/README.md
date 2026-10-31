# Q48: the original axial tail in native coordinates

Parent: 4a9deeb3e0b942b8241f0ff745284ed657c9398d (Q47). This concerns
only the initial projected source, before the initial means and signed
wave. We give a sufficient bound for its shifted axial remainder and
save two source derivatives over the previous Fourier count. We also
correct an omission in applying the historical raw estimate: factorial
conversions within one coordinate system do not convert between R and
U. This changes coefficients even when it does not change epsilon powers.
No historical records, physical parameters or fields are modified.

## Native estimate

Use the original notation of [the fast inverse, equations (2) and (8)](../core_fast_inverse_budget/README.md).
At fixed band parameters let H be the smooth, compactly zero-extended
source of auxiliary Haar mean zero. Take integer p>=1 and output order n.
For v=(1,1-sqrt(2)), omega_k=v.k, put K=2^h, h>=0, and take
K<=|k|_1<2K. The conjugate quadratic norm gives
1/|omega_l| < 3|l|_1 for every nonzero integer l. Distinct divisors
on the shell, including zero as a comparison point, are separated by
1/(12K): differences have norm <4K. On each side of zero, the j-th
nearest divisor therefore has magnitude >=j/(12K). Consequently

    sum_shell |omega_k|^(-2p) <= 4 (12K)^(2p).

Here sum j^(-2p)<=2. This uses the actual two-dimensional direction,
not an assumed uniform divisor gap. For any smooth torus coefficient G,
Parseval on the unit torus, splitting by the larger coordinate |k_i|>=K/2,
gives the square-sum norm of its shell coefficients at most
sqrt(2) (2/K)^(p+1) times max_i ||partial_yi^(p+1) G||_infty
(dropping favorable factors of 2pi). Cauchy--Schwarz with the preceding
divisor sum leaves K^-1. Summing over dyadic K costs two. In particular
8 * 2^(p+1) * 12^p is a sufficient constant independent of the band.

Differentiate the exact whole-line integral at fixed shift a, as in
Q47's pressure convention: y+Ma v is independent of the output U.
After Fourier expansion, integrate each nonzero mode by parts p times
in a. This costs (2 pi M |omega_k|)^(-p) and p native radial
derivatives of its coefficient. Extension by zero removes the endpoint
terms. Apply the shell bound pointwise in U+a to this differentiated
coefficient and integrate over the native support length L. The summable
bound also justifies summing the integrals. Zero Haar is essential:
a nonzero constant mode would have no M gain.

Let J=n+2p+1, B_n=N_n(k1 chi_U), A=N_(n+p)(f1), in the historical
raw native maximum norm. Only the first n+p derivatives can hit the
radial multiplier f1; the additional p+1 torus derivatives act only
on the source. Product rules at output and inside the integral give

    N_n(A1 H) <= K M^(-p) N_J(H),
    K = 8 * 2^(p+1) * 12^p * 2^(2n+p) * L * B_n * A.

The old n+2p+3 estimate is still sufficient in native coordinates.
The new estimate makes no claim of optimality. Use Q44's Umax-Umin
as an upper bound on L; a logarithmic length alone is not a native
length. For B_n use 2^n N_n(k1) N_n(chi_U) with Q46's moving X
derivatives. Pressure and axial plateaus must each describe their
actual cutoff on the original shell.

## From the original covariance to a weighted tail

Let v_in=12(J+3)+18, w_v=zeta delta^(-v), and W the original primary
covariance bound at N=J+3 in normalized R jets. Use the initial positive
hull and the graph and reciprocal-R sums of FreshMeans._data(N,v_in).
In the convention of the [original projected row](../core_initial_cancellation_fast_source/README.md),
with D=J+4, form

    Gr=(Dg+D+2r)W,   Pm=KT0_(J+1,v_in) Gr,
    Ez=(Dg+r)W+D(W+Pm).

The physical power epsilon^B, B=1-kappa, is factored out as in that row.
KT0 is the full shifted two-branch pressure at the source pole, with
both Q43 chart transfers. Sigma and rho P drop out here only because
they are auxiliary independent and Qcirc commutes with T0. This
argument does not apply to the later mean differences. The projected
original source H has weighted coefficient 2 Ez.

Use F_J from Q43 to take normalized R input to raw U at each input
point and G_n to take raw U output to normalized R at each output
point. For example at U=2, R=U^2, the germ R-4 has normalized
through-one R norm 1 there and raw U derivative 4: a same-chart J!
alone cannot justify that step. Transfer the inequality, never the
differentiated weight. Q45's E(v_in) bounds w_v globally. At the
output every derivative of the axial cutoff factor is supported on
its transition. There 1/zeta <= Z=exp((width^2+4)/gap_ax^2).
Outside that transition the output jets vanish. Thus the high input
pole need not be applied to the rest of Q47's rows: the resulting
output has pole zero and embeds at every nonnegative pole.

For c_i^-1 <= I0 S and the original 0<rho_g<1 graph floor, a sufficient
polynomial for the weighted output is

    P_n(S) = 2 K F_J E(v_in) Ez G_n Z (1+I0 S)^p,
    ||A1 Qcirc Ez_prefast||_(normalized R,n,w_0)
        <= epsilon^(B+p kappa) P_n(S).

Use I0 from the original initial compact archive at order zero; there
is no second temporal inversion. All products refer to the same initial
field, physical chart, annulus and width. The two witness intervals
need not agree. The endpoints come from the initial budget's own geometry: its stored
rational gaps h=2(t-1)/(t+1) invert to t=(2+h)/(2-h), where t is
l^2/(2Xa) or 2Xb/u^2 respectively. Thus we need not identify its
shell with a later pressure-profile hull. The supplied witnesses and
shape jets must still apply to the actual cutoffs. The
program returns coefficient and floor power separately. Q47 applies
another factor two to its supplied tail; if given this projected bound
that is conservative, not another projection needed in this estimate.

## Exact remaining domain and evidence

With kappa=1/100000 and p=99999 the conditional order is 2-2 kappa.
For row zero J=199999, N=200002, v_in=2400042. For postslow output
m the last row n=m+10 requests covariance through m+200012,
KT0 through m+200010 and its fixed shape through m+200011. The
axial output cutoff needs only n+1. The hard remaining input is the
same-state inequality for each complete original covariance entry:
for every point on the whole source hull of every accepted band,

    sum_(|alpha|<=N) |partial_R,Z,T,y^alpha Cov_ij| / alpha!
        <= epsilon W_N(S) zeta delta^(-v_in).

Here Cov_ij denotes the complete original angular covariance. Coverage
through a maximum pole means availability at the requested v_in; a
bound only at a larger pole would not imply this inequality.
This requires the original canonical root, target, wave and complete
curl seminorms through their corresponding nested orders, and the
fixed pressure cutoff seminorms on that shell. A low-order ledger
or a larger declared pole does not prove those premises.

The program computes from the original primary and uses separate
coverage declarations; it does not replace Q47's tail callback with
an unverified number. Only small-order manufactured coefficients were
evaluated. The exact ledger, divisor-spacing diagnostics, shifted-mode
integration and comparison against the independently constructed Q47 KT0
check the implementation. They do not enclose the physical high jets.
Next bind those high-order inputs or find a sharper same-source estimate.
Uniform actual pulses, the complete viscous/cutoff/interaction residual
and blow-up remain unproved.

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
