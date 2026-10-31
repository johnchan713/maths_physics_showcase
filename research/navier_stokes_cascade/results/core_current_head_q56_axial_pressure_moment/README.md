# Q56: the post-signed axial obstruction is a primary covariance moment

On Q54's canonical first-real prefix, the tail T=A₁ Dz QΦ can be
tested without high derivatives of the reconstructed pressure, signed
covariance, or generated axial potential. Under the inherited fixed-order
field, base, operator and whole-support bounds we reduce T, modulo
order 3−2κ, to one explicit shifted moment of the **original complete**
primary tensor. We do not estimate that last moment at the target rate.

## Integrate the actual pressure first

Freeze one common band, U=R^d. Write brackets for integration in dR
on a fixed positive hull [l,h] at the phase y+M(R^d−Uout)v_r. They
act on auxiliary fields, leaving the output slow parameters fixed.
Thus [D_r F]=0 for compact F; along the path D_r is the full derivative.
Take both ends outside all sources and transitions, with zero extensions.
Q is identity minus Haar. The exact [compact maps](../core_fast_inverse_budget/README.md)
give f₀ dU=dR, f₁ dU=R dR, and the output multiplier B₁=k₁χ₁,U.
Use χ₀ for the pressure cutoff expressed at R and set

    cχ = h²/2 − ∫_l^h R χ₀(R) dR = ∫_l^h (R²/2) ∂Rχ₀ dR,
    w(R)=cχ−R²/2.

Indeed T₀g along a single path is its running integral minus χ₀
times its complete integral. Fubini gives [R T₀g]=[w g]:
the first integral weights g(s) by (h²−s²)/2. This uses the same
shift in both maps, including through gaps. Extending the hull in
either plateau leaves cχ unchanged. Projection removes ρP only in
this channel, since it is independent of y. Hence exactly

    T = B₁ Dz Q [R(−2G_R Ψ + γ² + Wzz) + w gr].                  (1)

Dz=ε∂Z at fixed native coordinates. It acts on w, the base and all
fields INSIDE the moment, but not on B₁. No cutoff is divided out.
In particular ∂Z cχ=−∫R ∂Zχ₀ dR on the frozen hull; suppressing
this term would lose pressure memory. This is the full same-state
pressure pm=T₀(gr−ρP), equivalent to the local-plus-primitive form
of the [postwave axial identity](../core_axial_potential_pressure_gain/README.md).

Use the [radial source decomposition](../core_angular_pressure_normal_form/README.md)
at that state. Write S=2bβ+β²+Wrr, H=bγ+Gβ+βγ+Wrz and
J=2Vv+v²+Wθθ. Here v is angular mean velocity. Then gr equals
−tstar β−D₁S−Dz H+J/R+ε(D_r D₁β+Dz²β).
Integration by parts gives [−w D₁S]=[(w'−w/R)S]. The radial vector
viscosity disappears from this moment **exactly**:
[w D_r D₁β]=[(w''−w'/R)β]=0. Both integrations require the actual
smooth zero extension of the velocity, including its axial correction.
Neither viscosity in the axial row nor ε Dz²β in gr has been discarded
by this identity.

## Retain the linear mean response, then substitute its original source

Use O_a for the established fixed-jet bounds at ε^a up to the allowed
polynomial and weight/pole factors. The following error estimates are
inside the brackets, before the exterior Dz. The whole-line integral
and multiplication by B₁ cost no negative power: differentiate at
fixed shift in U as in Q48, using the positive hull, chart conversions
and transition weight floor. No derivative of a moving phase is taken
at fixed absolute integration coordinate for this norm estimate.

Put B=1−κ. The old tangential means have order B and β order B+1;
tstar on β costs at most polynomial factors. Thus γ² and v² are
O_(2−2κ), bβ and β² are smaller, Dz H is O₂ or better, and ε Dz²β
is smaller. Q51 gives the full signed covariance change O₂. After
projection the only sub-order-(2−2κ) terms still needed in (1) are
the original W entries and the linear terms in QΨ and Qv. Slow means
drop out of these two linear terms only; their correlations were
included in the square bounds. The signed wave did not change the
means at this pre-new-fast state.

Let I=c_i^-1 Nt^-1 on zero-Haar fields. The original fast update is
Qv=−I fθ and QΨ=−I T₁ fz; its axial velocity is D₁QΨ, including
the A₁ correction. By the [original projected sources](../core_initial_scalar_seed_budget/README.md)
and [full initial pressure bound](../core_current_head_q49_initial_flux_tail/README.md),

    Qv = −I D₂ QW_rθ + O₂,
    QΨ = −I QW_rz + O_(2−κ).

For the second substitution T₁ D₁ equals identity on the compact
flux: both directed integrals recover R times the flux with the
appropriate sign. T₁ Dz Q(Wzz+pm₀) stays in the error at 2−κ;
no derivative was commuted through T₁. I commutes with radial
integration and translations on the frozen torus. With a=2Vw/R,
[−a D₂ F]=[(a'−2a/R)F]. Substitution in (1) now gives

    T = B₁ Dz [ Q{R Wzz + (w'−w/R) Wrr + (w/R) Wθθ}
                 + 2R G_R I QWrz + (a'−2a/R) I QW_rθ ]
        + O_(3−2κ).                                               (2)

Since w′=−R, the diagonal part is R(Wzz−Wrr)+(w/R)(Wθθ−Wrr).
A common isotropic addition cancels exactly, even before integration.
Every W on this line is original, in common normalized units including
the curl square. Prime differentiates all slow radial factors at fixed
Z. Dz still acts on the entire integrand. The result uses finite extra
jets for I and for the substitutions at the requested output order.
The error bound uses absolute norms only on the higher-order terms;
no smallness of the leading moment is assumed to derive (2).

At 3−5κ the regular axial row meets its required bound if and only if
the displayed moment **with B₁ Dz** does, by (2) and [Q54](../core_current_head_q54_axial_row/README.md).
This isolates what must be checked on the original generated amplitudes.
Q55's resonant example concerns generic inputs; it neither proves nor
disproves cancellation between these five physical entries. Next bound
their combined resonant integral, or supply their high weighted jets
and the corresponding base/profile factors. There is no unconditional
all-scale pulse estimate, radial A₀ estimate, or full residual gain.
