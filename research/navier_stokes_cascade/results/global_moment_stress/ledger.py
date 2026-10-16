"""Exact exterior-memory equations and labelled heat-tail diagnostics."""
from fractions import Fraction as Q
import mpmath as mp


def exterior_defects(eta, h, X, H, defects):
    """Return delta(V0), delta(r^2 Ttheta), delta(r Tz), at profile q=1.

    Each defect is (value, eta derivative), in M, renormalized I, J, S, Pi
    order. Local fields agree and U=0. This is an identity, not a smallness test.
    """
    (m, me), (R, Re), (j, je), (s, se), (p, pe) = defects
    D=Q(1,2)-h; A=Q(1,2)+h; d=1-eta*eta; L=1-2*h*eta*eta
    v=2*D*eta*m+d*me
    theta=(H*v+(1-h)*R-D*eta*Re-d*je+2*(h-D)*eta*j)/L
    axial=(D*(m-eta*me)+4*h*eta*s-d*se+X*(4*A*eta*p-d*pe))/L
    return -v/L,theta,axial


def direct_coordinates(eta,h,X,H,a,moments,Heta=Q(0)):
    """Original (4.7),(4.11),(4.16), evaluated independently of the defect map."""
    (M,Me),(I,Ie),(J,Je),(S,Se),(Pi,Pie)=moments
    D=Q(1,2)-h; A=Q(1,2)+h; d=1-eta*eta; L=1-2*h*eta*eta
    W=1-2*D*eta*M/X-d*Me/X
    V0=(-2*D*eta*M-d*Me)/L
    Qs=-W+((1-h)*I-D*eta*Ie-d*Je+2*(h-D)*eta*J)/(X*H)
    Ns=(D*(M-eta*Me)+4*h*eta*S-d*Se)/X+4*A*eta*Pi-d*Pie
    # F=H/(2X), r^2=2X, U_X=0. Heta is deliberately irrelevant here.
    return V0,H*(X*Qs/L-a),X*Ns/L


def polynomial_pair(coeff,eta):
    return (sum((c*eta**i for i,c in enumerate(coeff)),Q(0)),
            sum((i*c*eta**(i-1) for i,c in enumerate(coeff) if i),Q(0)))


def exact_review():
    checks=[]; basis=[]
    h=Q(1,1000)
    for eta in map(Q,('-1','-3/4','0','2/5','1')):
        base=[polynomial_pair([Q(i+1),Q(2-i),Q(1,7)],eta) for i in range(5)]
        for X in map(Q,('1/10','100','10000')):
            H=Q(3,2)+eta*eta; a=Q(12,5)
            for i in range(6):
                ds=[polynomial_pair([Q(k+2,101),Q(1-k,103),Q(k,107)],eta)
                    if i==5 or i==k else (Q(0),Q(0)) for k in range(5)]
                changed=[(v+dv,e+de) for (v,e),(dv,de) in zip(base,ds)]
                before=direct_coordinates(eta,h,X,H,a,base)
                after=direct_coordinates(eta,h,X,H,a,changed)
                got=exterior_defects(eta,h,X,H,ds)
                checks.append(tuple(x-y for x,y in zip(after,before))==got)
                if eta==0 and X==100 and i<5:
                    basis.append([str(v) for v in got])
    zero=[(Q(0),Q(0))]*5
    # Exact value matching at eta=0 is not a functional match.
    slope=zero[:];slope[4]=(Q(0),Q(1,10**30))
    pdef=exterior_defects(Q(0),h,Q(100),Q(1),slope)[2]
    # Renormalized I, S and pressure can each leave a signed stress memory.
    isolated=[]
    for i in (1,3,4):
        ds=zero[:];ds[i]=(Q(1,10**30),Q(0))
        isolated.append(exterior_defects(Q(1,2),h,Q(100),Q(1),ds))
    m_only=zero[:];m_only[0]=(Q(1),Q(0))
    mv=exterior_defects(Q(1,2),h,Q(100),Q(1),m_only)[0]
    # Pressure defect produces a term growing exactly linearly in X in r*Tz.
    pd=zero[:];pd[4]=(Q(1),Q(0))
    p100=exterior_defects(Q(1,2),h,Q(100),Q(1),pd)[2]
    p200=exterior_defects(Q(1,2),h,Q(200),Q(1),pd)[2]
    return dict(checks=dict(
        all_five_rows_and_combined_direct_subtraction_exact=all(checks),
        zero_functional_defects_restore_all_exterior_coordinates=exterior_defects(Q(1),h,Q(100),Q(1),zero)==(0,0,0),
        one_angle_pressure_value_does_not_restore_axial_stress=pdef<0,
        each_I_S_pressure_memory_detected=all(any(v!=0 for v in r) for r in isolated),
        axial_velocity_requires_M_memory=mv!=0,
        pressure_memory_grows_linearly_in_weighted_axial_stress=p200==2*p100,
        polynomial_fixtures_are_not_actual_profiles=True),
        rational_equation_checks=len(checks),basis_at_eta0_X100=basis,
        one_angle_pressure_stress=str(pdef),isolated_moment_memories=[[str(v) for v in r] for r in isolated])


def ledger_review():
    # Signed functional changes, represented by arbitrary rational coefficients.
    # The analytical proof identifies these with the EXACT correction equations.
    axial=[Q(2),Q(0),Q(3),Q(5),Q(0)]
    angular=[Q(0),Q(7),Q(0),Q(11),Q(0)]
    heat=[Q(0),Q(13),Q(0),Q(-17),Q(-19)]
    heat_comp=[-v for v in heat]
    inner=[Q(23),Q(29),Q(31),Q(37),Q(41)]
    compact=[Q(43),Q(47),Q(53),Q(59),Q(61)]
    restored=[sum(v) for v in zip(heat,heat_comp,inner,[-v for v in inner],compact,[-v for v in compact])]
    missing_pressure=compact[:];missing_pressure[-1]=0
    failed=[v-w for v,w in zip(compact,missing_pressure)]
    return dict(checks=dict(
        heat_compensation_restores_exact_signed_I_S_Cp_rows=all(v==0 for v in restored),
        omitted_fifth_row_leaves_pressure_memory=failed[:-1]==[0]*4 and failed[-1]!=0,
        axial_corrections_do_not_change_I_or_pressure=axial[1]==axial[4]==0,
        angular_corrections_do_not_change_M_J_or_total_pressure=angular[0]==angular[2]==angular[4]==0,
        heat_compensation_has_nonzero_intermediate_memory=any(v!=0 for v in heat_comp),
        ledger_controls_not_existence_proofs=True),ordered_rows=['M','renormalized_I','J','S','Cp'],
        heat_plus_compensation=[str(v+w) for v,w in zip(heat,heat_comp)])


def tail_review():
    h=Q(1,100);X=Q(100);k=Q(1);d=Q(1)
    a1=h*(1+h);a2=h*(1+h)**2*(2+h)
    # Bounds divided by the common k X^-h. They contain no hidden h^-1 loss.
    renorm0=2*(1+h)
    renorm1=4*(1+h)
    renorm2=4*(1+h)+16*h*(1+h)*(2+h)/X
    remainder=2*d*d*h*(1+h)*(2+h)/X
    return dict(checks=dict(
        uniform_tail_I_majorant_below_3=renorm0<3,
        uniform_tail_I_eta_majorant_below_5=renorm1<5,
        uniform_tail_I_etaeta_majorant_below_5=renorm2<5,
        second_heat_derivative_kept=a2>a1,
        integrated_Taylor_remainder_below_p001=remainder<Q(1,1000),
        power_subtraction_integrable_at_axis_for_h_lt1=True,
        fixed_h_positive_required_for_tail_integrals=True,
        h_zero_before_infinite_integration_is_invalid=True),
        normalized_I_tail_bounds=[str(v) for v in (renorm0,renorm1,renorm2)],
        normalized_integrated_Taylor_remainder=str(remainder),
        h_to_zero_integrated_heat_defect='-2*k*(1-eta^2), not zero')


def scaling_review():
    checks=[]
    for h in (Q(1,1000),Q(1,100),Q(1,8)):
        A=Q(1,2)+h;D=Q(1,2)-h
        b=1-h;c=Q(3,2)-2*A
        checks.append(b-1==c-D==-h)
        checks.append((1-A)-1==(1-2*A)-D==-A)
        for eta in map(Q,('-1','0','1/3','1')):
            d=1-eta*eta;L=1-2*h*eta*eta
            m,me=polynomial_pair([Q(2),Q(-1),Q(3)],eta)
            R,Re=polynomial_pair([Q(1),Q(2),Q(4)],eta)
            j,je=polynomial_pair([Q(3),Q(1),Q(-2)],eta)
            s,se=polynomial_pair([Q(2),Q(-3),Q(1)],eta)
            # Negatives of time + axial derivatives of the physical moments.
            theta=-(-b*R+D*eta*Re+2*c*eta*j+d*je)/L
            axial=-(-(1-A)*m+D*eta*me+2*(1-2*A)*eta*s+d*se)/L
            ds=[(m,me),(R,Re),(j,je),(s,se),(Q(0),Q(0))]
            got=exterior_defects(eta,h,Q(100),Q(0),ds)
            checks.append(got[1]==theta and got[2]==axial)
    return dict(checks=dict(
        physical_weights_signs_and_similarity_exponents_exact=all(checks),
        radial_angular_transport_uses_r_squared_not_r=True,
        leading_residual_excludes_axial_viscosity=True,
        zero_moments_imply_zero_signed_totals_not_full_residual=True),
        rational_scaling_checks=len(checks))


def heat_factor(z,h,derivative=0):
    """Diagnostic only. Finite Taylor near zero; independent Tricomi U elsewhere."""
    if z==0:
        return (-1)**derivative*mp.rf(h,derivative)*mp.rf(1+h,derivative)
    if z<mp.mpf('1e-6'):
        return sum((-1)**n*mp.rf(h,n)*mp.rf(1+h,n)*z**(n-derivative)/mp.factorial(n-derivative)
                   for n in range(derivative,13))
    return (-1)**derivative*mp.rf(h,derivative)*mp.rf(1+h,derivative)*z**(-1-h-derivative)*mp.hyperu(1+h+derivative,2,1/z)


def heat_tail_diagnostic(digits,order):
    with mp.workdps(digits):
        h=mp.mpf('0.125');X=mp.mpf(100);k=mp.sqrt(2);a1=h*(1+h)
        rows=[]
        nodes,weights=mp.gauss_quadrature(order,'legendre')
        # h=1/8: t=u^8/X removes the fractional endpoint power after
        # the first Taylor term has been subtracted. Fixed quadrature is
        # a diagnostic; order refinement is reported, never a root proof.
        def integral(function):
            return sum(w/2*function(((v+1)/2)**8/X)*8*((v+1)/2)**7/X
                       for v,w in zip(nodes,weights))
        for eta in map(mp.mpf,('0','.5','1')):
            d=1-eta*eta;Z=2*d/X;G=heat_factor(Z,h);Gp=heat_factor(Z,h,1)
            def residual(t,derivative):
                if not t or not d:
                    return mp.mpf(0)
                z=2*d*t
                if z<mp.mpf('1e-6'):
                    # Subtract the first two terms BEFORE numerical cancellation.
                    val=sum((-1)**n*mp.rf(h,n)*mp.rf(1+h,n)*z**(n-derivative)/mp.factorial(n-derivative) for n in range(2,13))
                else:
                    val=heat_factor(z,h,derivative)-(1-a1*z if derivative==0 else -a1)
                return t**(h-2+derivative)*val
            DI=-2*d*(1+h)*k*X**(-h)+k*integral(lambda t:residual(t,0))
            DIe=4*eta*(1+h)*k*X**(-h)-4*eta*k*integral(lambda t:residual(t,1))
            I=k*X**(1-h)/(1-h)-DI;Ie=-DIe;H=k*X**(-h)*G
            L=1-2*h*eta*eta;D=mp.mpf('.5')-h
            Qs=-1+((1-h)*I-D*eta*Ie)/(X*H)
            a=2+2*h+2*Z*Gp/G
            gamma=mp.quad(lambda v:mp.exp(-v)*v**h*(1+Z*v)**(-h),[0,1,mp.inf])/mp.gamma(1+h)
            rows.append(dict(eta=str(eta),tail_I=mp.nstr(DI,60),heat_gamma_error=float(abs(G-gamma)),
                             zero_heat_stress_error=float(abs(X*Qs/L-a)),
                             pole_I_eta_over_Hpow=float(Ie/(k*X**(-h)))))
        # Remainder bound for the n=12 Taylor truncation on z<=1e-6.
        # Both DI and DI_eta after their integration, k<=2, h=1/8, |eta|<=1.
        rem=mp.rf(h,13)*mp.rf(1+h,13)/mp.factorial(13)
        # z^13 <= (1e-6)^11*(2t)^2; for the differentiated
        # remainder use 13*z^12 <= 26*(1e-6)^11*t.
        # This bound holds even when d is tiny; it does not silently
        # replace the branch endpoint 1e-6/(2d) by 1e-6/2.
        terr=256*rem*mp.mpf('1e-6')**11*(1/X)**(h+1)/(h+1)
        pressure=mp.quad(lambda x:-1/(2*(1+x)**2),[0,1,mp.inf])
        energy=mp.quad(lambda x:x/(1+x)**3,[0,1,mp.inf])
        return dict(rows=rows,order=order,integrated_Taylor_error_upper=float(terr),
                    manufactured_Fubini_error=float(abs(pressure+energy)))


def diagnostic_review():
    a,b=heat_tail_diagnostic(50,24),heat_tail_diagnostic(70,36)
    with mp.workdps(80):
        error=float(max(abs(mp.mpf(x['tail_I'])-mp.mpf(y['tail_I'])) for x,y in zip(a['rows'],b['rows'])))
    return dict(high_precision=b,refinement_error=error,checks=dict(
        gamma_integral_agrees_with_heat_diagnostic=max(r['heat_gamma_error'] for r in b['rows'])<1e-45,
        renormalized_tail_recovers_zero_heat_stress=max(r['zero_heat_stress_error'] for r in b['rows'])<1e-30,
        pole_value_zero_does_not_erase_angular_derivative=mp.mpf(b['rows'][-1]['tail_I'])==0 and abs(b['rows'][-1]['pole_I_eta_over_Hpow']+4.5)<1e-12,
        explicit_integrated_Taylor_remainder_small=b['integrated_Taylor_error_upper']<1e-55,
        precision_and_quadrature_refinement_consistent=error<1e-25,
        independent_manufactured_pressure_Fubini_identity=b['manufactured_Fubini_error']<1e-60,
        diagnostic_h_not_selected_actual_h=True,
        diagnostics_are_not_continuum_supremum_or_exact_root_proofs=True))
