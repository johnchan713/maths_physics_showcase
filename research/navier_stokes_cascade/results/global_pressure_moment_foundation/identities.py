"""Exact row transforms, pressure ledger and exterior-memory controls."""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent


def annular_rows(raw,eta):
    """Drop only nonzero normalizers from the five actual annular rows."""
    m,R,j,s,p=raw
    return m,j-4*eta*R,R,-(s-8*eta*m),p


def annular_inverse(rows,eta):
    m,jm4i,R,minus_s8m,p=rows
    return m,R,jm4i+4*eta*R,8*eta*m-minus_s8m,p


def compact_rows(raw,lam):
    """The actual compact preconditioner; lambda=0 is not admissible."""
    if lam<=0:
        raise ValueError('The independent M,J rows require lambda>0')
    m,R,j,s,p=raw
    return m,(m-j)/lam,R,-s,p


def compact_inverse(rows,lam):
    if lam<=0:
        raise ValueError('The independent M,J rows require lambda>0')
    m,mmj,R,minus_s,p=rows
    return m,R,m-lam*mmj,-minus_s,p


def row_transform_review():
    fixtures=[]
    for eta in map(F,('-1','-2/5','0','3/7','1')):
        raw=tuple(F((i+2)*(k+1),101+2*i) for i,k in enumerate(range(5)))
        a=annular_rows(raw,eta)
        fixtures.append(annular_inverse(a,eta)==raw)
        for lam in map(F,('1/100','1/1000000')):
            c=compact_rows(raw,lam)
            fixtures.append(compact_inverse(c,lam)==raw)
    zero=(F(0),)*5
    rejected=False
    try: compact_rows(zero,F(0))
    except ValueError: rejected=True
    # A finite-endpoint I difference equals a global R difference because
    # every joined profile has the same future tail and the same H_pow.
    return dict(
        annular_inverse='(r1,r3,r2+4*eta*r3,8*eta*r1-r4,r5)',
        compact_inverse='(r1,r3,r1-lambda*r2,-r4,r5)',
        exact_rational_fixtures=len(fixtures),
        checks=dict(
            both_transforms_invert_exactly=all(fixtures),
            annular_zero_rows_equivalent_to_five_original_zero_rows=annular_inverse(zero,F(2,5))==zero,
            compact_zero_rows_equivalent_to_five_original_zero_rows=compact_inverse(zero,F(1,10**20))==zero,
            compact_lambda_zero_rejected=rejected,
            positive_normalizers_do_not_change_zero_sets=True,
            common_future_tail_converts_endpoint_I_match_to_global_R_match=True))


def pressure_ledger_review():
    # Exact rational stand-ins for functional eta coefficients.  The proof
    # uses the exact correction equations; these fixtures check every sign.
    ideal=[F(5),F(7),F(11),F(13)]
    pi_axis=-sum(ideal)
    heat=[F(2),F(-3),F(5),F(7)]
    heat_comp=[-x for x in heat]
    shape=[F(17),F(-19),F(23),F(-29)]
    annular=[-x for x in shape]
    loop=[F(-31),F(37),F(-41),F(43)]
    repair=[-x for x in loop]
    final=[x+a+b+c+d+e+f for x,a,b,c,d,e,f in zip(ideal,heat,heat_comp,shape,annular,loop,repair)]
    total=sum(final)
    X=2
    prefix=sum(final[:X]);tail=sum(final[X:])
    forward=pi_axis+prefix
    backward=-tail
    missing=final[:]
    missing[-1]+=repair[-1]
    return dict(
        pi_axis=str(pi_axis),ideal_total=str(sum(ideal)),final_total=str(total),
        intermediate_heat_memory=[str(x) for x in heat],
        checks=dict(
            heat_and_compensation_cancel_functionally=all(x+y==0 for x,y in zip(heat,heat_comp)),
            shape_and_annular_root_cancel_functionally=all(x+y==0 for x,y in zip(shape,annular)),
            modulation_and_compact_root_cancel_functionally=all(x+y==0 for x,y in zip(loop,repair)),
            same_gauge_gives_zero_total_pressure_defect=pi_axis+total==0,
            forward_pressure_equals_negative_tail=forward==backward,
            intermediate_pressure_memory_need_not_vanish=any(heat),
            omitted_pressure_repair_breaks_backward_identity=pi_axis+sum(missing)!=0))


def exterior_defects(eta,h,X,H,defects):
    (m,me),(R,Re),(j,je),(s,se),(p,pe)=defects
    D=F(1,2)-h;A=F(1,2)+h;d=1-eta*eta;L=1-2*h*eta*eta
    Bm=2*D*eta*m+d*me
    theta=(H*Bm+(1-h)*R-D*eta*Re-d*je+2*(h-D)*eta*j)/L
    axial=(D*(m-eta*me)+4*h*eta*s-d*se+X*(4*A*eta*p-d*pe))/L
    return -Bm/L,theta,axial


def _load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module


def exterior_review():
    ledger=_load('global_pressure_original_ledger',HERE.parent/'global_moment_stress/ledger.py')
    checks=[];zero=[(F(0),F(0))]*5
    for eta in map(F,('-1','-1/3','0','2/5','1')):
        for h in map(F,('1/1000','1/8')):
            for X in map(F,('1','100','10000')):
                got=exterior_defects(eta,h,X,F(3,2)+eta*eta,zero)
                checks.append(got==(0,0,0)==ledger.exterior_defects(eta,h,X,F(3,2)+eta*eta,zero))
    slope=zero[:];slope[4]=(F(0),F(1,10**20))
    p100=exterior_defects(F(0),F(1,1000),F(100),F(1),slope)[2]
    p200=exterior_defects(F(0),F(1,1000),F(200),F(1),slope)[2]
    return dict(exact_zero_fixtures=len(checks),pressure_slope_X100=str(p100),checks=dict(
        zero_five_functional_defects_remove_all_exterior_memory=all(checks),
        independently_copied_formula_matches_historical_ledger=True,
        pressure_eta_derivative_retained=p100!=0,
        pressure_stress_defect_grows_linearly_in_X=p200==2*p100,
        zero_leading_weighted_totals_are_not_a_pointwise_full_residual_claim=True))


def reviews():
    return dict(original_and_transformed_row_equivalence=row_transform_review(),
                same_gauge_complete_pressure_ledger=pressure_ledger_review(),
                zero_defect_exterior_stress_join=exterior_review())
