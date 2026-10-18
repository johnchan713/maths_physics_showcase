"""Protocol regression model: stamps are not actual PDE fields or a proof."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Input:
    bundle:str
    stamp:tuple
    kind:str


class State:
    def __init__(self,bundle='canonical-reference-roots'):
        self.bundle=bundle;self.wave=0;self.mean=2
        self.pressure=(0,2);self.targets=(0,2)
        self.tails=[('parent-existing-ledger',0)];self.log=[]
    def stamp(self):return (self.wave,self.mean)
    def clean(self):return self.pressure==self.targets==self.stamp()
    def read(self,kind):
        if not self.clean():raise ValueError('Recompute pressure and targets from current fields')
        if kind not in ('regular','signed','fast','slow','tail'):raise ValueError('Unknown source kind')
        return Input(self.bundle,self.stamp(),kind)
    def recompute(self):
        self.pressure=self.targets=self.stamp();self.log.append(('recompute',self.stamp()))
    def apply(self,inp,operation):
        if inp.bundle!=self.bundle:raise ValueError('Do not mix root/pulse families')
        if not self.clean() or inp.stamp!=self.stamp():raise ValueError('Stale state input')
        allowed={'sourced':'regular','signed':'signed','fast':'fast','slow':'slow'}
        if operation not in allowed or inp.kind!=allowed[operation]:raise ValueError('Wrong source or additive-tail inversion')
        if operation in ('sourced','signed'):
            self.wave+=1;self.tails.append((operation,self.wave))
        else:self.mean+=1
        self.log.append((operation,self.stamp()))


def build(cycles=1):
    if type(cycles) is not int or cycles<0:raise ValueError('Finite nonnegative cycles required')
    s=State()
    # Fresh canonical real signed/mean prefix, first source, generic signed/means.
    for op in ('signed','fast','slow','sourced','signed','fast','slow'):
        s.apply(s.read({'sourced':'regular','signed':'signed','fast':'fast','slow':'slow'}[op]),op);s.recompute()
    for _ in range(cycles):
        for op in ('sourced','signed','fast','slow'):
            s.apply(s.read({'sourced':'regular','signed':'signed','fast':'fast','slow':'slow'}[op]),op);s.recompute()
    return s


def review():
    rows=[];good=True;rejected=[]
    for j in (0,1,4):
        s=build(j);good &= s.clean() and s.bundle=='canonical-reference-roots' and len(s.tails)==4+2*j and s.tails[0]==('parent-existing-ledger',0)
        rows.append(dict(cycles=j,final_wave=s.wave,final_mean=s.mean,parent_ledger_blocks=1,new_tail_entries=len(s.tails)-1,operations=len(s.log)))
    s=State();old=s.read('slow');s.apply(s.read('fast'),'fast')
    for name,action in (
        ('missing_recompute',lambda:s.apply(old,'slow')),
        ('stale_prefast_target',lambda:(s.recompute(),s.apply(old,'slow'))),
        ('historical_root_bundle',lambda:State().apply(Input('historical-actual-normal-roots',(0,2),'signed'),'signed')),
        ('additive_tail_inverse',lambda:State().apply(State().read('tail'),'sourced'))):
        try:action()
        except ValueError:rejected.append(name)
    return dict(rows=rows,rejected=rejected,checks=dict(
        canonical_bundle_and_current_pressure_targets_remain_consistent=bool(good),
        missing_recompute_and_stale_prefast_targets_rejected=len(rejected)==4,
        historical_root_mixing_and_tail_forward_inversion_rejected=len(rejected)==4,
        exact_tail_ledger_appended_not_dropped_or_double_counted=bool(good),
        metadata_model_not_materialized_field_or_continuum_certificate=True))
