"""Current first-real weighted wave composition, conditional on native providers."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
from functools import lru_cache
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q92_full_linear_envelope/q92.py'
s=importlib.util.spec_from_file_location('q93_linear',p);linear=importlib.util.module_from_spec(s);s.loader.exec_module(linear)
full=linear.full;EP=full.EP;convert=full.convert

def rescale(v):return EP({k:c.scale(Q(25,16)**k) for k,c in convert(v).c.items()})
def clone_wave(root,wave):
    return type(root)(H=root.H,target=root.target,wave=wave,curl=root.curl,loghalf=root.loghalf,
        K0=root.K0,beta=root.beta,R=root.R,mass=root.mass,width=root.width,trace=root.trace)

class CurrentWaves:
    def __init__(self,postslow,*,regular_wave,phase,inverse_normal_square):
        self.means=linear.meanwave.means.FinalMeans(postslow);self.state=self.means.state
        self.root=self.state.target.root
        if self.root is not self.state.primary.primary.root:raise ValueError('Original root identity required')
        if not callable(regular_wave):raise TypeError('Original regular pulse norm provider required')
        self.regular=clone_wave(self.root,regular_wave);self.cut=self.root;self.phase=phase
        if any(not callable(getattr(phase,k,None)) for k in ('normal','matrix','moving_normal')):
            raise TypeError('Full native normal, matrix and moving-normal jets required')
        if isinstance(inverse_normal_square,(float,bool)):raise ValueError('Exact positive normal floor required')
        self.q=Q(inverse_normal_square)
        if self.q<=0:raise ValueError('Exact positive normal floor required')
    @lru_cache(None)
    def signed(self,j,cut):
        root=self.cut if cut else self.regular;target=self.state.target
        def V(b):return sum((convert(root.K(a))*convert(target.total(b-a)) for a in range(b+1)),EP())
        return root.mass*Q(1,2)*EP({1:1})*sum((convert(root.E(a))*V(b)*convert(root.reciprocal(c))*convert(root.data('wave',j-a-b-c))
            for a in range(j+1) for b in range(j-a+1) for c in range(j-a-b+1)),EP())
    @lru_cache(None)
    def principal(self,j,cut):return convert((self.cut if cut else self.regular).amplitude(j))+self.signed(j,cut)
    @lru_cache(None)
    def curl(self,j,cut):
        root=self.cut if cut else self.regular
        return sum((convert(root.data('curl',i))*(self.principal(j-i,cut)+(j-i+1)*self.principal(j-i+1,cut)) for i in range(j+1)),EP())
    def jet(self,name,j):return convert(getattr(self.phase,name)(j))
    def convolution(self,a,b,j):return sum((comb(j,i)*a(i)*b(j-i) for i in range(j+1)),EP())
    @lru_cache(None)
    def inverse(self,j):
        if not j:return EP(self.q)
        normal=lambda i:self.jet('normal',i)
        return self.q*sum((comb(j,i)*self.convolution(normal,normal,i)*self.inverse(j-i) for i in range(1,j+1)),EP())
    def row(self,j):
        numerator=lambda i:self.convolution(lambda k:self.jet('normal',k),lambda k:self.jet('matrix',k),i)+self.jet('moving_normal',i)
        return self.convolution(self.inverse,numerator,j)
    def through(self,N):
        full.order(N)
        if linear.meanwave.means.FinalMeans(self.means.postslow).state is not self.state:raise ValueError('Current pressure state changed')
        if N>self.state.max_output:raise ValueError('Outside original current-state reserve')
        A=sum((self.principal(j,False) for j in range(N+1)),EP())
        C=sum((self.curl(j,False) for j in range(N+1)),EP())
        U=sum((self.principal(j,True)+self.curl(j,True) for j in range(N+1)),EP())
        row=sum((Q(5**j,factorial(j))*self.row(j) for j in range(N+1)),EP())
        normal=sum((Q(5**j,factorial(j))*self.jet('normal',j) for j in range(N+1)),EP())
        factor=2*convert(self.state.primary.chart.transfer(N))*16*188**N*factorial(N)
        vals={k:factor*rescale(v) for k,v in dict(T=A,R=C,Pi=row*A,U=U,normal=normal).items()}
        return dict(vals,order=N,pole=max(6*N+9,6*N+self.state.target.pole+6),
                    target_request=N+1,wave_state=self.state,chart='factorial-through-R',
                    actual_inputs_certified=False)

def record():
    return dict(status='current-target-original-root-weighted-wave-composition',
        independent_T_R_Pi_U_callbacks_removed=True,regular_and_cut_norm_slots_separate=True,
        current_target_object_retained=True,coarse_target_rate='2-3*kappa',
        Q51_stronger_coefficient_premise_required_for_this_coarse_bound=False,
        Gaussian_is_per_original_label=True,actual_native_provider_inequalities_certified=False,
        full_source_complete_residual_or_blowup=False)
