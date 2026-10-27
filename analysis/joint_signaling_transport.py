"""Sparse collective transport, equivalent to COMMOT's joint COT component.

Uses spatially sparse costs and CSR Sinkhorn updates to avoid dense all-cell
distance matrices. This is the collective component, not COMMOT's default
average of collective, row, column and pairwise solvers. The original code is
vendored for numerical equivalence checks and attribution.
"""
import time
import numpy as np
from scipy import sparse
from scipy.optimize import minimize

def collective(S,D,M,ligand_index,receptor_index,radii,eps=.1,rho=10.,maxiter=5000,tolerance=1e-6):
    start=time.monotonic();N=len(S);a=S.flatten('F');b=D.flatten('F')
    maxmass=max(a.sum(),b.sum())
    if maxmass<=0:return None,{'converged':True,'iterations':0,'edges':0,'maximum_supply_excess':0,'maximum_demand_excess':0}
    a=a/maxmass;b=b/maxmass
    rows=[];cols=[];costs=[];pair_codes=[]
    maximum_radius=max(radii);distance_scale=float(M.data[M.data<=maximum_radius].max()) if M.nnz else 1
    for p,(li,ri,radius) in enumerate(zip(ligand_index,receptor_index,radii)):
        si=np.flatnonzero(S[:,li]>0);di=np.flatnonzero(D[:,ri]>0)
        if not len(si) or not len(di):continue
        q=M[si][:,di].tocoo();take=q.data<=radius
        rows.append(si[q.row[take]]+li*N);cols.append(di[q.col[take]]+ri*N)
        costs.append(q.data[take]/max(distance_scale,1e-8));pair_codes.append(np.full(take.sum(),p,np.int16))
    if not costs:return None,{'converged':True,'iterations':0,'edges':0,'maximum_supply_excess':0,'maximum_demand_excess':0}
    row=np.concatenate(rows).astype(np.int32);col=np.concatenate(cols).astype(np.int32)
    cost=np.concatenate(costs).astype(np.float64);pair=np.concatenate(pair_codes)
    # Reduce zero-mass species/positions without changing the joint optimization.
    ai=np.flatnonzero(a>0);bi=np.flatnonzero(b>0)
    ar=np.full(len(a),-1,np.int32);br=np.full(len(b),-1,np.int32);ar[ai]=np.arange(len(ai));br[bi]=np.arange(len(bi))
    rr=ar[row];cc=br[col];aa=a[ai];bb=b[bi]
    K=sparse.coo_matrix((np.exp(-cost/eps),(rr,cc)),shape=(len(ai),len(bi))).tocsr();KT=K.T.tocsr()
    loga=np.log(aa);logb=np.log(bb);f=np.zeros(len(ai));g=np.zeros(len(bi));tiny=np.exp(-rho/eps)
    error=np.inf
    for iteration in range(100):
        oldf=f;oldg=g
        updated=eps*(loga-np.log(K@np.exp(np.clip(g/eps,-700,700))+tiny))
        f=updated
        updated=eps*(logb-np.log(KT@np.exp(np.clip(f/eps,-700,700))+tiny))
        g=updated
        if iteration%10==0:
            ef=np.max(abs(f-oldf))/max(np.max(abs(f)),np.max(abs(oldf)),1.)
            eg=np.max(abs(g-oldg))/max(np.max(abs(g)),np.max(abs(oldg)),1.)
            error=.5*(ef+eg)
            if error<tolerance:break
    # Minimize the same convex dual with diagonal mass preconditioning. This
    # avoids slow Sinkhorn drift in spatially constrained mass bottlenecks.
    sf=np.sqrt(eps/aa);sg=np.sqrt(eps/bb)
    last={'x':None,'relative_KKT':np.inf}
    def objective(z):
        ff=z[:len(ai)]*sf;gg=z[len(ai):]*sg
        # Algebraically identical sparse matrix form avoids exponentiating all
        # ligand/receptor edges on each objective evaluation.
        u=np.exp(ff/eps);v=np.exp(gg/eps)
        ss=u*(K@v);dd=v*(KT@u)
        uu=np.exp((ff-rho)/eps);vv=np.exp((gg-rho)/eps)
        value=eps*(ss.sum()+uu.sum()+vv.sum())-aa@ff-bb@gg
        gradf=ss+uu-aa
        gradg=dd+vv-bb
        last['x']=z.copy();last['relative_KKT']=max(np.max(abs(gradf)/aa),np.max(abs(gradg)/bb))
        return value,np.r_[gradf*sf,gradg*sg]
    def callback(xk):
        if last['x'] is None or not np.array_equal(xk,last['x']):objective(xk)
        if last['relative_KKT']<=max(tolerance,1e-8):raise StopIteration
    initial=np.r_[f/sf,g/sg]
    upper=np.r_[(rho+eps*loga)/sf,(rho+eps*logb)/sg]
    result=minimize(objective,initial,method='L-BFGS-B',jac=True,callback=callback,bounds=list(zip(np.full(len(initial),-np.inf),upper)),
        options={'maxiter':int(maxiter),'maxfun':int(maxiter*4),'maxls':40,'maxcor':15,'ftol':1e-15,'gtol':1e-11})
    f=result.x[:len(ai)]*sf;g=result.x[len(ai):]*sg
    # Polish resource equations with ordinary Sinkhorn, keeping the objective unchanged.
    for polish in range(2001):
        supply_norm=np.exp(f/eps)*(K@np.exp(g/eps)+tiny)
        demand_norm=np.exp(g/eps)*(KT@np.exp(f/eps)+tiny)
        error=max(np.max(abs(supply_norm-aa)/aa),np.max(abs(demand_norm-bb)/bb))
        if error<=max(tolerance,1e-8):break
        f=eps*(loga-np.log(K@np.exp(g/eps)+tiny))
        g=eps*(logb-np.log(KT@np.exp(f/eps)+tiny))
    mass=np.exp((f[rr]+g[cc]-cost)/eps)*maxmass
    supply=np.bincount(rr,weights=mass,minlength=len(ai));demand=np.bincount(cc,weights=mass,minlength=len(bi))
    excess_s=float(np.max(np.maximum(supply-aa*maxmass,0))/max(np.max(aa*maxmass),1e-8))
    excess_d=float(np.max(np.maximum(demand-bb*maxmass,0))/max(np.max(bb*maxmass),1e-8))
    diagnostics={'converged':bool(error<=max(tolerance,1e-8)),'iterations':int(result.nit),'warmup_sinkhorn_iterations':iteration+1,
        'polish_sinkhorn_iterations':polish,'maximum_relative_KKT_residual':float(error),'optimizer_success':bool(result.success),
        'optimizer_message':str(result.message),'edges':len(mass),'solver':'Mass-preconditioned L-BFGS on original joint-COT dual, with Sinkhorn polish',
        'transported_mass':float(mass.sum()),'supply_mass':float(S.sum()),'demand_mass':float(D.sum()),
        'maximum_supply_excess':excess_s,'maximum_demand_excess':excess_d,
        'maximum_relative_supply_excess':float(np.max(np.maximum(supply-aa*maxmass,0)/np.maximum(aa*maxmass,1e-12))),
        'maximum_relative_demand_excess':float(np.max(np.maximum(demand-bb*maxmass,0)/np.maximum(bb*maxmass,1e-12))),
        'seconds':time.monotonic()-start}
    return {'source_cell':row%N,'receiver_cell':col%N,'pair_index':pair,'mass':mass},diagnostics

def equivalence_check():
    from scipy.spatial.distance import cdist
    from vendor_commot._cot import cot_sparse
    rng=np.random.default_rng(178);N=24;xy=rng.uniform(0,240,(N,2))
    S=rng.gamma(1,2,(N,3));D=rng.gamma(1,2,(N,3));S[rng.random(S.shape)<.3]=0;D[rng.random(D.shape)<.3]=0
    A=np.full((3,3),np.inf);li=np.array([0,0,1,2]);ri=np.array([0,1,1,2]);radii=np.array([200,50,200,50])
    for l,r in zip(li,ri):A[l,r]=1
    cutoff=np.full((3,3),50.);cutoff[li,ri]=radii
    dm=cdist(xy,xy);np.fill_diagonal(dm,np.inf)
    old=cot_sparse(S,D,A,dm,cutoff,nitermax=15000,stopthr=1e-10)
    rr,cc=np.where(dm<=200);M=sparse.csr_matrix((dm[rr,cc],(rr,cc)),shape=(N,N))
    new,diag=collective(S,D,M,li,ri,radii,maxiter=15000,tolerance=1e-10)
    errors=[]
    for p,(l,r) in enumerate(zip(li,ri)):
        take=new['pair_index']==p
        q=sparse.coo_matrix((new['mass'][take],(new['source_cell'][take],new['receiver_cell'][take])),shape=(N,N))
        errors.append(np.max(abs(q.toarray()-old[l,r].toarray())))
    assert max(errors)<1e-7,(errors,diag)
    assert diag['converged'] and diag['maximum_supply_excess']<1e-4 and diag['maximum_demand_excess']<1e-4
    return {'test':'Sparse collective solver versus original COMMOT cot_sparse','maximum_absolute_error':float(max(errors)),**diag}

if __name__=='__main__':
    import json
    print(json.dumps(equivalence_check(),indent=2))
