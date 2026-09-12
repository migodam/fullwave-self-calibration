"""Finite class-C competing-world search, explicitly NOT a lower certificate.

Both worlds use the SAME fixed covariance. Values are upper-distance numerical
candidates; passivity residual radii are NOT yet available for these solutions.
"""
from run_experiment import *
from scipy.special import ndtr


def main():
    target=HERE/'results/finite_pair.json'
    if target.exists():raise FileExistsError(target)
    rx=mx.receivers(12,.6)
    A=np.array([3.8,4.8,.4]);gA=np.exp(.2j)
    model=Cluster(CENTERS,RADII,K,mx.illuminations(),order=3)
    def field(m,t):return m.field(np.array(t[:2])+1j*LOSS,rx+np.array([t[2]*.001,0,0])).ravel()
    fA=field(model,A);sigma=.01*np.linalg.norm(gA*fA)/np.sqrt(len(fA))
    output=[]
    # eps_2 differs by .21, so the two .1-success boxes are disjoint.
    for ref in (False,True):
        def residual(v,m=model):
            B=np.array([v[0],4.59,v[1]]);gB=v[2]*np.exp(1j*v[3])
            d=(gA*field(m,A)-gB*field(m,B))/sigma
            if ref:d=np.r_[d,(gA-gB)/.01]
            return np.sqrt(2)*np.r_[d.real,d.imag]
        fits=[least_squares(residual,s,bounds=([1.5,-2,.75,-np.pi],[4,2,1.25,np.pi]),
             ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=150)
             for s in ([3.6,.4,1.02,.2],[3.8,0,1,.2],[3.2,-1,1.1,.25])]
        best=min(fits,key=lambda v:v.fun@v.fun);v=best.x
        B=np.array([v[0],4.59,v[1]]);gB=v[2]*np.exp(1j*v[3])
        ds=[]
        for L in (3,5,7):
            m=model if L==3 else Cluster(CENTERS,RADII,K,mx.illuminations(),order=L,nt=24,nphi=48)
            ds.append({'L':L,'whitened_real_distance':float(np.linalg.norm(residual(v,m)))})
        d=ds[-1]['whitened_real_distance']
        output.append({'with_reference':ref,'world_A':A.tolist(),'gain_A':encode(gA),
            'world_B':B.tolist(),'gain_B':encode(gB),'complex_sigma_shared':sigma,
            'complex_reference_sigma':.01 if ref else None,'distances':ds,
            'binary_risk_at_numerical_mean':float(ndtr(-d/2)),
            'allowed_combined_whitened_mean_error_for_1percent_risk':float(-2*scipy.special.ndtri(.01)-d),
            'starts':[{'x':v.x.tolist(),'cost':float(v.fun@v.fun),'status':int(v.status)} for v in fits],
            'status':'realistic numerical candidate, NOT certified continuum counterexample'})
    # Explicit passive Maxwell residual bound constants for this domain.
    alpha=min(.03/(3**2+.03**2),.05/(4**2+.05**2))
    volume=4*np.pi/3*sum(RADII**3);dmin=.6-.06-.035-.002
    S=np.sqrt(12*volume/(16*np.pi**2*dmin**2)*(2*K**4+2*K**2/dmin**2+6/dmin**4))
    result={'scope':'exploratory finite candidates, no covering lower bound',
        'pairs':output,'passivity_constants_diagnostic':{'alpha':alpha,'radius_min':dmin,
            'S_norm_Hilbert_Schmidt_upper':S,'output_residual_amplification_gmax_S_over_alpha':1.25*S/alpha},
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
