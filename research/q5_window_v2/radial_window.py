"""Exact radial-window audit: ideal angular electric-dipole readout only.

Deterministic analysis, not a new acquisition-policy experiment. Raw projected
error is fixed at 1e-6; it is NOT the inherited normalized Gaussian noise.
"""
from pathlib import Path
import json
from modal_interval import Q,iv,rat,bounds,export,enclose

def c2(z):
    z=Q(z)
    if z<=0: raise ValueError('positive radial size required')
    return 1/z**2-1/z**4+1/z**6

def main():
    iv.dps=60;results=[];Ba=Q(3291,10**6);raw=Q(1,10**6)
    for x,hi in [(Q(1,5),Q(4)),(Q(3,20),Q(5))]:
        H=enclose(rat(hi),rat(x))['h'];hm=enclose(rat(hi-Q(1,5)),rat(x))['h']
        C=rat(Q(3,4))*(H-hm)-rat(2*Ba)*hm
        lowC,highC=bounds(C);assert lowC>0
        threshold_low=(2*raw/highC)**2;threshold_high=(2*raw/lowC)**2
        lo=Q(1,5);hiR=Q(64)
        assert c2(lo)>threshold_high and c2(hiR)<threshold_low
        for _ in range(55):
            mid=(lo+hiR)/2;v=c2(mid)
            if v>=threshold_high: lo=mid
            elif v<=threshold_low:hiR=mid
            else:break
        assert c2(lo)>=threshold_high and c2(hiR)<=threshold_low
        results.append({'size':str(x),'finite_C':export(C),'root_enclosure':[str(lo),str(hiR)],'root_display':[float(lo),float(hiR)]})
    result={'scope':'ideal full-angle separately read electric dipole; known geometry; amplitude box experiment only',
        'raw_projected_error':str(raw),'real_reference_error':str(Ba),'formula':'c(z)^2=z^-2-z^-4+z^-6',
        'derivative':'-2*((z^2-1)^2+2)/z^7 < 0 for z>0','results':results,
        'no_interior_optimum':True,'no_class_C_transfer':True}
    (Path(__file__).parent/'results/radial_window.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
