"""Full coverage certificate for the restricted lossless A5 class M.

Every epsilon endpoint is an exact rational. Size parameters are enclosed from
1/5 and 3/20, not Python floats. All exported decimal bounds round outward.
No endpoint-minimum assertion and no midpoint-Lipschitz shortcut are used.
"""
import argparse, hashlib, json, platform, time
from pathlib import Path
from dyadic_interval import I, BITS, SCALE, modal

HERE=Path(__file__).resolve().parent

def certify(out, step_denominator=2000, size_relative_ppm=0):
    if out.exists():
        raise FileExistsError(out)
    start=time.perf_counter(); rows=[]
    for name,begin,end,xp,xq in [('sphere1',3000,8000,1,5),('sphere2',4000,10000,3,20)]:
        # Here begin/end are specified in units 1/2000, fixed exact cover.
        if step_denominator != 2000:
            raise ValueError('Coverage endpoints are registered for denominator 2000')
        ppm=size_relative_ppm
        x=I.box(xp*(1000000-ppm),xq*1000000,xp*(1000000+ppm),xq*1000000)
        totals={}; digest=hashlib.sha256(); failures=[]
        for k in range(begin,end):
            epsilon=I.box(k,2000,k+1,2000)
            quantities=modal(epsilon,x)
            for key,val in quantities.items():
                old=totals.get(key,val)
                totals[key]=I(min(old.lo,val.lo),max(old.hi,val.hi))
            if quantities['q'].lo<=0 or quantities['qprime'].lo<=0 or quantities['denominator'].lo<=0:
                failures.append(k)
            digest.update((str(k)+json.dumps({key:v.record() for key,v in quantities.items()},sort_keys=True)).encode())
        row={'name':name,'epsilon_domain':[[begin,2000],[end,2000]],
             'x_rational':[xp,xq],'size_relative_ppm':ppm,'x_enclosure':x.record(),
             'coverage':{'first_index':begin,'exclusive_last_index':end,'denominator':2000,
                         'count':end-begin,'all_consecutive':True},
             'bounds':{k:v.record() for k,v in totals.items()},
             'per_box_digest_sha256':digest.hexdigest(),'nonpositive_boxes':failures}
        rows.append(row)
    result={'kind':'integer_outward_interval_cover','bits':BITS,
            'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'dyadic_interval.py']},
            'python':platform.python_version(),'wall_seconds':time.perf_counter()-start,
            'status':'certified_restricted_modal' if all(not r['nonpositive_boxes'] for r in rows) else 'unresolved',
            'claims_not_made':['hardware accuracy','class C uniform recovery','endpoint attainment of derivative infimum'],
            'rows':rows}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=HERE/'results/modal_certificate.json');p.add_argument('--size-relative-ppm',type=int,default=0)
    args=p.parse_args();certify(args.out,size_relative_ppm=args.size_relative_ppm)
