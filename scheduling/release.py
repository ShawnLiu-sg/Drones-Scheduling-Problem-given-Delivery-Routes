"""Fixed-plan release uncertainty on a user-supplied, physically validated schedule."""
import hashlib
import math
import numpy as np
from scipy.special import ndtr, ndtri, gammainc, gammaincinv
from scipy.stats import laplace

DISTS=("uniform","normal","laplace","gamma")
KAPPAS=(0.05,0.10,0.20)
N=10000
SEED=20261124
def uniform_quantile(u,law):
    if law=="uniform": return (2*u-1)*np.sqrt(3)
    if law=="normal": return ndtri(u)
    if law=="laplace": return laplace.ppf(u,scale=1/np.sqrt(2))
    if law=="gamma": return (gammaincinv(4,u)-4)/2
    raise ValueError(law)
def cdf(x,law):
    if law=="uniform": return np.clip((x/np.sqrt(3)+1)/2,0,1)
    if law=="normal": return ndtr(x)
    if law=="laplace": return laplace.cdf(x,scale=1/np.sqrt(2))
    if law=="gamma": return gammainc(4,np.maximum(4+2*x,0))
    raise ValueError(law)
def prob(s,scale,law):
    return float(np.exp(np.log(cdf(s/scale,law)).sum()))
def released(r,scale,e):
    return np.maximum(0,np.ceil(r+scale*e))
def wilson(h,n):
    z=1.959963984540054;p=h/n
    m=(p+z*z/(2*n))/(1+z*z/n)
    d=z*math.sqrt((p*(1-p)+z*z/(4*n))/n)/(1+z*z/n)
    return (0 if h==0 else m-d,1 if h==n else m+d)

def evaluate(instance, plan, identifier):
    """Common random numbers across distributions and strengths, with exact CDF checks."""
    releases, departures, batches = [], [], []
    for k, (route, schedule) in enumerate(zip(instance["routes"], plan, strict=True)):
        starts = {int(i): t for i, t in schedule["starts"].items()}
        for parcel, site in zip(route["parcels"], schedule["assign"], strict=True):
            releases.append(parcel["r"])
            departures.append(starts[site])
            batches.append((k,site))
    r, R = np.array(releases), np.array(departures)
    if any(t != int(t) for t in departures):
        raise ValueError("Release uncertainty expects integer ticks in the replayed plan.")
    seed = int.from_bytes(hashlib.sha256(f"{SEED}|{identifier}".encode()).digest()[:16], "big")
    uniform = np.random.Generator(np.random.PCG64(seed)).random((N,len(r)))
    uniform = np.maximum(uniform, np.nextafter(0.,1.))
    records = []
    for law in DISTS:
        residual = uniform_quantile(uniform,law)
        for kappa in KAPPAS:
            scale=np.maximum(1,kappa*r)
            actual=released(r,scale,residual)
            passed=(actual<=R)
            count=int(passed.all(axis=1).sum())
            lo,hi=wilson(count,N)
            failed={f"{k}:{site}":int((~passed[:,[i for i,b in enumerate(batches) if b==(k,site)]].all(axis=1)).sum())
                    for k,site in sorted(set(batches))}
            records.append(dict(distribution=law,kappa=kappa,samples=N,pass_count=count,
                                probability=count/N,analytical_probability=prob(R-r,scale,law),
                                wilson95=[lo,hi],batch_failures=failed))
    return dict(seed=seed, package_slacks=(R-r).tolist(), metrics=records)
