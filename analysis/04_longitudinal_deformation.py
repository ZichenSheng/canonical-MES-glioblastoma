"""Frozen patient-exclusive all-k deformation arithmetic on prepared paired scores.
Supported dimensions are frozen inputs; this module does not rerun their selection,
finite-cell null generation, or component decomposition.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
from scipy.stats import beta
EPS=1e-12
def zfit(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mu, sd = x.mean(0), x.std(0, ddof=1)
    valid = np.isfinite(mu) & np.isfinite(sd) & (sd > EPS)
    return (x[:, valid]-mu[valid])/sd[valid], mu, sd, valid
def orient_loadings(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float).copy()
    for j in range(v.shape[1]):
        idx = int(np.argmax(np.abs(v[:, j])))
        if v[idx, j] < 0:
            v[:, j] *= -1
    return v

def run(xp,xr,label):
    n,p=xp.shape
    if xr.shape!=xp.shape or not np.isfinite(xp).all() or not np.isfinite(xr).all():raise ValueError("Paired finite matrices required")
    kmax=min(20,n-2,p,n-1)
    xp_df=pd.DataFrame(xp)
    spec=SimpleNamespace(dataset_id=label,label=label)
    deformation_rows=[]

    for i, pid in enumerate(xp_df.index):
        mask=np.arange(n)!=i
        trainz,fmu,fsd,fvalid=zfit(xp[mask])
        pv=(xp[i,fvalid]-fmu[fvalid])/fsd[fvalid]
        rv=(xr[i,fvalid]-fmu[fvalid])/fsd[fvalid]
        delta=rv-pv
        denergy=float(delta@delta)
        _,s,vt=np.linalg.svd(trainz,full_matrices=False)
        vfull=orient_loadings(vt[:kmax].T)
        eigfold=(s**2)/max(len(trainz)-1,1)
        for k in range(1,kmax+1):
            v=vfull[:,:k]
            parallel=v@(v.T@delta); orth=delta-parallel
            rp=pv-v@(v.T@pv); rr=rv-v@(v.T@rv)
            pe=float(parallel@parallel); oe=float(orth@orth)
            frac=pe/denergy if denergy>EPS else np.nan
            random_med=float(beta.ppf(.5,k/2,(len(delta)-k)/2)) if k<len(delta) else 1.0
            random_p=float(1-beta.cdf(frac,k/2,(len(delta)-k)/2)) if k<len(delta) else np.nan
            deformation_rows.append(dict(dataset_id=spec.dataset_id,dataset_label=spec.label,patient_id=pid,k=k,
                training_primary_n=n-1,feature_n=len(delta),parallel_energy=pe,orthogonal_energy=oe,total_delta_energy=denergy,
                parallel_energy_fraction=frac,orthogonal_energy_fraction=oe/denergy if denergy>EPS else np.nan,
                heldout_primary_normalized_residual=float(rp@rp/max(pv@pv,EPS)),recurrent_normalized_residual=float(rr@rr/max(rv@rv,EPS)),
                recurrent_minus_heldout_primary_residual=float(rr@rr/max(rv@rv,EPS)-rp@rp/max(pv@pv,EPS)),
                reference_explained_variance_fraction=float(eigfold[:k].sum()/eigfold.sum()),
                isotropic_random_reference_median_fraction=random_med,isotropic_random_reference_p_one_sided=random_p,
                real_minus_random_reference_median_fraction=frac-random_med,reference_fit="LOPO_PRIMARY_ONLY"))

    d=pd.DataFrame(deformation_rows)
    return {str(k):float(q.orthogonal_energy_fraction.median()) for k,q in d.groupby("k")}

if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("input");p.add_argument("--output",required=True)
    a=p.parse_args();z=np.load(a.input,allow_pickle=False)
    result={label:run(z[label+"_primary"],z[label+"_recurrent"],label) for label in ["CARE","GLASS","GSE174554"]}
    Path(a.output).write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result))
