#!/usr/bin/env python3
import os
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GBM_MES_DATA_ROOT", "data"))
RESULT_ROOT = Path(os.environ.get("GBM_MES_OUTPUT_ROOT", "results"))
REPO_ROOT = Path(os.environ.get("GBM_MES_REPO_ROOT", Path(__file__).resolve().parents[2]))
import numpy as np,pandas as pd
def cohen_kappa_score(a,b,labels):
    tab=pd.crosstab(pd.Categorical(a,categories=labels),pd.Categorical(b,categories=labels),dropna=False).to_numpy(float);n=tab.sum();po=np.trace(tab)/n;pe=(tab.sum(axis=1)*tab.sum(axis=0)).sum()/n**2;return (po-pe)/(1-pe) if pe<1 else np.nan
RUN=RESULT_ROOT/'provenance';V13=RESULT_ROOT/'interpretability';D=RUN/'stability'
p=pd.read_csv(V13/'03_gene_provenance/CANONICAL_95_PROVENANCE_MAP.tsv',sep='\t');p['gene']=p.gene.str.upper();layers=['BAYESPRISM','GBMAP_GSE174554_SINGLE_CELL_SINGLE_NUCLEUS','PSEUDOBULK','SPATIAL','IVY_GAP'];rows=[]
for _,r in p.iterrows():
 for layer in layers:
  deconv=pd.notna(r.malignant_median_contribution) and layer!='BAYESPRISM';gmal=r.GBmap_malignant_pseudobulk if layer!='GBMAP_GSE174554_SINGLE_CELL_SINGLE_NUCLEUS' else np.nan;smal=r.GSE174554_tumor_pseudobulk if layer!='GBMAP_GSE174554_SINGLE_CELL_SINGLE_NUCLEUS' else np.nan
  direction=r.single_cell_direction if layer!='GBMAP_GSE174554_SINGLE_CELL_SINGLE_NUCLEUS' else 'UNKNOWN'
  if deconv and r.malignant_median_contribution>=.60 and r.malignant_ci_low>.50 and direction=='MALIGNANT':cl='MALIGNANT_DOMINANT'
  elif deconv and r.ecological_median_contribution>=.60 and r.ecological_ci_low>.50 and direction=='ECOLOGICAL':cl='ECOLOGICAL_DOMINANT'
  elif deconv:cl='SHARED_MIXED_ORIGIN'
  elif pd.notna(gmal) and pd.notna(smal):cl='SHARED_MIXED_ORIGIN'
  else:cl='LOW_INFORMATION_OR_UNSTABLE'
  rows.append(dict(gene=r.gene,excluded_evidence_layer=layer,primary_class=r.primary_origin,leaveout_class=cl,retained_primary_class=cl==r.primary_origin,transition=f'{r.primary_origin}-->{cl}',algorithm_note='v1.3 thresholds unchanged; pseudo/spatial/IVY are annotation-only for primary classification'))
lo=pd.DataFrame(rows);lo.to_csv(D/'LOEO_PROVENANCE_CLASSIFICATION.tsv',sep='\t',index=False);chg=lo.groupby(['excluded_evidence_layer','primary_class','leaveout_class']).size().rename('n_genes').reset_index();chg.to_csv(D/'LOEO_CLASS_CHANGE_MATRIX.tsv',sep='\t',index=False)
order=['MALIGNANT_DOMINANT','SHARED_MIXED_ORIGIN','ECOLOGICAL_DOMINANT','LOW_INFORMATION_OR_UNSTABLE'];summary=[]
for layer,q in lo.groupby('excluded_evidence_layer'):
 a=q.primary_class;b=q.leaveout_class;k=cohen_kappa_score(a,b,labels=order);w=[]
 for x,y in zip(a,b):
  if x==y:w.append(1)
  elif {x,y}<={'MALIGNANT_DOMINANT','SHARED_MIXED_ORIGIN'} or {x,y}<={'ECOLOGICAL_DOMINANT','SHARED_MIXED_ORIGIN'}:w.append(.5)
  elif 'LOW_INFORMATION_OR_UNSTABLE' in {x,y}:w.append(.25)
  else:w.append(0)
 summary.append(dict(result_level='LEAVEOUT',excluded_evidence_layer=layer,metric='agreement',estimate=(a==b).mean(),cohen_kappa=k,weighted_agreement=np.mean(w),n_genes=len(q)))
 for cl in order:
  sens=((b[a==cl])==cl).mean() if (a==cl).any() else np.nan;ppa=((a[b==cl])==cl).mean() if (b==cl).any() else np.nan;summary.append(dict(result_level='CLASS',excluded_evidence_layer=layer,metric=cl,estimate=np.nan,cohen_kappa=np.nan,weighted_agreement=np.nan,n_genes=len(q),class_sensitivity=sens,positive_predictive_agreement=ppa))
cnt=lo.groupby('gene').retained_primary_class.sum();pol=lo.assign(polarity=lambda x:((x.primary_class=='MALIGNANT_DOMINANT')&(x.leaveout_class=='ECOLOGICAL_DOMINANT'))|((x.primary_class=='ECOLOGICAL_DOMINANT')&(x.leaveout_class=='MALIGNANT_DOMINANT'))).groupby('gene').polarity.any()
cov=pd.read_csv(V13/'03_gene_provenance/CANONICAL_95_COVERAGE.tsv',sep='\t').rename(columns={'harmonized_symbol':'gene'}).set_index('gene');thr=pd.read_csv(V13/'03_gene_provenance/CANONICAL_95_THRESHOLD_SENSITIVITY.tsv',sep='\t');cross_thr=thr.groupby('gene').origin_class.apply(lambda x:('MALIGNANT_DOMINANT'in set(x)) and ('ECOLOGICAL_DOMINANT'in set(x)))
plat=p.set_index('gene').annotations.str.contains('PLATFORM_SENSITIVE',na=False);out=[]
for _,r in p.iterrows():
 g=r.gene;coverage=bool(cov.loc[g,'GBmap_present'] and cov.loc[g,'GSE174554_present']);n=int(cnt[g]);flip=bool(pol[g] or cross_thr.get(g,False))
 if n<=1 or flip or not coverage or r.primary_origin=='LOW_INFORMATION_OR_UNSTABLE':conf='UNRESOLVED_PROVENANCE'
 elif n>=4:conf='HIGH_CONFIDENCE_PROVENANCE'
 else:conf='MODALITY_DEPENDENT_PROVENANCE'
 out.append(dict(gene=g,primary_class=r.primary_origin,retained_primary_count_of_5=n,malignant_ecological_polarity_flip=flip,coverage_sufficient=coverage,platform_sensitive=bool(plat[g]),threshold_sensitive=bool(r.threshold_sensitive_origin),final_provenance_confidence=conf,changed_layers=';'.join(lo[(lo.gene==g)&(~lo.retained_primary_class)].excluded_evidence_layer)))
final=pd.DataFrame(out);final.to_csv(D/'FINAL_PROVENANCE_CONFIDENCE.tsv',sep='\t',index=False);cs=final.final_provenance_confidence.value_counts();summary.append(dict(result_level='OVERALL',excluded_evidence_layer='ALL_FIVE',metric='completely_unchanged_genes',estimate=(final.retained_primary_count_of_5==5).sum(),n_genes=95));summary.append(dict(result_level='OVERALL',excluded_evidence_layer='ALL_FIVE',metric='malignant_to_shared_changes',estimate=((lo.primary_class=='MALIGNANT_DOMINANT')&(lo.leaveout_class=='SHARED_MIXED_ORIGIN')).sum(),n_genes=95));summary.append(dict(result_level='OVERALL',excluded_evidence_layer='ALL_FIVE',metric='ecological_to_shared_changes',estimate=((lo.primary_class=='ECOLOGICAL_DOMINANT')&(lo.leaveout_class=='SHARED_MIXED_ORIGIN')).sum(),n_genes=95));summary.append(dict(result_level='OVERALL',excluded_evidence_layer='ALL_FIVE',metric='polarity_flips',estimate=pol.sum(),n_genes=95));pd.DataFrame(summary).to_csv(D/'PROVENANCE_STABILITY_SUMMARY.tsv',sep='\t',index=False,na_rep='NA')
status='PROVENANCE_ROBUST_ACROSS_EVIDENCE_LAYERS' if cs.get('HIGH_CONFIDENCE_PROVENANCE',0)>=76 else 'PROVENANCE_PARTLY_MODALITY_DEPENDENT' if cs.get('HIGH_CONFIDENCE_PROVENANCE',0)>=48 and cs.get('UNRESOLVED_PROVENANCE',0)<=19 else 'PROVENANCE_SUBSTANTIALLY_UNSTABLE'
(D/'PROVENANCE_STABILITY_SUMMARY.md').write_text(f'''# Provenance stability analysis summary

Status: `{status}`

- Five leave-one-evidence-layer-out classifications used the unchanged v1.3 algorithm and thresholds.
- Confidence counts: HIGH={cs.get('HIGH_CONFIDENCE_PROVENANCE',0)}, MODALITY_DEPENDENT={cs.get('MODALITY_DEPENDENT_PROVENANCE',0)}, UNRESOLVED={cs.get('UNRESOLVED_PROVENANCE',0)}.
- Malignant↔ecological polarity flips: {int(pol.sum())} genes.
- Pseudo-bulk, spatial and IVY leave-outs are algorithmically unchanged because those layers supplied annotations, not the v1.3 primary class decision.
- BayesPrism and single-cell/single-nucleus leave-outs expose the expected dependence of dominant-polarity classes on those two decision layers; thresholds were not adjusted to increase stability.
''')
print('MODULE_B2_LOEO_COMPLETE',status,cs.to_dict())
