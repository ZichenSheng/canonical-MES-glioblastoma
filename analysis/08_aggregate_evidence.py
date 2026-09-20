"""Recalculate manuscript summaries from safe frozen aggregate inputs.
No patient rows, source expression, null generation or figure assembly is implied.
"""
from pathlib import Path
import argparse,json,itertools
from decimal import Decimal
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
def read(name):return pd.read_csv(ROOT/'results_reference'/name,sep='\t')
def close(actual,expected,tolerance=1e-12):
    if not np.isfinite(actual) or abs(Decimal(str(actual))-Decimal(str(expected)))>Decimal(str(tolerance)):
        raise ValueError(f'Aggregate mismatch: {actual} vs {expected} (+/-{tolerance})')
def run():
    out={}
    m=pd.read_csv(ROOT/'resources/signatures/Primary201_canonical.tsv',sep='\t')
    sets=m.groupby('model_instance_id').gene_symbol.agg(set)
    j=np.array([len(a&b)/len(a|b) for a,b in itertools.combinations(sets,2)])
    out['zero_overlap']=int((j==0).sum());out['pair_n']=len(j);out['median_jaccard']=float(np.median(j))
    close(out['zero_overlap'],19370,0);close(len(j),20100,0);close(np.median(j),0,0)
    a=read('matched_panel_summaries.tsv')
    a=a[a.world_type=='RANDOM_MATCHED_PANEL']
    if len(a)!=1000:raise ValueError('Expected1000 matched panels')
    out['matched_p_structure']=float((1+(a.score_correlation_structure>=.7309225783071415).sum())/1001)
    out['matched_p_CGGA693_variance']=float((1+(a.PC1_variance_CGGA693>=.5804900544696167).sum())/1001)
    close(out['matched_p_structure'],1/1001);close(out['matched_p_CGGA693_variance'],1)
    a=read('localization_classes.tsv');out['localization']=a.final_class.value_counts().to_dict()
    if out['localization']!={'malignant-dominant':130,'TME-dominant':35,'mixed':16,'unstable':20}:raise ValueError('Localization mismatch')
    if not (a.evaluable_patient_n==59).all():raise ValueError('Paired patient denominator mismatch')
    a=read('ivy_model_summary.tsv');rho=float(spearmanr(a.partial_r2,a.balanced_accuracy).statistic)
    close(rho,.804,.0005);out['ivy_rho']=rho
    a=read('geomx_model_summary.tsv');close(len(a),198,0);close(a.balanced_accuracy.median(),.602040816327)
    out['geomx_ba']=float(a.balanced_accuracy.median())
    if not ((a.patient_n==26)&(a.observation_n==98)).all():raise ValueError('GeoMx denominator mismatch')
    a=read('residual_mode_gates.tsv')
    passes=(a.eigenvalue>a.permutation_eigenvalue_95)&(a.bootstrap_loading_cosine_ci_low>=.60)&(a.loo_loading_cosine_min>=.60)
    close(len(a),10,0);close(passes.sum(),0,0);out['stable_modes']=int(passes.sum())
    a=read('dimension_support.tsv');a=a[a.method=='supported_range_union_prespecified_methods']
    out['supported_dimensions']={r.dataset_label:[int(r.estimate_low),int(r.estimate_high)] for r in a.itertuples()}
    if out['supported_dimensions']!={'CARE':[7,16],'GLASS':[3,12],'GSE174554':[2,11]}:raise ValueError('Dimension mismatch')
    a=read('null_calibration.tsv');a=a[a.q.notna()]
    for r in a.itertuples():
        close(r.blocked_rejections/r.null_world_n,r.blocked_size,5e-12)
        close(r.naive_rejections/r.null_world_n,r.naive_size,5e-12)
    out['null_calibration_rows']=len(a)
    a=read('reference_comparison.tsv');cols=['naive_only','blocked_only','both','neither']
    if not (a[cols].sum(axis=1)==a.world_n).all():raise ValueError('World partition mismatch')
    out['historical_world_n']=int(a.world_n.sum());close(out['historical_world_n'],640,0)
    out['historical_counts']=a[cols].sum().astype(int).tolist()
    if out['historical_counts']!=[140,114,4,382]:raise ValueError('Historical categories mismatch')
    q=a[(a.q==1)&a.alpha_reachable];out['formal_q100_counts']=q[cols].sum().astype(int).tolist()
    if out['formal_q100_counts']!=[51,56,4,29] or q.world_n.sum()!=140:raise ValueError('Formal q100 mismatch')
    out['protagonists']={}
    for id,expected in [('PF4C-B4-007',.003),('PF4C-B4-029',.107),('PF4C-FINAL-079',.240)]:
        r=q[q.signature_id==id].iloc[0];close(r.blocked_p_median,expected,.0005)
        out['protagonists'][id]={'blocked_median_p':float(r.blocked_p_median),'blocked_recovery':int(r.blocked_only+r.both),'naive_recovery':int(r.naive_only+r.both)}
    a=read('crc04_inference.tsv').iloc[0]
    close(a.N_PATIENTS,9,0);close(a.N_SEGMENTS,85,0);close(a.TARGET_GENE_COUNT,11,0)
    close(a.EXACT_P_VALUE,16/512,0);close(a.MIN_ATTAINABLE_P,1/512,0)
    close(a.OBSERVED_TEST_STATISTIC,1.596,.0005)
    a=read('crc04_reproduction_summary.tsv');failed=set(a.loc[~a.passed,'check'])
    if failed!={'segment_scores','core_patient_region_scores','invasive_front_patient_region_scores','paired_differences'}:raise ValueError('Historical CRC04 discrepancy changed')
    out['crc04_intermediate_failures']=len(failed)
    a=read('representation_publications.tsv')
    close(a.model_instance_id.nunique(),201,0);close(a.publication_id.nunique(),197,0)
    a=read('supplement_S1.tsv')
    close(float(a.iloc[2]['estimate']),73,0);close(float(a.iloc[3]['estimate']),108,0)
    if a.iloc[14]['estimate']!='structure concordance = 0.691':raise ValueError('Overlap-removed structure changed')
    a=read('spatial_axis_summary.tsv')
    if len(a)!=3 or not (pd.to_numeric(a.value)>.8).all():raise ValueError('Spatial common-axis summary changed')
    a=read('care_transfer_baseline.tsv').iloc[0]
    close(a.median_mean_cv_r2,.09031,1e-12);close(a.patient_n,56,0);close(a.repeat_n,100,0);close(a.fold_n,5,0)
    a=read('matched_score_descriptive.tsv')
    if len(a)!=4 or not ((a.median_ratio>0)&(a.median_ratio<1)).all():raise ValueError('Matched-score summary changed')
    if not (a.status=='DEPRECATED_AS_INFERENTIAL_EVIDENCE').all():raise ValueError('Do not promote matching to inference')
    a=read('source_change_descriptive.tsv')
    if not (a.value.min()<0<a.value.max()):raise ValueError('Source-change sign summary changed')
    a=read('matched_ecology_descriptive.tsv').iloc[0]
    close(a.matched_record_n,171,0);close(a.positive_ecology_distance_n,171,0)
    close(a.median_score_distance,.0122222886281552,1e-12);close(a.median_ecology_distance,2.3949644758960313,1e-12)
    a=read('reference_examples.tsv').set_index('signature_id')
    for id,naive,blocked in [('PF4C-B4-007',.172,.002),('PF4C-FINAL-079',.025,.123)]:
        close(a.loc[id,'naive_p'],naive,0);close(a.loc[id,'blocked_p'],blocked,0)
    a=read('protagonist_evaluability.tsv');a=a[(a.dataset=='Ivy')&(a.claim_scope=='Malignant-cell-specific MES state')]
    if len(a)!=3 or not (a.endpoint_compatibility=='NO').all() or not (a.score_computable=='YES').all():raise ValueError('Claim-scope boundary changed')
    a=read('claim_records.tsv')
    if len(a)!=4 or a.context.nunique()!=4 or a.signature_id.nunique()!=1:raise ValueError('Claim records changed')
    out['additional_current_claim_checks']='PASS: universe,cohorts,overlap-control,spatial,CARE,descriptive matching,reference examples,claim records'
    a=read('cohort_geometry.tsv').set_index(['dataset_label','geometry_component'])
    for label,translation,angle,trace in [('CARE',3.329,76.6,1.012),('GLASS',3.233,32.3,1.349),('GSE174554',4.784,85.1,.531)]:
        close(a.loc[(label,'centroid_translation'),'estimate'],translation,.0005)
        close(a.loc[(label,'subspace_rotation_max_principal_angle_deg'),'estimate'],angle,.05)
        close(a.loc[(label,'variance_expansion_topk_trace_ratio'),'estimate'],trace,.0005)
    a=read('supplement_S5.tsv')
    if len(a)!=3 or not (a.single_cell_localization=='malignant-dominant').all():raise ValueError('Protagonist localization changed')
    if not a.mathematically_equivalent_to_original.str.startswith('NO').all():raise ValueError('Original predictors are not equivalent')
    if list(a.DS03_blocked_recovery_q100)!=['20/20','0/20','0/20'] or list(a.DS03_naive_recovery_q100)!=['0/20','0/20','20/20']:raise ValueError('Protagonist recovery changed')
    if not a.spatial_support.str.contains('matrix-rich:').all() or not a.spatial_support.str.contains('/17 patients positive').all():raise ValueError('Spatial provenance missing')
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output');args=p.parse_args();result=run()
    if args.output:Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: aggregate evidence checks.');print(json.dumps(result))
