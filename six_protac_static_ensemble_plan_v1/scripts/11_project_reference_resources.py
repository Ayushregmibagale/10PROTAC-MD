"""Reference-based scenarios, not a production capacity certificate."""
import argparse,csv,json
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--resources',type=Path,required=True);p.add_argument('--transfer',type=Path,required=True)
    p.add_argument('--pairs',type=int,default=3600);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.pairs<1:raise ValueError('Positive pair count required')
    report=json.loads(a.resources.read_text());transfer=json.loads(a.transfer.read_text())
    if report['status']!='SOURCE_INPUT_OUTPUT_INTEGRITY_PASS':raise ValueError('Resource receipts not verified')
    denominator=sum(x['run_logical_bytes'] for x in report['runs'])
    ratio=transfer['archive_bytes']/denominator
    rows=[]
    for run in report['runs']:
        if not run['name'].startswith('global_'):continue
        for concurrent in (1,2,4,8,16,24):
            rows.append({'reference':run['reference'],'nominal_pairs':a.pairs,'concurrent_pairs':concurrent,
             'allocated_CPU_slots':run['cpu_affinity_limit']*concurrent,
             'reference_based_wall_days':a.pairs*run['wall_seconds']/86400/concurrent,
             'raw_logical_storage_TB_decimal':a.pairs*run['run_logical_bytes']/1e12,
             'archive_storage_TB_decimal_at_observed_ratio':a.pairs*run['run_logical_bytes']*ratio/1e12,
             'summed_peak_process_RSS_GB_decimal':run['peak_summed_group_RSS_bytes']*concurrent/1e9})
    a.output.mkdir(parents=True,exist_ok=True)
    with (a.output/'REFERENCE_RESOURCE_SCENARIOS.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result={'status':'REFERENCE_SCENARIOS_ONLY_NOT_RELEASED','nominal_pairs':a.pairs,
      'observed_archive_ratio':ratio,'source_archive_sha256':transfer['archive_sha256'],
      'rows':rows,'limitations':['Two bound-reference systems at50global×10local×20conformers, one seed each',
       'Other kinases/PROTACs may need more time/storage; these are observed reference scenarios, not worst-case bounds',
       'CPU slots are affinity allocations, not measured CPU-time utilization','No reserved96CPU entitlement is inferred',
       'No scheduler delays, preprocessing, replica/seed expansion, backups, feature extraction or GNN training included',
       'Archive ratio is from a combined gzip of these raw references; per-pair compression and future datasets can differ',
       'Uncompressed staged data plus archives require additional simultaneous storage; full raw data retention is preserved losslessly',
       'The local500×1×20comparator uses different hardware; wall-time differences are not a paired speed benchmark',
       'Nominal3600pairs reflects requested6×600scope, not3600currently qualified input pairs'],
       'MD':False,'production_submission_authorized':False}
    (a.output/'REFERENCE_RESOURCE_SCENARIOS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'archive_ratio':ratio,'one_pair_at_a_time_days':[{k:v for k,v in row.items() if k not in ('summed_peak_process_RSS_GB_decimal',)} for row in rows if row['concurrent_pairs']==1]},indent=2))

if __name__=='__main__':main()
