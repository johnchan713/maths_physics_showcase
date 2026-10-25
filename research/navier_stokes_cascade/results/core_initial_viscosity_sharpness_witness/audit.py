import argparse,hashlib,json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=('core_initial_final_row_obstruction','real_signed_covariance_update','post_signed_mean_relaxation')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
    bound={};checks={}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,h in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==h
    r=review();paths=[HERE/n for n in ('README.md','witness.py','controls.py','audit.py','test_audit.py','requirements.txt','.gitignore')]
    paths += [HERE.parent/'core_initial_final_row_obstruction/bounds.py',HERE.parent/'real_signed_covariance_update/signed.py',HERE.parent/'post_signed_mean_relaxation/relax.py']
    return dict(status='sharp-local-source-to-viscosity-witness-and-premise-obstruction',review=r,
      gates=dict(exact_controls=all(r['checks'].values()),three_frozen_records=len(bound)==3,historical_sources=all(checks.values())),
      dependency_evidence_sha256=bound,inherited_source_checks=checks,
      source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in paths},
      flags=dict(local_source_rate_D_and_viscosity_rate_nu_jointly_attained=True,
        inherited_size_premises_imply_positive_gain=False,actual_pulse_special_identity_excluded=False,
        actual_high_jets_support_complete_residual_and_blowup_certified=False))
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Sharpness witness, provenance or scope changed')
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args()
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print('initial viscosity sharpness witness:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'historical source checks')
if __name__=='__main__':main()
