"""Rebuild RMUT ESS analysis-set views from saved formal summaries.

This only applies an algebraic transformation to saved summaries. It does not
generate patient records, draw model parameters, or fit a model.
"""
from pathlib import Path
import sys
_SUPPORT_CODE_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(_SUPPORT_CODE_ROOT / 'src/common'))
from output_paths import require_external_output
import csv
import argparse
import math

NMAX={'IMpower131':680,'KEYNOTE-598':568,'OAK':850,'survival_colon':620}
RHO_LEVELS={
    'IMpower131':{'rho_0p3':'low','rho_0p7':'high'},
    'KEYNOTE-598':{'rho_0p3':'low','rho_0p7':'high'},
    'OAK':{'rho_0p5':'low','rho_0p9':'high'},
    'survival_colon':{'rho_0p7':'low','rho_0p98':'high'},
}


def summarize(data_directory):
    root=Path(data_directory)
    src=root/'RMUT_Run_Summaries.csv'
    result=[]
    with src.open(encoding='utf-8-sig',newline='') as f:
        for r in csv.DictReader(f):
            if r['stage']=='calibration':
                continue  # Calibration deliberately retains its all-attempts rule.
            case={'IMpower':'IMpower131','KEYNOTE598':'KEYNOTE-598'}.get(r['case'],r['case'])
            scenario={'H1':'H1_standard','H1_calibration':'H1_standard','Realcase':'H1_observed'}.get(r['scenario'],r['scenario'])
            rho_level=RHO_LEVELS[case][r['rho_level']] if r['stage']=='rho' else ''
            total,valid,failed=(int(float(r[k])) for k in ['total','valid','failures'])
            assert total==valid+failed and 0<valid<=total
            all_ess=float(r['ESS'])
            # The formal summary assigns Nmax to each failed replicate.
            valid_ess=(all_ess*total-failed*NMAX[case])/valid
            assert math.isfinite(valid_ess) and 0<valid_ess<=NMAX[case]
            result.append({'analysis':r['stage'],'case':case,'scenario':scenario,'rho_level':rho_level,
                'n_total':total,'n_valid':valid,'n_failed':failed,'Nmax':NMAX[case],
                'ESS_all_attempts':all_ess,'ESS_valid_only':valid_ess,
                'ESS_fraction_valid_only':valid_ess/NMAX[case]})
    assert len(result)==48
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_directory',type=Path)
    parser.add_argument('output_csv',type=Path)
    args=parser.parse_args()
    args.output_csv = require_external_output(args.output_csv, _SUPPORT_CODE_ROOT)
    if args.output_csv.exists():
        parser.error('Choose a new output file to preserve existing records.')
    rows=summarize(args.data_directory)
    with args.output_csv.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    print(f'Wrote {len(rows)} saved-summary views to {args.output_csv.name}')
