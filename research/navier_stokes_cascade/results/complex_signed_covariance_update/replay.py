#!/usr/bin/env python3
"""Compatibility replay: only three near-zero manufactured roundoff fields."""
import argparse
from functools import lru_cache
import importlib.util
import math
from pathlib import Path
import runpy
import sys

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/'first_sourced_wave_update'
ROUND_FIELDS=frozenset('full_viscous_potential_pressure_cutoff.'+name for name in
    ('maximum_relative_commutator_error','maximum_curl_product_error','maximum_divergence_error'))
ABS_ROUNDOFF=2e-12  # Stricter than every unchanged underlying acceptance gate.


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def compatible_compare(a,b,path,strict):
    if type(a) is not type(b): raise ValueError('Changed type at '+path)
    if isinstance(a,dict):
        if set(a)!=set(b): raise ValueError('Changed keys at '+path)
        for key in a: compatible_compare(a[key],b[key],path+'.'+key,strict)
    elif isinstance(a,list):
        if len(a)!=len(b): raise ValueError('Changed list at '+path)
        for i,(x,y) in enumerate(zip(a,b)): compatible_compare(x,y,path+'.'+str(i),strict)
    elif isinstance(a,float) and any(path==f or path=='root.reviews.'+f for f in ROUND_FIELDS):
        if not (math.isfinite(a) and math.isfinite(b) and 0<=a<=ABS_ROUNDOFF and 0<=b<=ABS_ROUNDOFF):
            raise ValueError('Roundoff certificate exceeded at '+path)
    else: strict(a,b,path)


@lru_cache(None)
def parent_module():
    sys.path.insert(0,str(PARENT))
    try: parent=load('complex_signed_frozen_first_source',PARENT/'audit.py')
    finally: sys.path.pop(0)
    strict=parent.compare_records
    parent.compare_records=lambda a,b,path='root': compatible_compare(a,b,path,strict)
    return parent


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tests',action='store_true')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify-record',type=Path,default=PARENT/'evidence.json')
    args=parser.parse_args(argv);parent=parent_module()
    if args.tests:
        if args.output: parser.error('--tests does not accept --output')
        sys.modules['audit']=parent;sys.argv=[str(PARENT/'test_audit.py'),'-v']
        runpy.run_path(str(PARENT/'test_audit.py'),run_name='__main__')
    else:
        if args.output is None: parser.error('--output is required for evidence replay')
        parent.main(['--verify-record',str(args.verify_record),'--output',str(args.output)])


if __name__=='__main__': main()
