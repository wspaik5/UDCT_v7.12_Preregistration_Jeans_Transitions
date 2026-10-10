#!/usr/bin/env python3
"""UDCT v7.12 prospective follow-up. Default: plan; --selftest: synthetic only.
Galaxy calculations require --run and an operator-attested public registration.
See the accompanying PDF for the complete interpretation and decision rules.
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parent
FROZEN = {
    'UDCT_v7_11_Preregistration_Anisotropic_Jeans.py': '3f00abae251d7be64a4ca5abb48b4f550644b21c946206571fd448464b13ade2',
    'UDCT_v7_10_Postregistration_run.py': 'a4d1a07ebcaa06348be3022d0f1f6ca463347e00194f7fde59f2bd9dbcd521b7',
    'UDCT_v7_06_Hessian_Term_Sign_Audit_reproduce.py': 'a64440e1350ae8bf10eb3188ec34829b6d586abf275de5cdc841a50ad3478284',
}
BETAS = tuple(round(i * .05, 2) for i in range(11))
ELLS = (1., 1.2, 1.5, 2., 2.5, 3.)
LAWS = ('BASE', 'A-f1/2', 'A-f2', 'A-Qs', 'A-Qd')
CONTROLS = ('QUMOND-nu0', 'QUMOND-nu1', 'QUMOND-nu2')
GRID_INNER = tuple((g, inn) for inn in (.002, .0005) for g in ('smoke', 'default', 'refined'))
DECIDING = (('default', .002), ('refined', .002), ('default', .0005), ('refined', .0005))
PLAN = dict(version='v7.12', author='Won Shik Paik', preparation_date='2026-10-10',
    status='prospective follow-up after known v7.11 results; not publicly registered yet',
    target='Hydrus I', profile='unsoftened Tapered cusp', betas=BETAS,
    screening_lengths_pc=ELLS, laws=LAWS, controls=CONTROLS,
    outer_pressure_si=0., outer_radius_rh=40., common_radii_rh=[.002, 10.], common_points=1200,
    fixed_field_evaluations=198, fixed_pressure_records=2178, fixed_primary_decisions=363,
    stage_1='ell=1 pc, all five screened laws and three QUMOND controls, all 11 betas',
    stage_2='remaining five lengths, all five screened laws, all 11 betas; controls not repeated',
    adaptive='all adjacent raw sign-changing brackets; beta width <=0.001 or ell width <=0.01 pc; max 20 steps; stop unresolved on numerical margin',
    adaptive_max_brackets=575, adaptive_max_midpoints=11500,
    scope='sampled signed Jeans second moments, not equilibrium, stability or observational fit',
    deferred='core softening and local-u analysis require separate pre-execution specification',
    known_records={'v7.11 preregistration':'https://doi.org/10.5281/zenodo.23244522',
                   'v7.11 results':'https://doi.org/10.5281/zenodo.23263099'})

def now():
    return datetime.now(timezone.utc).isoformat()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_frozen():
    for name, expected in FROZEN.items():
        if not (ROOT/name).is_file() or sha(ROOT/name) != expected:
            raise RuntimeError('Missing or modified frozen source: '+name)
    spec = importlib.util.spec_from_file_location('udct712_frozen711', ROOT/'UDCT_v7_11_Preregistration_Anisotropic_Jeans.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module

def native_min(rec):
    return rec['min_sigma2_kms2'] * 1e6 / rec['vref2_si']

def classify(old, records):
    chosen = [records[f'{g}|inner={inn:g}'] for g, inn in DECIDING]
    decision = old.decide([r['common_sigma2_over_vref2'] for r in chosen])
    if 'tolerance_over_vref2' not in decision:
        decision.update(F=None, raw_sign=None, numerical_sign=None)
        return decision
    minima = [min(float(np.min(r['common_sigma2_over_vref2'])), native_min(r)) for r in chosen]
    F = min(minima)
    tau = decision['tolerance_over_vref2']
    if decision['status'] == 'P' and min(minima) <= tau:
        decision.update(status='I', reason='native minimum fails positive margin')
    decision.update(F=F, raw_sign=(-1 if F < 0 else 1 if F > 0 else 0),
        numerical_sign=(-1 if F < -tau else 1 if F > tau else 0),
        physical_candidate_eligibility=True,
        eligibility_note='cusp beta<=0.5 necessary screen only; no positive distribution function established')
    return decision

def bisect_bracket(evaluate, lo, hi, width, max_steps=20):
    """Refine a raw zero bracket; stop at uncertainty, never bisect N/P labels.
    F is the minimum normalized variance across four common/native comparisons.
    No monotonicity or uniqueness is assumed. Return endpoint metadata.
    """
    left, right = evaluate(lo), evaluate(hi)
    if left['F'] is None or right['F'] is None or left['F'] * right['F'] >= 0:
        return dict(status='no_strict_sign_bracket', lo=lo, hi=hi, steps=0)
    history = []
    if not left['numerical_sign'] or not right['numerical_sign']:
        return dict(status='unresolved_endpoint_margin', lo=lo, hi=hi, left=left, right=right, steps=0)
    status = 'iteration_limit'
    for _ in range(max_steps):
        if hi-lo <= width:
            status = 'coordinate_width_reached'; break
        mid = (lo+hi)/2
        middle = evaluate(mid)
        history.append(dict(x=mid, decision=middle))
        if middle['F'] is None or not middle['numerical_sign']:
            status = 'unresolved_midpoint_margin'; break
        if left['F'] * middle['F'] < 0:
            hi, right = mid, middle
        else:
            lo, left = mid, middle
    if hi-lo <= width:
        status = 'coordinate_width_reached'
    return dict(status=status, lo=lo, hi=hi, left=left, right=right, steps=len(history), history=history,
        interpretation='numerical sampled-zero bracket only; may contain multiple roots')

def selftest():
    old = load_frozen()
    old.selftest()  # analytic/synthetic arrays and frozen imports only
    r = np.geomspace(.1, 3., 20000)
    for beta in BETAS:
        pressure, _, q = old.solve(r, np.ones_like(r), np.full_like(r, 2.), beta, .7)
        exponent = 2*beta+1
        exact = (2*(r[-1]**exponent-r**exponent)/exponent/r[-1]**(2*beta)+.7)/q
        assert np.max(abs(pressure-exact))/np.max(abs(exact)) < 2e-6
    def analytic(x):
        value = x-.317
        return dict(F=value, numerical_sign=-1 if value<0 else 1 if value>0 else 0)
    result = bisect_bracket(analytic, 0., .5, .001)
    assert result['lo'] <= .317 <= result['hi'] and result['hi']-result['lo'] <= .001
    def uncertain(x):
        d = analytic(x)
        if abs(d['F']) < .1: d['numerical_sign'] = 0
        return d
    assert bisect_bracket(uncertain, 0., .5, .001)['status'] == 'unresolved_midpoint_margin'
    assert bisect_bracket(analytic, .4, .5, .001)['status'] == 'no_strict_sign_bracket'
    records = {f'{g}|inner={inn:g}':dict(common_sigma2_over_vref2=[1.]*1200,
        min_sigma2_kms2=-1e-6, vref2_si=1.) for g, inn in DECIDING}
    assert classify(old, records)['status'] == 'I'
    records = {k:dict(common_sigma2_over_vref2=[0.]*1200, min_sigma2_kms2=0., vref2_si=1.) for k in records}
    assert classify(old, records)['status'] == 'I'
    # Compact compressed output round trip, without a galaxy or persistent file.
    encoded = gzip.compress(json.dumps(PLAN, allow_nan=False).encode())
    assert json.loads(gzip.decompress(encoded))['fixed_primary_decisions'] == 363
    print('v7.12 SELFTEST PASSED: synthetic only; no galaxy cases evaluated.')

def execute(args, journal):
    old = load_frozen()
    base = old.load_baseline()
    obj = base.DWARFS['Hydrus I']
    grids = {g:base.make_grid(g) for g in ('smoke','default','refined')}
    environments = {g:base.V706.galactic_field(obj, grid)[:2] for g,grid in grids.items()}
    law_map = {law.name:law for law in (base.BASE_LAW,)+base.PART_A_LAWS}
    for kind in ('nu0','nu1','nu2'):
        law_map['QUMOND-'+kind] = base.Law('QUMOND-'+kind, 1., kind)
    field_cache, evaluated = {}, {}
    all_gates = []
    counts = dict(field_evaluations=0, pressure_records=0, case_decisions=0)
    def event(kind, **data):
        journal.write(json.dumps(dict(event=kind, utc=now(), **data), separators=(',',':'), allow_nan=False)+'\n')
        journal.flush()
    def field(name, ell, g, inn):
        key = (name,ell,g,inn)
        if key not in field_cache:
            grid = grids[g]
            e,h = environments[g]
            st = base.model_state(obj, 'Tapered cusp', base.mesh('Tapered cusp', grid, rmin=inn), e,h,grid,law_map[name],ell,0.)
            work, inward = base.V706.work_from_state(st)
            _, old_sigma = base.V706.aperture_jeans(st, inward, obj.re_pc)
            p,_,_ = old.solve(st['r'],st['rho'],inward,0.)
            sigma = np.divide(p,st['rho'],out=np.zeros_like(p),where=st['rho']>1e-280)
            vref2 = float(np.interp(1.,st['r']/st['rh'],st['gn']))*st['rh']
            error = float(np.max(abs(sigma-old_sigma))/vref2)
            passed = bool(np.isfinite(error) and error<=1e-10)
            all_gates.append(passed)
            counts['field_evaluations'] += 1
            event('field_gate', law=name, ell_pc=ell, grid=g, inner_rh=inn, error=error, passed=passed)
            if not passed: raise RuntimeError('Gate R failed; no outcome claims permitted')
            field_cache[key] = (st,inward,work)
        return field_cache[key]
    def evaluate(name,ell,beta,stage):
        key = f'{name}|ell={ell}|beta={beta:.12g}'
        if key in evaluated: return evaluated[key]['decision']
        records = {}
        for g,inn in GRID_INNER:
            st,inward,work = field(name,ell,g,inn)
            rec = old.summary(st['r'],st['rho'],inward,beta,st['rh'],st['gn'])
            rec.pop('common_sigma2_kms2')  # redundant: reconstruct from normalized array and vref2
            rec['v710_work_ratio'] = work['virial']
            records[f'{g}|inner={inn:g}'] = rec
            counts['pressure_records'] += 1
        decision = classify(old,records)
        evaluated[key] = dict(stage=stage, law=name, ell_pc=ell, beta=beta, decision=decision, records=records)
        counts['case_decisions'] += 1
        event('case_complete', key=key, **evaluated[key])
        print(f'{stage}: {key} {decision["status"]}',file=sys.stderr)
        return decision
    stage1 = [(name,1.) for name in LAWS] + [(name,None) for name in CONTROLS]
    stage2 = [(name,ell) for ell in ELLS[1:] for name in LAWS]
    for stage,specs in (('stage1',stage1),('stage2',stage2)):
        for name,ell in specs:
            for beta in BETAS: evaluate(name,ell,beta,stage)
    fixed_counts = counts.copy()
    assert fixed_counts == dict(field_evaluations=198,pressure_records=2178,case_decisions=363)
    brackets = []
    # Enumerate from the entire fixed matrix BEFORE adaptation. Each original
    # adjacent sign bracket is refined once; no outcome-driven extra search.
    jobs = []
    for name in LAWS:
        for ell in ELLS:
            for lo,hi in zip(BETAS[:-1],BETAS[1:]):
                jobs.append(('beta',name,ell,lo,hi,.001))
        for beta in BETAS:
            for lo,hi in zip(ELLS[:-1],ELLS[1:]):
                jobs.append(('ell_pc',name,beta,lo,hi,.01))
    assert len(jobs) == 575
    for axis,name,fixed,lo,hi,width in jobs:
        evaluator = (lambda x,n=name,f=fixed:evaluate(n,f,x,'adaptive_beta')) if axis=='beta' else (lambda x,n=name,f=fixed:evaluate(n,x,f,'adaptive_ell'))
        left,right = evaluator(lo),evaluator(hi)
        if left['F'] is None or right['F'] is None or left['F']*right['F']>=0: continue
        result = bisect_bracket(evaluator,lo,hi,width)
        bracket = dict(axis=axis,law=name,fixed_coordinate=fixed,original_lo=lo,original_hi=hi,**result)
        brackets.append(bracket)
        event('bracket_complete',**bracket)
    return dict(protocol=PLAN, completed_utc=now(),fixed_counts=fixed_counts,total_counts=counts,
        isotropic_gate_passed=all(all_gates),common_r_over_rh=old.COMMON_X.tolist(),
        cases=evaluated,brackets=brackets,
        interpretation='No full equilibrium, universal-u law, observational fit or MOND exclusion established.')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--selftest',action='store_true')
    mode.add_argument('--run',action='store_true')
    parser.add_argument('--registration-url')
    parser.add_argument('--registration-published-utc')
    parser.add_argument('--confirm-public-registration',action='store_true')
    parser.add_argument('--output',type=Path,help='new filename ending .json.gz')
    args = parser.parse_args()
    if args.selftest: selftest(); return
    if not args.run:
        print(json.dumps(PLAN,indent=2)); return
    parsed = urlparse(args.registration_url or '')
    if parsed.scheme!='https' or not parsed.netloc or not args.confirm_public_registration:
        parser.error('public https registration URL and --confirm-public-registration are required')
    try:
        published = datetime.fromisoformat((args.registration_published_utc or '').replace('Z','+00:00'))
        if published.utcoffset() is None or published.utcoffset().total_seconds()!=0 or published>=datetime.now(timezone.utc): raise ValueError()
    except ValueError:
        parser.error('--registration-published-utc must be a past timezone-aware UTC timestamp')
    if not args.output or not str(args.output).endswith('.json.gz'):
        parser.error('--output must end .json.gz')
    partial = Path(str(args.output)+'.partial')
    journal_path = Path(str(args.output)+'.journal.jsonl.gz')
    if any(p.exists() for p in (args.output,partial,journal_path)):
        parser.error('output or attempt files exist; choose a new output to preserve all attempts')
    selftest()
    provenance = dict(started_utc=now(),registration_url=args.registration_url,
        registration_published_utc=published.isoformat(),registration_publication_attested_by_operator=True,
        registration_url_independently_verified=False,runner_sha256=sha(__file__),frozen_sources=FROZEN,
        command=sys.argv,environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()))
    # Reserve attempt files exclusively; preserve partial data and journal on failure.
    with journal_path.open('xb') as raw:
        with gzip.open(raw,'wt',encoding='utf-8') as journal:
            journal.write(json.dumps(dict(event='start',provenance=provenance))+'\n');journal.flush()
            try:
                result = execute(args,journal)
                result['provenance'] = provenance
                with partial.open('xb') as output_raw:
                    with gzip.open(output_raw,'wt',encoding='utf-8') as out:
                        json.dump(result,out,separators=(',',':'),allow_nan=False)
                with gzip.open(partial,'rt',encoding='utf-8') as check:
                    verified = json.load(check)
                if verified['fixed_counts']['case_decisions']!=363: raise RuntimeError('output completeness check failed')
                os.link(partial,args.output)  # exclusive final installation; never overwrite
                partial.unlink()
                journal.write(json.dumps(dict(event='success',utc=now(),output_sha256=sha(args.output)))+'\n')
            except Exception as exc:
                journal.write(json.dumps(dict(event='failed',utc=now(),error=repr(exc),claims_permitted=False))+'\n')
                raise
    print('Complete compressed result saved to '+str(args.output))

if __name__=='__main__':
    main()
