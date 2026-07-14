#!/usr/bin/env python3
"""Verify archived benchmark evidence or reproduce the portable 4D suite."""
from __future__ import annotations
import argparse,json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/'results'
REQUIRED=[
 RESULTS/'benchmark_summary.json',RESULTS/'adaptive_fem/results.json',RESULTS/'stability_bridge/results.json',
 RESULTS/'axisymmetric_2d/results.json',RESULTS/'triadic_3d/results.json',
 RESULTS/'observer_manifold/conditioned_results.json',RESULTS/'observer_manifold/manifold_results.json',
 RESULTS/'evolution6_4d/results.json']
def verify():
 bad=[]
 for p in REQUIRED:
  try:
   if not json.loads(p.read_text(encoding='utf-8')):bad.append(str(p))
  except Exception as e:bad.append(f'{p}: {e}')
 return bad
def reproduce_4d():
 out=RESULTS/'evolution6_4d'/'reproduced';out.mkdir(parents=True,exist_ok=True)
 env=os.environ.copy();env['TRIAD_BENCH_OUT']=str(out)
 return subprocess.run([sys.executable,str(ROOT/'scripts/evolution6_4d.py')],env=env,check=False).returncode
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true');ap.add_argument('--reproduce-4d',action='store_true');a=ap.parse_args()
 bad=verify();print(f'Archived benchmark evidence: {len(REQUIRED)-len(bad)}/{len(REQUIRED)} valid')
 if bad:
  print('Invalid evidence:');[print(' -',x) for x in bad];return 2
 if a.reproduce_4d or not a.quick:
  rc=reproduce_4d()
  if rc:return rc
 print('TRIAD benchmark verification: PASS');return 0
if __name__=='__main__':raise SystemExit(main())
