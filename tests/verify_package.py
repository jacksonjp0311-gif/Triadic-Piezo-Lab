#!/usr/bin/env python3
from pathlib import Path
import json,re,sys
ROOT=Path(__file__).resolve().parents[1];errors=[]
for rel in ['OPEN_LAB.html','README.md','manifest.json','benchmarks/run_benchmarks.py','benchmarks/results/benchmark_summary.json','source/ORIGINAL_CLASSIC_STYLE_REFERENCE.html']:
 if not (ROOT/rel).exists():errors.append(f'Missing {rel}')
html=(ROOT/'OPEN_LAB.html').read_text(encoding='utf-8')
for token in ['displayPanel','uiBrightness','sceneBloom','evidenceOverlay','axis2dOverlay','runFem','runCondense','runBridge','three@0.185.1']:
 if token not in html:errors.append(f'HTML token missing: {token}')
ids=re.findall(r'\bid="([^"]+)"',html);dups=sorted({x for x in ids if ids.count(x)>1})
if dups:errors.append('Duplicate HTML IDs: '+', '.join(dups[:20]))
for p in (ROOT/'benchmarks/results').rglob('*.json'):
 try:json.loads(p.read_text(encoding='utf-8'))
 except Exception as e:errors.append(f'Invalid JSON {p.relative_to(ROOT)}: {e}')
print('Package root:',ROOT);print('HTML bytes:',(ROOT/'OPEN_LAB.html').stat().st_size);print('Unique IDs:',len(set(ids)))
if errors:
 print('VERIFY: FAIL');[print(' -',e) for e in errors];raise SystemExit(1)
print('VERIFY: PASS')
