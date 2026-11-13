"""amend3_triage.py -- list (default) or move aside (--move) the files written under older code that take a NEW amend3/amend4 path
and so must be recomputed (BAR-MET3-amend3 and -amend4, data policy):
  part B (clwB_e*): status other than 'ok'; chart 'xchart' (amend4: the chart switch now fires only on targeted crude bounds, so
                    every xchart file may take a different path; this includes amend3's 'xchart with infinite M');
                    status 'ok' with a non-finite or missing targeted value (amend4 status rule).
  part A (clw_e*):  B3 status other than 'ok' (e.g. 'FAIL torus': xchart retry also infinite -> amend3 splits).
python amend3_triage.py A|B DIR [--move]     moved files go to DIR_superseded_amend3/ (nothing is deleted)."""
import sys, os, json, glob, math, shutil

part, d = sys.argv[1], sys.argv[2]
move = '--move' in sys.argv
pat = 'clw_e*_*_*.json' if part == 'A' else 'clwB_e*_*_*.json'
dest = d.rstrip('/\\') + '_superseded_amend3'
NEED = ['K1', 'K2', 'kappa', 'Rc1', 'Rc2', 'Rc3', 'R0c1', 'R0c2', 'ric_lower']
hit = []
for fn in sorted(glob.glob(os.path.join(d, pat))):
    r = json.load(open(fn))
    rec = r['B3'] if part == 'A' else r
    st = rec.get('status')
    if st != 'ok':
        hit.append((fn, 'status %s' % st))
    elif part == 'B' and rec.get('chart') == 'xchart':
        hit.append((fn, 'xchart (amend4 trigger); infinite M %s' % [k for k, v in (rec.get('M') or {}).items() if not math.isfinite(float(v))]))
    elif part == 'B':
        ca = rec.get('cert_a') or {}
        bad = [k for k in NEED if ca.get(k) is None or not math.isfinite(float(ca[k]))]
        if bad:
            hit.append((fn, 'targeted value not finite %s' % bad))
for fn, why in hit:
    print('%s: %s%s' % (os.path.basename(fn), why, ' -> moved' if move else ''))
    if move:
        os.makedirs(dest, exist_ok=True)
        shutil.move(fn, os.path.join(dest, os.path.basename(fn)))
print('%d file(s) %s' % (len(hit), 'moved to ' + dest if move else 'to recompute (dry run)'))
