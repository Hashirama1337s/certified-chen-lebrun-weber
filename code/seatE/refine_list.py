"""refine_list.py -- the sub-squares to refine (BAR-MET5): every SRC_DIR/clw_e*_*_*.json whose A.integrals.lap2.upper is more than
THRESH (default 0.005 = 0.5%) of the total over all files.  Prints the shares (stderr) and the k:i:j list (stdout, one line).
python refine_list.py SRC_DIR [THRESH]"""
import sys, os, json, glob

d = sys.argv[1]; th = float(sys.argv[2]) if len(sys.argv) > 2 else 0.005
R = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(d, 'clw_e*_*_*.json')))]
tot = sum(r['A']['integrals']['lap2']['upper'] for r in R)
sel = sorted((r['A']['integrals']['lap2']['upper'] / tot, r['k'], r['i'], r['j']) for r in R)[::-1]
for sh, k, i, j in sel:
    if sh > th:
        print('e%d [%d,%d] lap2 share %.2f%%' % (k, i, j, 100 * sh), file=sys.stderr)
print(','.join('%d:%d:%d' % (k, i, j) for sh, k, i, j in sel if sh > th))
