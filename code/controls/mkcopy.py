# BAR-CTRL1: makes a mutable COPY of the pinned seatE/cm.py in controls/mut/ (pinned file is only read, never written).
# The two mutation edits are the ones used by the referee (REFEREE-CM); the copy is written to controls/mut/.
import hashlib, os
SE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'seatE')
PIN = '8697432728c08f732418af0f17b4a3cc123096cbd796c2f9dad989c040542025'
raw = open(os.path.join(SE, 'cm.py'), 'rb').read()
assert hashlib.sha256(raw).hexdigest() == PIN, 'pinned cm.py hash mismatch'
src = open(os.path.join(SE, 'cm.py')).read()   # text mode, as in the referee's script
a = "                tail += x.abs_upper()"
assert src.count(a) == 3
src = src.replace(a, "                tail += (x.abs_upper() if MUT != 'notail' else 0)")
b = "        return CM(q.p, (q.norm() * delta / (1 - delta)).abs_upper(), deg=q.deg)"
assert src.count(b) == 1
src = src.replace(b, "        return CM(q.p, arb(0) if MUT == 'noinvrem' else (q.norm() * delta / (1 - delta)).abs_upper(), deg=q.deg)")
src = src.replace("D = 28\n", "import os\nMUT = os.environ.get('CM_MUT', 'none')\nD = 28\n", 1)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mut', 'cm.py')
open(out, 'w', newline='').write(src)
print('copy written', out, hashlib.sha256(src.encode()).hexdigest())
