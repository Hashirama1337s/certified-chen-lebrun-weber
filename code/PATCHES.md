# PATCHES: shipped scripts that differ from data/SHA256SUMS-PROD-met6.txt

Only absolute lab paths were replaced by relative ones (and the lab memory-cap wrapper by plain `python`). No other change.

## code/seatE/prod_A.py (lab: seatE/prod_A.py)
- lab sha256: `1d8d77b7862e2cf317f35e0415f048e74ea90658b9ee6ff908296f3f475fd1fa`
- shipped sha256: `e7cf2e9c9a49e5b0e189c19cfcb0a9e8ee0a16be3e80f100b795cded92964250`
- was: `np.load('<lab>/PAGE/vectors/v13_clw/seatW/bg_bach_M30.npz')`
- now: `np.load(os.path.join(common.HERE, '..', 'seatW', 'bg_bach_M30.npz'))`

## code/seatE/run_prod3.sh (lab: seatE/run_prod3.sh)
- lab sha256: `6cdd8a274dff28cdc1d2432ddbd12ad3f228dc2aeb9bb0add0fa7c5cb720fb0d`
- shipped sha256: `2a1af6b9a2e4c0d5c5e329b08be8d21477219f318adebf6b06c721c0036681b0`
- was: `MC="python <lab>/POCKETLAB/lab/tools/memcap.py 8 python"`
- now: `MC="python"`

## code/seatE/run_refine.sh (lab: seatE/run_refine.sh)
- lab sha256: `1e76dbc4baee2a5ce0fac0e98470376ea664c094f1295a28cc2c0664b61e5128`
- shipped sha256: `7a2daf689ddf23b8152955d4ffc613a2bffabcd10d05e0d3589af0c9bcb3a404`
- was: `MC="python <lab>/POCKETLAB/lab/tools/memcap.py 8 python"`
- now: `MC="python"`

## code/seatB2/assemble.py (lab: seatE/../seatB2/assemble.py)
- lab sha256: `616a852fd03631b3184df68ae0c8f5407bdf1a95f484fb2c556154fd4287c67b`
- shipped sha256: `4596f00d1c63c2ea2e7825d6ea36d05bce587800de1525300db62020f91b8697`
- was: `sys.path.insert(0, '<lab>/seatX')`
- now: `import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'seatX'))`

## controls/ (lab: release-work/ctrl1/)
The control scripts were adapted from the referee scratch versions before the BAR-CTRL1 run; their lab hashes are in data/ctrl1/CTRL1.json.
- cmcheck.py: lab `f580d2e3...`, shipped `2e1e5ff9...` (path lines only)
- mkcopy.py: lab `4b6caebc...`, shipped `c9cd3f85...` (path lines only)
- run_ctrl1.sh: lab `5ddc05c8...`, shipped `8e2d114b...` (path lines only)

Note: the docstring of make_cert.py still says lambda1 comes from cert_sampled.json; the code (BAR-MET5-amend1) reads the certified lower bound from G4_eps_m48.json.
