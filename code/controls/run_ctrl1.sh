#!/bin/bash
# BAR-CTRL1 runner: one process at a time. Outputs go to ../../out/.
cd "$(dirname "$0")"
MC="python"
$MC mkcopy.py
CM_MUT=none     $MC cmcheck.py square 28 7 ../../out/ctrl1_out_pos_square_d28.json > ../../out/ctrl1_log_pos_square_d28.txt 2>&1
CM_MUT=none     $MC cmcheck.py square 12 7 ../../out/ctrl1_out_pos_square_d12.json > ../../out/ctrl1_log_pos_square_d12.txt 2>&1
CM_MUT=none     $MC cmcheck.py clw    28 5 ../../out/ctrl1_out_pos_clw_d28.json > ../../out/ctrl1_log_pos_clw_d28.txt 2>&1
CM_MUT=noinvrem $MC cmcheck.py square 28 7 ../../out/ctrl1_out_neg_inv_d28.json > ../../out/ctrl1_log_neg_inv_d28.txt 2>&1
CM_MUT=notail   $MC cmcheck.py square 12 7 ../../out/ctrl1_out_neg_tail_d12.json > ../../out/ctrl1_log_neg_tail_d12.txt 2>&1
echo ALLDONE
