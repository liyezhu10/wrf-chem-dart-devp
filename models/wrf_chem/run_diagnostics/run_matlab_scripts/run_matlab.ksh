#!/bin/ksh -x
#
# Script to run matlab
#bsub -Is -q geyser -W 1:00 -n 1 -P P19010000 matlab
bsub -Is -q caldera -W 1:00 -n 1 -P P19010000 matlab
exit
