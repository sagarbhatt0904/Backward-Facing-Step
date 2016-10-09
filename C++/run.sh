#!/bin/sh

echo 'Running Navier-stokes Solver for Backward-Facing Step problem'
echo '\n'
make -C source/

source/./bkwd

rm source/bkwd