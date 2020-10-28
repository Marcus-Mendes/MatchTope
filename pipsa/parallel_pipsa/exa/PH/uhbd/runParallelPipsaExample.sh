#!/bin/bash
../../../bin/npotsim -pg ../../../bin/2potsim_skin -fn names -lg sims.log -pr 3 -sk 4
mv sims.log sims.log.conventional

../../../bin/npotsim -pg ../../../bin/2potsim_skin -fn names -lg sims.log -pr 3 -sk 4 -ma 4 && make -j 2 -k
