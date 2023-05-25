#! /bin/sh
#
# Invoke with
# runtargets.sh path/to/targets

./laplace_problem &
PID=$!

./boundary_condition
RET=$?

set -e
wait $PID
exit $RET
