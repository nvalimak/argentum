#!/bin/bash

set -o nounset
set -o errexit
set -o pipefail

i=1
while true
do
    ./rand-matrix $i 10 10 > rand.input
    ./main -S -i rand.input -v --debug --no-prediction --enumerate | ./enumerate-example
    i=$((i+1))
done
