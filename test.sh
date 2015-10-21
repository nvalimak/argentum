#!/bin/bash

set -o nounset
set -o errexit
set -o pipefail

while true
do
    ./rand-matrix 20 20 > rand.input
    ./main -S -i rand.input -v --debug 1 --rewind
done
