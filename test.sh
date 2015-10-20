#!/bin/bash

set -o nounset
set -o errexit
set -o pipefail

while true
do
    ./rand-matrix 10 10 > rand.input
    ./main -S -i rand.input -v --debug --rewind
done
