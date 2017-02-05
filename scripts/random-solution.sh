#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo 'Please, provide a SCALE number'
    exit 1
fi

cd $(dirname $0)
cd ../bin/

./gen_random -s $1
./solution -in random-$1
./validation -in random-$1 -res random-$1.res

rm random-$1
rm random-$1.res
