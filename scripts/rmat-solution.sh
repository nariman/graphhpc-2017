#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo 'Please, provide a SCALE number'
    exit 1
fi

cd $(dirname $0)
cd ../bin/

./gen_rmat -s $1
./solution -in rmat-$1
./validation -in rmat-$1 -res rmat-$1.res

rm rmat-$1
rm rmat-$1.res
