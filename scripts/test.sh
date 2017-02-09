#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo 'Please, provide:'
    echo '  - a graph type (random OR rmat) as an first argument'
    echo '  - a SCALE number as an second argument'
    echo
    echo 'The rest of arguments will be passed to the solution executable'
    exit 1
fi

cd $(dirname $0)
cd ../

cd ./bin/
./gen_$1 -s $2

cd ./../
./solution -in ./bin/$1-$2 ${@:3}

cd ./bin
./validation -in $1-$2 -res $1-$2.res
rm $1-$2
rm $1-$2.res
