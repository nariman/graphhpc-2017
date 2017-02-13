#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo 'Test generation script'
    echo 'Please, provide:'
    echo '  - a graph type (`random` or `rmat`) as an first argument'
    echo '  - a `scale` number as an second argument'
    echo
    echo 'The rest of arguments will be passed to the answer generator executable'
    exit 1
fi

cd $(dirname $0)
cd ./../

cd ./tests/
./../bin/gen_$1 -s $2
./../bin/gen_valid_info -in $1-$2 ${@:3}
