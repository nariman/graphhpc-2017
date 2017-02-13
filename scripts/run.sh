#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo 'Solution runner script'
    echo 'Please, provide:'
    echo '  - a graph type (`random` or `rmat`) as an first argument'
    echo '  - a `scale` number as an second argument'
    echo
    echo 'The rest of arguments will be passed to the solution executable'
    exit 1
fi

cd $(dirname $0)
cd ./../

./solution -in ./tests/$1-$2 ${@:3}

./bin/validation -ans ./tests/$1-$2.ans -res ./tests/$1-$2.res
rm ./tests/$1-$2.res
