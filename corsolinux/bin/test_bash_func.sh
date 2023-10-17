#!/usr/bin/env bash

# Use on the command line as
# bash test_bash_func.sh my_f1 pippo
# or
# bash test_bash_func.sh my_f2 pippo


my_f1() {
    echo "$1"
    echo "number of arguments $#"
}

#my_f1 pluto

my_f2() {
    echo "$1 $1"
}


"$@"
