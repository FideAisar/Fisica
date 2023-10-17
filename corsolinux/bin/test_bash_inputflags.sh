#!/bin/bash

# check reading input using flags
# case could be more useful
if [ $1 = -a ]
then
    var1=$2
elif [ $1 = -b ]
then
    var2=$2
fi

if [ $3 = -a ]
then
    var1=$4
elif [ $3 = -b ]
then
    var2=$4
fi

echo $var1 $var2
