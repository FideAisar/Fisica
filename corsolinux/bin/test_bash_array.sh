#!/bin/bash
declare -a arrayname=("element1" "element2" "element3")
for el in  "${arrayname[@]}"
do
    echo "item: $el"
done
       
