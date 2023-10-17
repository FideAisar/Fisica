#!/bin/bash
#
# Shows only the subdirectories in the current dir
# It doesn't show symbolic links
for dir in *;
do
	if [ -d "$dir" ] && [ ! -L "$dir" ] 
	then
		echo "$dir"
	fi
done
