#!/bin/bash

path=$1"_filt"
mkdir $path

for i in $(ls $1); do grep -v ":-:" $1/$i > $path/$i; done

mv $1 $1_padded
mv $path $1

