#!/bin/bash

path=$1"/pop/"

for i in $(ls $path); do grep -v ":-:" $path/$i > $1/$i; done

rm -r $path
