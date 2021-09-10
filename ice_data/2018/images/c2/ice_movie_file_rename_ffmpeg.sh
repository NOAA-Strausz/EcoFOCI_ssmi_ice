#!/bin/bash

count=1
for file in $*
do
    mv $file ice_plot_${count}.png
    ((count++))
done
