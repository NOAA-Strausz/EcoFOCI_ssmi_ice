#!/bin/bash

for file in $*
do
    echo "Working on file $file"
    /home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/plot_ssmi.py -ex chukchi $file
done
    
