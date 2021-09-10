#!/bin/bash

for file in $*
do
    /home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/plot_ssmi.py -m ck2 -ex chukchi $file
done
    
