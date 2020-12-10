#!/bin/bash

# For Linux
rm lib*
sed -i "s,local_dir = '.*',local_dir = \'`pwd`\'," igmabs.f90
echo "initpackage igmabs lmodel_igmabs.dat "`pwd`"\nquit\ny" | xspec


rm *~ *.o
rm *FunctionMap.* lpack_* 
rm -f *.mod Makefile
