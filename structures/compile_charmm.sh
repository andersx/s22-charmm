#!/bin/bash

cd /home/andersx/dev/charmm_cpe
export MAKE_COMMAND="make -j4"
./install.com gnu T gfortran +DFTBMKL > /dev/null 2> /dev/null
