#!/bin/bash

cpwd=$(pwd -P)

cd ~
source set664.sh
cd ${cpwd}

cd cmt
cmt br cmt config
cmt br make clean
cmt br make
source setup.sh

cd ${cpwd}
