#!/bin/bash
mkdir dist

tar xzf openbabel-2.4.1.tar.gz
cd openbabel-2.4.1/
mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=../../dist -DWITH_INCHI=OFF
make -j4
make install
cd ../..

tar xzf qhull-2015-src-7.2.0.tgz
cd qhull-2015.2/
mkdir cbuild
cd cbuild
cmake .. -DCMAKE_INSTALL_PREFIX=../../dist
make -j4
make install

