#!/bin/bash

# download libs
INSTALL_DIR="install_build"
mkdir ${INSTALL_DIR}
cd ${INSTALL_DIR}

# nauty
wget http://pallini.di.uniroma1.it/nauty27rc1.tar.gz
tar -xzf nauty27rc1.tar.gz
mv nauty27rc1 nauty
cd nauty
./configure
sed -i '/CFLAGS=/s/$/ -fPIC/' makefile
make
sed -i '/#define _FILE_OFFSET_BITS 0/c\#ifndef _FILE_OFFSET_BITS\n#define _FILE_OFFSET_BITS 0' nauty.h
sed -i '/#undef _FILE_OFFSET_BITS/c\#undef _FILE_OFFSET_BITS\n#endif' nauty.h
cd ..

# pybind11
git clone https://github.com/pybind/pybind11 pybind11

# aboria
git clone https://github.com/martinjrobins/Aboria Aboria
cd Aboria
sed -i 's/\<WORDSIZE\>/WORDSIZE_NANNOFLANN/g' third-party/nanoflann/nanoflann.hpp
cd ..

# mover para libs
SRC_LIBS=../src/cpp/libs
mv nauty Aboria pybind11 ${SRC_LIBS}

# excluir pasta
cd ..
rm -rf ${INSTALL_DIR}
