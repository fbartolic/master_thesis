#!/bin/bash
echo "Compiling binary_c..."
cd ../src
make clean
make
make libbinary_c.so
echo "\nLinking library to /usr/lib ..."
sudo ln -sf ~/Drive/Work/master_thesis/binary_c/src/libbinary_c.so /usr/lib/libbinary_c.so
echo "\nCompiling libevolve.so\n"
cd ../python_API
rm libevolve.so
make
echo "...done"
