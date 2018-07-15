#!/bin/bash

python3 setup.py build
mkdir module
mv build/**/* ./module/
rm -rf build/
