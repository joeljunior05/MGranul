#!/bin/bash

if [ -d "module" ]; then
  rm -rf module
fi

python3 setup.py build
mkdir module
mv build/**/* ./module/
rm -rf build/
