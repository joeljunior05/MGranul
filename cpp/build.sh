#!/bin/bash

FOLDER_NAME=release

export ROOT_FOLDER=`pwd`

export LDFLAGS="-L$ROOT_FOLDER/local/lib"
export CPPFLAGS="-I$ROOT_FOLDER/local/include"
export PKG_CONFIG_PATH="$ROOT_FOLDER/local/lib/pkgconfig"
export OpenCV_DIR="$ROOT_FOLDER/local/share/OpenCV"
export PATH="$ROOT_FOLDER/local/bin:$PATH"

########### BUILD SCRIPT FOR MSGRANUL ###############

if [ ! -d "$FOLDER_NAME" ]; then
  mkdir $FOLDER_NAME

  cd $FOLDER_NAME
  cmake .. -DWITH_OMP=OFF
  make

else
  unset get_op
  echo "$FOLDER_NAME folder exists! Would you like erase it (y/n)?"
  read get_op

  if [[ ${get_op} = 'y' || ${get_op} = 'Y'  ]]; then
    rm -r $FOLDER_NAME
    mkdir $FOLDER_NAME
    cd $FOLDER_NAME
    cmake .. -DWITH_OMP=OFF
    make
  else
    echo "$FOLDER_NAME not erased."
  fi
fi

cd $ROOT_FOLDER

if [ -d "bin" ]; then
  rm -rf bin
fi

mkdir bin
cp release/exe_granul ./bin/
cp release/libgranul.a ./local/lib/
echo "export LD_LIBRARY_PATH=$ROOT_FOLDER/local/lib" >> ./bin/mgranul
echo "export DYLD_LIBRARY_PATH=$ROOT_FOLDER/local/lib" >> ./bin/mgranul
echo "$ROOT_FOLDER/bin/exe_granul \"\$@\"" >> ./bin/mgranul
chmod +x ./bin/mgranul

rm -rf release
