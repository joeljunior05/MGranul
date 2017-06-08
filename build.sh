#!/bin/bash

FOLDER_NAME=release

########### BUILD SCRIPT FOR MSGRANUL ###############

if [ ! -d "$FOLDER_NAME" ]; then
  mkdir $FOLDER_NAME

  cd $FOLDER_NAME
  cmake ..
  make
  mv mgranul ../mgranul

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
    mv mgranul ../mgranul
  else
    echo "$FOLDER_NAME not erased."
  fi
fi
