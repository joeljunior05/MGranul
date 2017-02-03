#!/bin/bash

FOLDER_NAME=release

########### BUILD SCRIPT FOR MSGRANUL ###############
mkdir $FOLDER_NAME
cd $FOLDER_NAME
cmake ..
make
mv mgranul ../mgranul
