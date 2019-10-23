###################################
###### BUILDING OPENCV 2 ##########
###################################
git clone --depth 1 --branch 2.4.13.7 https://github.com/joeljunior05/opencv.git
CUR_FOLDER=`pwd`
INSTALL=$CUR_FOLDER/local
cd opencv
mkdir ./build
cmake .  -DBUILD_opencv_calib3d=OFF -DBUILD_opencv_ml=OFF -DBUILD_opencv_video=OFF -DBUILD_opencv_legacy=OFF\
         -DBUILD_opencv_objdetect=OFF -DBUILD_opencv_photo=OFF -DBUILD_opencv_gpu=OFF -DBUILD_opencv_ocl=OFF\
         -DBUILD_opencv_nonfree=OFF -DBUILD_opencv_contrib=OFF -DBUILD_opencv_python=OFF -DBUILD_opencv_stitching=OFF\
         -DBUILD_opencv_superres=OFF -DBUILD_opencv_ts=OFF -DBUILD_opencv_videostab=OFF -DBUILD_opencv_core=ON\
         -DBUILD_opencv_features2d=ON -DBUILD_opencv_imgproc=ON -DBUILD_opencv_highgui=ON -DBUILD_opencv_flann=ON\
         -DWITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=$INSTALL . -B./build

cd build
make install
cd $CUR_FOLDER
rm -rf opencv
