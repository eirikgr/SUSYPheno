# This script compiles all the needed packages
export SUSYPHENO_PATH=$PWD
echo $SUSYPHENO_PATH
cd $SUSYPHENO_PATH/HDECAY
make
cd $SUSYPHENO_PATH/MICROMEGA/micromegas_3.5.5/
gmake
cd $SUSYPHENO_PATH/MICROMEGA/micromegas_3.5.5/MSSM/
gmake main=micromegas_mini_slha.c
cd $SUSYPHENO_PATH/SUSYHIT/
make
cd $SUSYPHENO_PATH/DARKSUSY/darksusy-5.1.1/
./configure
sed -i -- 's/higgsbounds d/d/g' makefile
make
#make install
cd $SUSYPHENO_PATH


