# This script compiles all the needed packages
export SUSYPHENO_PATH=$PWD
echo $SUSYPHENO_PATH
cd $SUSYPHENO_PATH/HDECAY
make
cd $SUSYPHENO_PATH/MICROMEGA/micromegas_5.0.8/
make
cd $SUSYPHENO_PATH/SUSYHIT/
make
cd $SUSYPHENO_PATH/DARKSUSY/darksusy-5.1.1/
./configure
make
cd $SUSYPHENO_PATH


