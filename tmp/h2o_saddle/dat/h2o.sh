#!/bin/csh

set DATA=/home/rtake/h2o_saddle
set WORK=/scr/rtake.h2o.$$
set GRRM=/home/GRRMdv0123b

setenv subgrr ${GRRM}/GRRM.out
setenv submpi mpirun
setenv g09root /usr/local/gaussian09d01
setenv GAUSS_SCRDIR /scr
setenv OMP_NUM_THREADS 1
source $g09root/g09/bsd/g09.login
setenv subgau g09
setenv subchk formchk
setenv submol "molpro -W${WORK} -d${WORK}"
setenv subsiesta siesta
setenv siestascr ${WORK}
setenv TURBODIR /usr/local/TURBOMOLE7.0/TURBOMOLE
set path = ( $TURBODIR/scripts $path)
setenv PARA_ARCH SMP
set path = ($TURBODIR/bin/`sysname` $path)
setenv subgms /usr/local/gamess/rungms
setenv gmstmp ~/scr
setenv gmsscr1 ${WORK}
setenv gmsscr2 /home/rtake/scr
if (! -d /scr/rtake ) mkdir -m 700 /scr/rtake
setenv subqchem /usr/local/qchem-4.1.2/exe/qcprog.exe
setenv qchemtmp ${WORK}
setenv QC /usr/local/qchem-4.1.2
setenv subdalton /home/grrm_admin/buchi/dalton/buchi_dalton/dalton
setenv daltontmp ${WORK}
setenv submndo /usr/local/mndo99/mndo99_20161021_intel64_ifort-13.1.3.192_mkl-11.1.4.211
setenv subvasp /usr/local/VASP-5.3.3-22May2013/bin/vasp
setenv subdftb_plus /usr/local/dftb_plus/bin/dftb+ 
alias cp 'cp -pf'

if ( -d ${WORK} ) rm -rf ${WORK}
mkdir -m 700 ${WORK}

mv ${DATA}/h2o.log ${WORK}
mv ${DATA}/h2o_* ${WORK}
cp ${DATA}/h2o.* ${WORK}
cp ${DATA}/*.inp ${WORK}
cp ${DATA}/*.psf ${WORK}
cp ${DATA}/*.hsd ${WORK}
cp ${DATA}/h2o/* ${WORK}

cd ${WORK}
${GRRM}/GRRMp h2o -p1 -s172800

cp ${WORK}/h2o.* ${DATA}
cp ${WORK}/h2o_* ${DATA}

cd ${DATA}
rm -rf ${WORK}



