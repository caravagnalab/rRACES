#!/usr/bin/env bash

PKGNAME=$(sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS=$(sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC=$(basename `pwd`)
PKG_PACKAGE_NAME="${PKGNAME}_${PKGVERS}.tar.gz"

CHECK_LOG="../${PKGNAME}.Rcheck/00check.log"

rm -rf "../${PKG_PACKAGE_NAME}"

(cd .. && R CMD build ${PKGSRC})

if [ ! -f ../${PKG_PACKAGE_NAME} ]
then
 echo "\"R CMD build\" has failed: commit not allowed!"
 exit 1
fi

(cd .. && R CMD check ${PKG_PACKAGE_NAME} --as-cran)

if [ ! -f ${CHECK_LOG} ]
then
 echo "\"R CMD check\" has failed: commit not allowed!"
 exit 1
fi

if [ $(grep ERROR ${CHECK_LOG} | wc -l) -ne 0 ]
then
 echo "\"R CMD check\" produces errors: commit not allowed!"
 echo "Please, check the file ${CHECK_LOG}"
 exit 1
fi

if [ $(grep WARNING ${CHECK_LOG} | wc -l) -ne 0 ]
then
 echo "\"R CMD check\" produces warnings: commit not allowed!"
 echo "Please, check the file ${CHECK_LOG}"
 exit 1
fi
