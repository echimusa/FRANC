#! /bin/bash

prefix="/usr/local/lamp"
export PLPLOT_LIB="$prefix/bin/fnt"
echo $PLPLOT_LIB
ancfile=$1

$prefix/bin/generategraph $ancfile
cp $prefix/bin/index.html .
