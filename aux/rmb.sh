#!/bin/sh

RMBDIR=/tmp/`whoami`_rmb

mkdir $RMBDIR || exit 1

for file in $*
do
  if [ -f $file ]
  then
    tmpfile=$RMBDIR/`basename $file`
    sed -e 's/ *$//' -e 's/	/        /g' $file > $tmpfile
    diff $file $tmpfile > /dev/null 2>&1
    if [ $? = 1 ]
    then
      echo $file
      mv $tmpfile $file
    fi
  fi
done

\rm -rf $RMBDIR

