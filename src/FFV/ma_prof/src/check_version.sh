#!/bin/sh

# Copyright (C) 2014 RIKEN AICS
# This library is released under the terms of the MIT license.
# http://fiber-miniapp.mit-license.org/

if [ $# -ne 1 ]
then
        echo "usage: $0 command" 1>&2
        exit 1
fi

if type $1 > /dev/null 2>&1
then
        :
else
        echo "error $0: command $1 not exist." 1>&2
        exit 1
fi

case $1 in
        *frt|*fcc|*FCC|*frtpx|*fccpx|*FCCpx)  # Fujitsu
                $1 -V 2>&1 | head -1
                ;;
        gfortran|g95|gcc|g++)  # GCC
                $1 --version | head -1
                ;;
        *ifort|*icc|*icpc)  # Intel
                $1 --version | head -1
                ;;
        pg*)  # PGI
                $1 --version | head -2 | tail -1
                ;;
        nagfor)  # NAG
                $1 -V 2>&1 | head -1
                ;;
        f90|f77|cc|c++|cxx|mpif77|mpif90|mpicc|mpiCC|mpic++|mpicxx)
                if tmp=`$1 --version 2> /dev/null`
                then
                        # echo first nonblank line
                        echo "$tmp" | while read line;
                        do
                                if [ -n "$line" ]
                                then
                                        echo $line
                                        break
                                fi
                        done
                elif tmp=`$1 -V 2> /dev/null`   # Sun
                then
                        echo "$tmp" | head -1
               #elif ...
               #then
               #        :
                else
                        echo "Unknown"
                fi 
                ;;
        *)
                echo "Unknown"
                ;;
esac
