#!/bin/bash
##############################################
#  this script tests single configurations
#  to be tested with all eight configurations
#
#  make sure supply bagel input job name and that
#  your distribution directories (smithdir and bageldir)
#  are correctly set.
#
#  to run, for example
#
#  ./test-singles.sh y.in > singles.log
#
##############################################

smithdir=~/smith
bageldir=~/bagel
input=$1

# must give input file name
if [[ -z "$@" ]]; then
    echo >&2 "You must supply a BAGEL input file name!"
    exit 1
fi

# check to see if input exists
if [[ ! -f "$bageldir/obj/src/$input" ]]; then
    echo >&2 "Your input file doesn't exist!"
    echo >&2 "File needed:"
    echo >&2 "$bageldir/obj/src/$input"
    exit 1
fi

# Single configuration test cases. Order is important, corresponds to original Roos ordering.
cases="ccaa xcaa xxaa ccxa cxxa ccxx xxxa xcxx"

## run script ###

echo "Testing input: "
echo $bageldir/obj/src/$input
cat $bageldir/obj/src/$input

for i in `echo "$cases"`
 do
   echo "Starting test case: $i"
   sed -i "/[/][/] *[*]test single configuration cases[*]/,/[/][/] *[*]end test single configuration cases[*]/s/\(^  \)\( *if.*\)/\/\/\2/;/[/][/] *${i} */s/^ *[/][/]/  /" $smithdir/prep/generate_main.cc
   cd $smithdir/obj
   make -j
   ./prep/Prep > ../main.cc
   make -j
   ./MulWick

   # copy output to bagel, compile and run
   cp CAS_test* $bageldir/src/smith
   cd $bageldir/obj/src
   make
   ./BAGEL "$input"
   echo "Energy was for case: $i"
   cd $smithdir/obj
 done

