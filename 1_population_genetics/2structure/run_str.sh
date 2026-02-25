#!/usr/bin/bash
num=$1;  #num of K

data="Gpen66_mis20_13129snp_str.txt"  #input
rdir="result"
logdir="resultlog"

mkdir -p $rdir   #result dir
mkdir -p $logdir   #log file dir
for i in {1..20}
{
     random=$(echo $RANDOM)     
     fileout="result_"$num"-"$i
     filelog="log"$num"-"$i
     nohup structure -K $1 -D $random -i $data -o ./$rdir/$fileout > ./$logdir/$filelog 2>&1 &
}
