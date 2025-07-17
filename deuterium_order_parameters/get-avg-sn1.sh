#!/bin/bash

# Script is generated to combine and calculate the average SCD by Xiaohong Zhuang on 7/25/2014
# The transpose code is from online source: http://stackoverflow.com/questions/1729824/transpose-a-file-in-bash

unset noclobber

#transpose the files of column radius to row
# ADJUST BOUNDARIES FOR LOOP
for (( i=1; i<=10; i++ ))

do

# select the pore radius column
 cat c3-$i.dat | awk {'print $2'} > c3-$i-scd.dat

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}'  c3-$i-scd.dat > 1-tr-c3-$i-scd.dat

rm -f c3-$i-scd.dat 
done

# combine the the dyn files
cat 1-tr-c3-*-scd.dat > 1-tr-c3-scd-alldyn.dat

rm -f 1-tr*scd.dat 

# Calculate average and stderr
awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} 
          END {for (i=1;i<=NF;i++) {
          print sum[i]/NR," ",sqrt((sumsq[i]-sum[i]^2/NR)/NR^2)}
         }' 1-tr-c3-scd-alldyn.dat >> 1-tr-c3-scd-alldyn_avg.dat

#transpose the files 
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}'  1-tr-c3-scd-alldyn_avg.dat > 1-tr2-c3-scd-alldyn_avg.dat

# add title line

cat header-sn1.dat 1-tr2-c3-scd-alldyn_avg.dat > 1-tr2-c3-scd-alldyn_avg_title.dat

#transpose the files 
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}'  1-tr2-c3-scd-alldyn_avg_title.dat > sn1-avg.dat

rm -f 1-tr2*.dat 
rm -f 1-tr*alldyn.dat 
rm -f 1-tr-c3-scd-alldyn_avg.dat 
