#!/bin/bash

#use this file to batch download bootstrap data from the nsidc website
#in order to use a cookie must be set for wget, do the following:
# wget --http-user=strausz --ask-password --save-cookies ~/.urs_cookies https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0079.003/1978.11.01/bt_19781101_n07_v3.1_n.bin
#this will download one file and your cookie should now be set
#so you can proceed to use this sript.  This makes it so you don't have to save passwords in a file


#wget options
options="--load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -nd -e robots=off"
url="https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0079.003/"

#specify start month and desired year
#month=11  doesn't work yet specify month in for loop below
year=2018

for i in {11..12}
do
  month=$(printf "%02d" $i)
  for x in {1..31}
  do
    day=$(printf "%02d" $x)
    filename="bt_${year}${month}${day}_f17_v3.1_n.bin"
    directory="${year}.${month}.${day}"
    echo "executing wget $options ${url}${directory}/${filename}"
    wget $options ${url}${directory}/${filename}
    
  done
  
done
