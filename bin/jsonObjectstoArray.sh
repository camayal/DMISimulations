#!/bin/bash

#Convert list of objects into an unique array in the JSON

for file in $@
do
echo "Organizing JSONs objects in one array in ${file}..."
sed -i '/\[{/d' $file
sed -i '/}]/d' $file
sed -i '/^$/d' $file
sed -i '1s/^{/[{/' $file
sed -i '$ s/},$/}]/' $file
done


