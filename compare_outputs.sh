#!/bin/sh

for FILE in test_dis.csv test_dis_interv.csv; do
    (grep -E '.*,(fe)?male,(2|42),2010' $FILE \
        | awk 'BEGIN { FS="," } { printf "%-6s %-3d %0.5f\n", $2, $3, $7 }' \
        | tac)
     echo
done
