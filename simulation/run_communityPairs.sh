#!/bin/bash
for value in {0..1000}
do
    echo "starting run '$value'"
    python3 run_communityPairs.py input_communityPairs.csv $value
    echo "finished run '$value'"
done

echo All Done
