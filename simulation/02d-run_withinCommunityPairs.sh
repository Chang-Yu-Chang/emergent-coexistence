#!/bin/bash
for value in {0..1000}
do
    echo "starting run '$value'"
    python3 run.py 01d-input_withinCommunityPairs.csv $value
    echo "finished run '$value'"
done

echo All Done
