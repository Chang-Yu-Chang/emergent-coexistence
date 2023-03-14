#!/bin/bash
for value in {0..100}
do
    echo "starting run '$value'"
    python3 run.py 04a-input_withinCommunityTrios.csv $value
    echo "finished run '$value'"
done

echo All Done
