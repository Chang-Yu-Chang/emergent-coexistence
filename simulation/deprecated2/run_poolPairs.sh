#!/bin/bash
for value in {0..1000}
do
    echo "starting run '$value'"
    python3 run_poolPairs.py input_poolPairs.csv $value
    echo "finished run '$value'"
done

echo All Done
