#!/bin/bash
for value in {0..100}
do
    echo "starting run '$value'"
    python3 run.py 03a-input_poolPairs.csv $value
    echo "finished run '$value'"
done

echo All Done
