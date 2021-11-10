#!/bin/bash
for value in {24..99}
do
    echo "starting run '$value'"
    python3 run_poolPairs.py experiment_list.csv $value
    echo "finished run '$value'"
done

echo All Done
