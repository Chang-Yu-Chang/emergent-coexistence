#!/bin/bash
for value in {0..1000}
do
    echo "starting run '$value'"
    python3 run.py input_independent.csv $value
    echo "finished run '$value'"
done

echo All Done
