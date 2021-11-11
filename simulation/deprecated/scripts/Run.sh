for value in {0..4}
do
    echo "starting run '$value'"
    python Run_GlobalPool.py Experiment_list.csv $value 
    echo "finished run '$value'"
done

echo All Done

