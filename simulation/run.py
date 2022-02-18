from user_tools import *

# input_csv = 'simulation/input_parameters.csv' # Input file name
# row_number = 3 # Which row of experiment to run
# input_row = pd.read_csv(input_csv).loc[row_number]
# run_simulations(input_row)

if __name__ == "__main__":
    input_csv = str(sys.argv[1]) # Input file name
    row_number = int(sys.argv[2]) # Which row of experiment to run
    input_row = pd.read_csv(input_csv).loc[row_number]
    run_simulations(input_row)

