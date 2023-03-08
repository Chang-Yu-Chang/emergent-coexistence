from user_tools import *

if __name__ == "__main__":
    input_csv = str(sys.argv[1]) # Input file name
    input_row = pd.read_csv(input_csv).loc[0]
    
    
    seed = int(input_row['seed'])
    exp_id = int(input_row['exp_id'])

    a = load_assumptions(input_row)
    
    D, c, l = sample_matrices(a)
    S = np.sum(input_row['sa'] * 2)
    # # Load only two species' D
    D[0].to_csv(input_row['output_dir'] + '01-matrices/' + 'D_S' + str(0) + '.csv')
    D[input_row['sa']].to_csv(input_row['output_dir'] + '01-matrices/' + 'D_S' + str(input_row['sa']) + '.csv')
    # Uncomment to load all species
    # for s in range(S):
    #     D[s].to_csv(input_row['output_dir'] + '01-matrices/'+ 'D_S' + str(s) + '.csv')
    
    c.to_csv(input_row['output_dir'] + '01-matrices/' + 'c.csv')
    l.to_csv(input_row['output_dir'] + '01-matrices/' + 'l.csv')


