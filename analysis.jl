base_path = "C:\\Users\\User\\Dropbox\\CLOexternalBrain\\projects\\MA Thesis\\code\\"
include(joinpath(base_path, "causal_stream_functions.jl"))

raw_data_path = joinpath(base_path, "data\\raw\\")
results_data_path = joinpath(base_path, "data\\results\\")

TPM_hard_path = joinpath(raw_data_path, "TPM_hard.csv")
TPM_easy_path = joinpath(raw_data_path, "TPM_easy.csv")
state_data_hard_path = joinpath(raw_data_path, "state_data_hard.csv")
state_data_easy_path = joinpath(raw_data_path, "state_data_easy.csv")

TPM_hard = DataFrame(CSV.File(TPM_hard_path))
TPM_easy = DataFrame(CSV.File(TPM_easy_path))
state_data_hard = DataFrame(CSV.File(state_data_hard_path))
state_data_easy = DataFrame(CSV.File(state_data_easy_path))

n_nodes = 8
trials = [6,24,38,57,70,91,102,122]

trial_streams_hard, trial_streams_strength_hard = get_trial_streams(n_nodes, trials, state_data_hard, TPM_hard)
trial_streams_easy, trial_streams_strength_easy = get_trial_streams(n_nodes, trials, state_data_easy, TPM_easy)

stream_perspectives_hard = get_stream_perspectives_IO(trial_streams_hard, trial_streams_strength_hard)
stream_perspectives_easy = get_stream_perspectives_IO(trial_streams_easy, trial_streams_strength_easy)

CSV.write(joinpath(results_data_path, "causal_streams_hard.csv"), stream_perspectives_hard)
CSV.write(joinpath(results_data_path, "causal_streams_easy.csv"), stream_perspectives_easy)

