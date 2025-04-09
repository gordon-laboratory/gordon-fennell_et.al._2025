# libraries
## standard libraries
import sys
import numpy as np
import pandas as pd
import math
from matplotlib import pyplot as plt
import os
#%%

## custom libraries
sys.path.append(r'.\functions')
from glm import get_ynull_wrapping
from glm import compute_pvalues_null_data


# GLM for multi-spout data
#%%
# notes
# current organization
#  indiviudal sessions located in ./data/sessions/
#  combined predictors  located in ./data/sessions_combined_subject/##/glm/
#  combined signals  located in ./data/sessions_combined_subject/##/glm/signals/


#%%
def ensure_directory_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

def get_folders(directory):
    return [name for name in os.listdir(directory) if os.path.isdir(os.path.join(directory, name))]

def get_file_names(directory):
    return [name for name in os.listdir(directory) if os.path.isfile(os.path.join(directory, name))]


def read_data(base_dir, file_names):
    data = {}
    for file_name in file_names:
        file_path = f"{base_dir}{file_name}"
        data[file_name[:-4]] = pd.read_csv(file_path, header=None).to_numpy()

    return data

def repeat_array_horizontally(array, n):
    # Replicate the array horizontally n times
    return np.tile(array, (1, n))



def padd_and_trim_null_vertically(Y_null, Y_all):
    Y_null_padded =  Y_null

    # stack Y_null to match or exceed rows in Y_all
    n = len(Y_all) / len(Y_null)  # Calculate how many times to stack Y_null

    for n_stack in np.arange(1, math.ceil(n)):

        # shuffle columns of Y_null and bind vertically
        shuffled_columns = np.random.permutation(Y_null.shape[1])

        Y_null_padded = np.vstack([Y_null_padded, Y_null[:, shuffled_columns]])

    # trim to length of Y_null to match Y_all
    Y_null_padded = Y_null_padded[np.arange(0, len(Y_all)),:]

    return(Y_null_padded)


#%%
dir_sets = './data/sessions_combined_subject/'
dir_nulls = './data/glm_null_distributions/'
subject_sets = get_folders(dir_sets)
toggle_filter_lick = 1

blockname_region_filtered = pd.read_csv('key_blockname_region_filtered.csv')



for subject_set in subject_sets:
    print(subject_set)

    # create output folder to store glm results
    dir_output = dir_sets + subject_set + '/analysis_output_glm/'
    ensure_directory_exists(dir_output)

    # read in predictors
    dir_predictors = dir_sets + subject_set + '/predictors/'
    fn_predictors = get_file_names(dir_predictors)

    data_x = read_data(dir_predictors, fn_predictors)


    predictor_ids = [
        'lick_kernal',
        'diagonal_true_trial',
        'diagonal_true_solution_conc_scaled_0_p1',
        'diagonal_true_solution_conc_scaled_0_p1_history03'
    ]

    # combine predictors defined in predictor_ids

    X_full = []
    for predictor_id in predictor_ids:
        X_full.append(data_x[predictor_id])

    # read in categorical info
    df_cat = pd.read_csv(dir_sets + subject_set + '/sample_info.csv')

    if(toggle_filter_lick):
        df_cat = df_cat[df_cat['trial_lick'] == 1]


    # get signal file names
    dir_signals = dir_sets + subject_set + '/signals/'
    fn_signals = get_file_names(dir_signals)

    # filter to zscoreblsub
    fn_signals = [file for file in fn_signals if 'zscoreblsub' in file]

    for signal in fn_signals:
        print(' -' + signal)

        Y_all = pd.read_csv(dir_signals + signal, header=None).to_numpy()

        # for signals with individual sessions filtered
        blocks_in_df_cat = df_cat['blockname'].unique()
        signal_id = signal[:signal.rfind("_")]

        # for str fibers, filter blockname_region_filtered based on signal_id
        if signal_id.find("lha") == -1:
            blockname_region_filtered_signal_id = blockname_region_filtered[blockname_region_filtered['region_original'] == signal_id]
        else:
            blockname_region_filtered_signal_id = blockname_region_filtered

        blocks_filtered = blockname_region_filtered_signal_id[blockname_region_filtered_signal_id['blockname'].isin(blocks_in_df_cat)]
        filt_index = df_cat['blockname'].isin(blocks_filtered['blockname'])

        # filter Y and X based on data included in signal
        X_filt = [matrix[filt_index, :] for matrix in X_full]

        # read in corresponding null matrix (produced in r)
        null_prefix = signal.split('.', 1)[0]
        Y_null = pd.read_csv(dir_nulls + null_prefix + '_null_data.csv', header=None).to_numpy()

        Y_null = padd_and_trim_null_vertically(Y_null, Y_all)

        Y_all_padded = repeat_array_horizontally(Y_all, Y_null.shape[1])

        Y_null_wrap = get_ynull_wrapping(Y_all_padded)

        if(len(X_filt[1]) == len(Y_all_padded)):

            pvals_wrap, delta_r2_wrap = compute_pvalues_null_data(
                X_list=X_filt,
                Y_full=Y_all_padded,
                Y_null=Y_null_wrap
            )

            pvals_wrap_df = pd.DataFrame(pvals_wrap, columns=np.array(predictor_ids))
            filename = os.path.join(dir_output, f'pvals_{signal}')
            pvals_wrap_df.head(1).to_csv(filename, index=False)

            delta_r2_wrap_df = pd.DataFrame(delta_r2_wrap, columns=np.array(predictor_ids))
            filename = os.path.join(dir_output, f'deltar2_{signal}')
            delta_r2_wrap_df.head(1).to_csv(filename, index=False)

        else:
            print('   * error, pred/sig length missmatch')



































