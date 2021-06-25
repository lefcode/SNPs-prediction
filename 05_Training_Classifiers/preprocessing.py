'''
Script that contains all necessary functions for the preprocessing.
'''

import pandas as pd
import numpy as np
import os
import errno
import random
import warnings
from sklearn.feature_extraction import FeatureHasher
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.impute import KNNImputer
from sklearn.metrics import mean_absolute_error

warnings.filterwarnings('always')  # "error", "ignore", "always", "default", "module" or "once"
pd.options.mode.chained_assignment = None  # default='warn'

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
project_path = os.path.abspath(os.path.join(script_path, os.pardir))
features_path = project_path + "/04_Features_computation/"

DATADIR = script_path + "/data/"
try:
    os.makedirs(DATADIR)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

try:
    os.makedirs(script_path + "/ensembleGASVR/datasets")
except OSError as e:
    if e.errno != errno.EEXIST:
        raise


def encodeCategoricalFeatures(dataset_file):
    """ Creates a list of dictionaries from the labels-values of every column that belongs to the
    categorical features. The dictionary is hashed with the FeatureHasher of sklearn and 2 new columns
    with hashed/encoded values is produced. These columns replace the single column of the feature. """

    print("Encoding categorical features")

    dataset = pd.read_csv(dataset_file, sep="\t", low_memory=False)

    features_to_be_hashed = ["PANTHER-PSEP", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
                             "S", "T",
                             "V", "W", "Y", "U", "O"]
    # list=	["nA","nC","nD","nE","nF","nG","nH","nI","nK","nL","nM","nN","nP",
    #						 "nQ","nR","nS","nT","nV","nW","nY","nU","nO","nX","struct_interactions"]

    for feat in features_to_be_hashed:
        # find the correct n_features value. the value must be multiplied with 2
        # and give result >= #distinct feature values
        # e.x. 3 distinct values 2*n_features>=3 -> n_features=2
        # e.x. 10 distinct values 2*n_features>=10 -> n_features=4
        n_feat = 2
        vals = dataset[feat].value_counts().drop("None", errors='ignore')

        while True:
            if 2 ** n_feat >= len(vals):  # len(dataset[feat].value_counts()):
                break
            else:
                n_feat += 1

        n_feat = 2
        feats_list = [feat + "_" + str(i + 1) for i in range(n_feat)]

        hasher = FeatureHasher(n_features=n_feat, input_type="dict")
        column = dataset[["labels", feat]].set_index('labels')  # .T.to_dict("records")

        # insert each index-column value to a dictionary and put them all in a list
        column_dict_list = list()
        for ind, c in zip(column.index, column.values):
            try:
                column_dict = dict({ind: int(c[0])})
            except:
                column_dict = dict({ind: "None"})
            column_dict_list.append(column_dict)

        column_encoded = hasher.fit_transform(column_dict_list).toarray()
        dataset = dataset.drop(feat, axis=1)
        dataset[feats_list] = column_encoded

    print("Encoding done!")
    return dataset


def mergeSets(pathogenic_dataset, non_pathogenic_dataset_train_set, output_file_name):
    """ Merges the pathogenic and the diminished non-pathogenic datasets into one dataset (merged_dataset)
    which will be used 2/3 for training and 1/3 for testing, along with the rest non-pathogenic records"""

    print("Merging two datasets into {}".format(output_file_name))

    pathogenic_dataset = pd.read_csv(pathogenic_dataset, sep="\t", low_memory=False, index_col=0)
    # pathogenic_dataset["output"] = 1

    non_pathogenic_dataset_train_set = pd.read_csv(non_pathogenic_dataset_train_set, sep="\t", low_memory=False,
                                                   index_col=0)
    # non_pathogenic_dataset_train_set["output"] = 0

    df_complete = pd.concat([pathogenic_dataset, non_pathogenic_dataset_train_set])
    df_complete.to_csv(output_file_name, sep="\t", index=True)

    print("Merging done!")


def diminishSet(full_data):
    """ Randomly selects the 1% of rows of dataset and deletes the value from a randomly selected column """

    full_data = pd.DataFrame(full_data)
    rowIndices = np.random.choice(full_data.shape[0], int(full_data.shape[0] / 100), replace=False)

    nan_replacement_dict = {}
    for i in rowIndices:
        col = random.randint(0, full_data.shape[1] - 1)
        if not pd.isnull(full_data.loc[i][col]):  # only if the value is known
            # keep old values to a dictionary with row,col as key and old values as values.
            nan_replacement_dict[str(i) + "," + str(col)] = full_data.loc[i][col]
            full_data.loc[i][col] = np.nan

    return full_data.to_numpy(), nan_replacement_dict


def applyScaleKnn(dataset, neighbors, indices, scale="standard"):
    """Applies Knn missing values imputation with the best k found in findBestK() function"""

    print("performing knn imputation on dataset with k = " + str(neighbors))

    new_dataset = dataset.set_index("labels").astype(object).replace('None', np.nan)
    cols = [x for x in dataset.columns if x != 'labels']

    output_values = new_dataset["output"].to_numpy()

    if scale == "standard":
        print("Using Standard scaler")
        scaled_dataset = StandardScaler().fit_transform(new_dataset.drop(["output"], axis=1).to_numpy())
    elif scale == "min_max":
        print("Using MinMax scaler")
        scaled_dataset = MinMaxScaler().fit_transform(new_dataset.drop(["output"], axis=1).to_numpy())

    output_values = output_values.reshape(1, output_values.shape[0])
    # reshape output_values np array so it can be concatenated to train_x
    # scaled_dataset = np.concatenate((scaled_dataset, output_values.T), axis=1)

    imputed = KNNImputer(n_neighbors=neighbors, metric="nan_euclidean").fit_transform(scaled_dataset)
    imputed = np.concatenate((imputed, output_values.T), axis=1)

    imputed = pd.DataFrame(imputed, columns=cols, index=indices)
    print(imputed)
    print("knn imputation completed!")
    print("New dataset created!")

    return imputed


def findBestK(data, default_k, scale):
    """ Sets labels ase indices, replaces "None" with np.nan, scales with StandardScaler and
    performs Knn missing values Imputation.
    After using this function the outcome is that k=7 results on the minimum RMSE """

    # assign labels as indexes and replace "None" with np.nan
    print("Applying Standard Scaling and testing possible K values for optimal KNN imputation")

    data = pd.read_csv(data, sep="\t", low_memory=False)  # .drop(["output"], axis=1)
    '''
    new_dataset = data.set_index("labels").astype(object).replace('None', np.nan)
    scaled_dataset = StandardScaler().fit_transform(new_dataset.drop(["output"], axis=1).to_numpy())
    scaled_dataset, nan_replacement_dict = diminishSet(scaled_dataset)

    # finding the K which results on the lowest error
    errors_dict = dict()
    for k in range(1, 21, 2): 
        print("Imputing with k = {}".format(str(k)))

        imputed = KNNImputer(n_neighbors=k, metric="nan_euclidean").fit_transform(scaled_dataset)

        cols = [x for x in data.columns if x != 'labels' and x != 'output']
        df_imputed = pd.DataFrame(imputed, columns=cols)  #, index=new_dataset.index

        predicted = list()
        for key, val in nan_replacement_dict.items():
            row, col = key.split(",")[0], key.split(",")[1]
            predicted.append(df_imputed.loc[int(row)][int(col)])

        mae = mean_absolute_error(predicted, list(nan_replacement_dict.values()))
        print({"K" : k, "Mean absolute error": mae})
        errors_dict[k] = mae

    minimum_error_k = min(errors_dict.keys(), key=(lambda k: errors_dict[k]))
    print('Minimum MAE with K:', minimum_error_k, "MAE:", errors_dict[minimum_error_k])
    '''
    df_imputed = applyScaleKnn(data, neighbors=default_k, indices=data["labels"], scale=scale)  # impute with the best K
    return df_imputed


def splitTrainTestSets(merged_dataset, training_set, testing_set):
    """Maintains 1/3 of pathogenic data and equal amount of non pathogenic data to the train set
     and keeps the rest data as the test set"""

    merged_dataset = pd.read_csv(merged_dataset, sep="\t", low_memory=False)

    path_df = merged_dataset.loc[merged_dataset['output'] == 1]
    non_path_df = merged_dataset.loc[merged_dataset["output"] == 0]

    path_length = len(path_df)
    non_path_length = len(non_path_df) + path_length

    train_count = 0
    # removing 2/3 * #path_length records from the non pathogenic set and insert them to the pathogenic so the training
    # set is created. The rest of the records are considered as test set

    for non_path in non_path_df.values:
        if train_count < 2 * int(path_length / 3):  # train set gets the 1/3 of the pathogenic

            non_path_df.loc[non_path_length + train_count] = path_df.loc[train_count] # insert in test set
            path_df = path_df.drop(labels=train_count, axis=0, inplace=False) # remove from train set

            if train_count % 2 == 1:  # in order to keep only half of the records for the train set
                # every 2 insert in train set
                path_df.loc[path_length + train_count] = non_path
                # remove from test set
                non_path_df = non_path_df.drop(labels=path_length + train_count, axis=0, inplace=False)

            train_count += 1
        else:
            break
    print("Train set consists of {} records".format(str(train_count)))

    path_df.to_csv(training_set, sep='\t', index=False)  # train set
    non_path_df.to_csv(testing_set, sep="\t", index=False)  # test set


def countMinMax(train_set):
    """Finds the min and max values of whole dataset"""

    train_set = pd.read_csv(train_set, sep="\t", low_memory=False).drop(['labels','output'], axis=1)

    max = -1000
    for m in train_set.max():

        if float(m) > max:
            max = float(m)

    print("max value: {}".format(max))
    min = 1000
    for m in train_set.min():
        if float(m) <min:
            min = float(m)
    print("min value: {}".format(min))


def ensembleSplit(train_set, training_set_name, training_labels_name):
    """Format files as required in ensembleGASVR algorithm. Transpose data thus have features in rows
    and proteins in columns and keep labels in separate file"""

    print("Splitting datasets to train and test sets")

    train_set = pd.read_csv(train_set, sep="\t", low_memory=False)#.drop(["labels"], axis=1)

    labels = train_set["output"].values
    features = train_set[train_set.columns.drop(["output"])].values

    data_columns = train_set[train_set.columns.drop(["output"])].columns
    train_y = labels.astype(int).reshape(1, len(labels))

    pd.DataFrame(train_y).to_csv(training_labels_name, sep="\t", index=False, header=False)

    train_x = features.T
    train_set = pd.DataFrame(train_x)
    train_set.index = data_columns

    train_set.to_csv(training_set_name, sep="\t", header=False)


def formatEnsemble(merge, scale, default_k, train_test_split, format, path_set, non_path_set, merged_set,
                   merged_imputed, dataset, train_set_name, train_labels_name):
    """
    Performs merging, scaling with MinMax, train test splitting and formats file as need in ensembleGASVR

    :param merge: use merge function for pathogenic and non pathogenic datasets
    :param scale: scale with MinMax scaler
    :param default_k: default k to perform MinMax scaling
    :param train_test_split: split train and test sets (only for the collected dataset)
    :param format: perform ensemble format
    :param path_set: the pathogenic set (only if to be merged)
    :param non_path_set: the  non pathogenic set (only if to be merged)
    :param merged_set: the merged dataset
    :param merged_imputed: the name of scaled merged set to be stored
    :param dataset: dataset to be split to train_set and train_labels
    :param train_set: the name of train set
    :param train_labels: the name of train set labels
    :return: None
    """
    if merge:
        mergeSets(path_set, non_path_set, merged_set)

    if scale:
        merged_dataset_imputed = findBestK(merged_set, default_k=default_k, scale="min_max")

        merged_dataset_imputed.to_csv(merged_imputed, sep="\t", index=True)

    if train_test_split: # only for the collected dataset
        splitTrainTestSets(merged_imputed,
                                          script_path +"/ensembleGASVR/ensemble_train_set.tsv",
                                          script_path +"/ensembleGASVR/ensemble_test_set.tsv")

    if format:
        """Format files as required in ensembleGASVR algorithm. Transpose data thus have features in rows
        and proteins in columns and keep labels in separate file"""
        print("Splitting dataset to data and labels")

        dataset = pd.read_csv(dataset, sep="\t", low_memory=False)  # .drop(["labels"], axis=1)

        labels = dataset["output"].values
        features = dataset[dataset.columns.drop(["output"])].values

        data_columns = dataset[dataset.columns.drop(["output"])].columns
        train_y = labels.astype(int).reshape(1, len(labels))

        pd.DataFrame(train_y).to_csv(train_labels_name, sep="\t", index=False, header=False)

        train_x = features.T
        dataset = pd.DataFrame(train_x)
        dataset.index = data_columns
        print(dataset.shape)
        dataset.to_csv(train_set_name, sep="\t", header=False)


def createEnsemble():
    """Calls format ensemble function for all datasets"""

    # train set
    formatEnsemble(merge=True, scale=True, default_k=1, train_test_split=True, format=True,
                   path_set=DATADIR + "path_dataset.tsv", non_path_set=DATADIR + "non_path_dataset.tsv",
                   merged_set=DATADIR + "merged_dataset.tsv",
                   merged_imputed=script_path +"/ensembleGASVR/ensemble_merged_dataset_imputed.tsv",
                   dataset=script_path + "/ensembleGASVR/ensemble_train_set.tsv",
                   train_set_name=script_path + "/ensembleGASVR/datasets/ensemble_training_set.tsv",
                   train_labels_name=script_path + "/ensembleGASVR/datasets/ensemble_training_set_labels.tsv")

    # test set
    formatEnsemble(merge=False, scale=False, default_k=1, train_test_split=False, format=True,
                   path_set="", non_path_set="", merged_set="", merged_imputed="",
                   dataset=script_path + "/ensembleGASVR/ensemble_test_set.tsv",
                   train_set_name=script_path + "/ensembleGASVR/datasets/ensemble_testing_set.tsv",
                   train_labels_name=script_path + "/ensembleGASVR/datasets/ensemble_testing_set_labels.tsv")

    # humvar
    formatEnsemble(merge=False, scale=True, default_k=3, train_test_split=False, format=True,
                   path_set="", non_path_set="",
                   merged_set=DATADIR + "humvar_dataset.tsv",
                   merged_imputed=script_path + "/ensembleGASVR/ensemble_humvar_dataset_imputed.tsv",
                   dataset=script_path + "/ensembleGASVR/ensemble_humvar_dataset_imputed.tsv",
                   train_set_name=script_path + "/ensembleGASVR/datasets/ensemble_humvar.tsv",
                   train_labels_name=script_path + "/ensembleGASVR/datasets/ensemble_humvar_labels.tsv")

    # humdiv
    formatEnsemble(merge=False, scale=True, default_k=7, train_test_split=False, format=True,
                   path_set="", non_path_set="",
                   merged_set=DATADIR + "humdiv_dataset.tsv",
                   merged_imputed=script_path + "/ensembleGASVR/ensemble_humdiv_dataset_imputed.tsv",
                   dataset=script_path + "/ensembleGASVR/ensemble_humdiv_dataset_imputed.tsv",
                   train_set_name=script_path + "/ensembleGASVR/datasets/ensemble_humdiv.tsv",
                   train_labels_name=script_path + "/ensembleGASVR/datasets/ensemble_humdiv_labels.tsv")

if __name__ == "__main__":
    createEnsemble()




















