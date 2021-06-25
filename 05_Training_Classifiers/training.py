'''
Necessary files: the created pathogenic and non pathogenic datasets from the previous step 04.Features_computation

Script that calls the necessary function that execute the preprocessing step and the training of Random Forest
classifier and the testing of the model
1. Encode categorical features for both pathogenic and non pathogenic sets, create output column and save file
2. Merge pathogenic and non pathogenic sets
3. Keep 1/3 of pathogenic and equal number of non pathogenic records to train set and store all the rest to test set
4. Find best k for KNN imputation of missing values and perform KNN with best found k and save file
5. Perform Hyper-parameter tuning for the Random Forest Classifier
6. Train Random Forest Classifier with best resulting values
7. Evaluate model, plot the 30 high-scoring features and create equivalent datasets with these features,
   plot the AUROC and visualize the forest
8. Predict test set labels and evaluate the model
'''

import pandas as pd
import numpy as np
import os
import errno
import pickle
from matplotlib import pyplot

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, KFold, learning_curve, cross_val_score, GridSearchCV
from sklearn.metrics import f1_score, accuracy_score, mean_absolute_error, matthews_corrcoef, balanced_accuracy_score,\
    precision_score, recall_score, roc_auc_score, fbeta_score, roc_curve, confusion_matrix, classification_report,\
    plot_roc_curve, auc
from imblearn.metrics import geometric_mean_score
from sklearn import tree
import preprocessing
import tuning

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
project_path = os.path.abspath(os.path.join(script_path, os.pardir))
features_path = project_path + "/04_Features_computation/"

DATADIR = script_path + "/data/"
RESULTSDIR = script_path + "/results/"
try:
    os.makedirs(RESULTSDIR)
except OSError as e:
    if e.errno != errno.EEXIST: raise


def trainRandomForest(train_set, model_name):
    """Trains Random Forest Classifier with the train set"""

    # tuning.tuneRandomForest(train_set)

    train_set = pd.read_csv(train_set, sep="\t", low_memory=False)

    train_output = train_set["output"].values
    train_features = train_set[train_set.columns.drop(["labels", "output"])].values

    kf = KFold(n_splits=5, shuffle=True, random_state=42)  # Use for 5-Fold classification

    n_estimators = 2000
    max_features = "sqrt"
    criterion = "gini"
    bootstrap = True # True only if oob_score = True else False
    oob_score = True # True only if bootstrap = True else False
    # not used
    min_samples_split = 2 # None or 2 (default value)
    min_samples_leaf = 1 # None or 1 (default value)
    warm_start = False  # or None
    ccp_alpha = 0 # or None

    rf = RandomForestClassifier(n_jobs=-1, n_estimators=n_estimators, max_features=max_features,
                                criterion=criterion, bootstrap=bootstrap, oob_score=oob_score, 
                                max_samples=None)

    print("""Training Random Forest with the following parameters: 
    n_estimators = {}
    max_features = {}
    criterion = {}
    bootstrap = {}
    oob_score = {}
    """.format(n_estimators, max_features, criterion, bootstrap, oob_score))

    mae_list = list()
    precision_list = list()
    recall_list = list()
    specificity_list=list()
    f1_score_list = list()
    area_under_curve_list = list()
    accuracy_list = list()
    f2_score_list = list()

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    i = 0
    for train_index, test_index in kf.split(train_features):
        X_train, X_test, y_train, y_test = train_features[train_index], train_features[test_index], \
                                           train_output[train_index], train_output[test_index]

        rf.fit(X_train, y_train)
        prediction = rf.predict(X_test)

        fpr, tpr, t = roc_curve(y_test, prediction)
        tprs.append(np.interp(mean_fpr, fpr, tpr))
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)

        mae, precision, recall, specificity, f1, area_under_curve, accuracy, f2 = \
            evaluateModel(prediction, y_test, show_results=False)

        mae_list.append(mae)
        precision_list.append(precision)
        recall_list.append(recall)
        specificity_list.append(specificity)
        f1_score_list.append(f1)
        area_under_curve_list.append(area_under_curve)
        accuracy_list.append(accuracy)
        f2_score_list.append(f2)

    print("RF classifier trained!")
    print("Average Mean absolute error (MAE): %f" % np.mean(mae_list))
    print("Average Precision Score : %f" % np.mean(precision_list))
    print("Average Recall Score : %f" % np.mean(recall_list))
    print("Average Specificity Score : %f" % np.mean(specificity_list))
    print("Average F1 Score : %f " % np.mean(f1_score_list))
    print("Average AUC :  %f " % np.mean(area_under_curve_list))
    print("Average Accuracy :  %f " % np.mean(accuracy_list))
    print("Average F2 score :  %f" % np.mean(f2_score_list))

    saveModel(rf, model_name, train_set.columns)
    plotAUC(tprs, mean_fpr, "Validation set", 'blue', merge=True)


def saveModel(model, model_name, columns):
    """Saves the model and plots features importance graph and a file"""

    pickle.dump(model, open(model_name, 'wb'))  # save the model via pickle

    try:

        if not columns.empty:
            # plot the most important features in the random forest
            pyplot.rcParams.update({'figure.figsize': (30, 30)})
            pyplot.rcParams.update({'font.size': 30})

            sorted_idx = model.feature_importances_.argsort()
            features = list()
            importances = list()
            with open(RESULTSDIR + "/features_importances.txt", "w+") as importances_file, \
                    open(RESULTSDIR + "/important_features.txt", "w+") as important:

                importances_file.write("Feature name\tImportance\n")
                important.write("Feature name\tImportance\n")
                important.write("labels\t0.0\n")

                for c, imp in zip(columns[sorted_idx], model.feature_importances_[sorted_idx]):
                    importances_file.write(c + "\t" + str(imp) + "\n")
                    if imp >= 0.0078:  # we keep the most important features not all!
                        features.append(c)
                        importances.append(imp)
                        important.write(c + "\t" + str(imp) + "\n")

                important.write("output\t0.0\n")

            pyplot.title('Most important features', fontsize=35)
            pyplot.barh(features, importances)
            pyplot.ylabel("Features")
            pyplot.xlabel("Significance")
            pyplot.savefig(RESULTSDIR + "/features_importance.png")
            # pyplot.show()
            pyplot.close()

            print("features_importance.png created!")
            # visualizeForest(model)
            # extractImportantFeaturesDatasets()
    except:
        pass


def visualizeForest(rf):
    """Creates visualization of the first 5 trees of the forest"""

    cn = ["0", "1"]
    fig, axes = pyplot.subplots(nrows=1, ncols=5, figsize=(20, 15), dpi=900)
    for index in range(0, 5):
        tree.plot_tree(rf.estimators_[index], class_names=cn, filled=True, ax=axes[index])
        axes[index].set_title('Estimator: ' + str(index), fontsize=11)

    fig.savefig(RESULTSDIR + '/rf_5trees.png')
    print("Trees visualization complete!")


def evaluateModel(prediction, y_test, show_results=True):
    """ Computes the metrics to evaluate the model's correctness """

    mae = mean_absolute_error(y_test, list(map(int, prediction)))
    precision = precision_score(y_test, list(map(int, prediction)), average='weighted', labels=np.unique(prediction))
    recall = recall_score(y_test, list(map(int, prediction)), average='weighted', labels=np.unique(prediction))

    tn, fp, fn, tp = confusion_matrix(y_test, prediction).ravel()
    specificity = tn / (tn + fp)

    f1 = f1_score(y_test, list(map(int, prediction)), average='weighted', labels=np.unique(prediction))
    try:
        area_under_curve = roc_auc_score(y_test, list(map(int, prediction)))
    except ValueError: # in case the testing set includes only data from one class
        area_under_curve = 0

    accuracy = accuracy_score(y_test, list(map(int, prediction)))
    f2 = fbeta_score(y_test, list(map(int, prediction)), average='weighted', beta=0.5)

    #bacc = balanced_accuracy_score(y_test, list(map(int, prediction)))

    if show_results:
        print("Mean absolute error (MAE): %f" % mae)
        print("Precision Score : ", precision)
        print("Recall Score : ", recall)  # sensitivity
        print("Specificity Score : ", specificity)
        print("F1 Score : ", f1)
        print("AUC : ", area_under_curve)
        print("Accuracy : ", accuracy)
        print("F2 score : ", f2)

    return mae, precision, recall, specificity, f1, area_under_curve, accuracy, f2


def plotAUC(tprs, mean_fpr, dataset_name, color, merge=False):
    """Plots all ROC for validation, test, humvar and humdiv sets"""

    mean_tpr = np.mean(tprs, axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)

    pyplot.plot(mean_fpr, mean_tpr, color=color,
                label=dataset_name + " ROC (AUC = %0.2f )" % (mean_auc), lw=2, alpha=1)
        # r'Mean ROC (AUC = %0.2f )'

    pyplot.rcParams.update({'figure.figsize': (30, 30)})
    pyplot.rcParams.update({'font.size': 13})
    pyplot.xlabel('False Positive Rate')
    pyplot.ylabel('True Positive Rate')
    pyplot.title('ROC')

    if merge:  # in case the user want to merge all ROCs to one plot
        pyplot.legend(['Train', 'Test', "HumVar", "HumDiv"], loc="lower right")
    else:
        pyplot.legend(dataset_name, loc="lower right")

    pyplot.savefig(RESULTSDIR + "AUROC.png")
    if "HumDiv" in dataset_name:
        pyplot.close()
    '''
        pyplot.plot(mean_fpr, mean_tpr, color=color,
                    label=dataset_name + " ROC (AUC = %0.2f )" % (mean_auc), lw=2, alpha=1)
        # r'Mean ROC (AUC = %0.2f )'
        pyplot.rcParams.update({'figure.figsize': (30, 30)})
        pyplot.rcParams.update({'font.size': 13})
        pyplot.xlabel('False Positive Rate')
        pyplot.ylabel('True Positive Rate')
        pyplot.title('ROC')
        pyplot.legend(dataset_name, loc="lower right")
        pyplot.savefig(RESULTSDIR + dataset_name + "_AUROC.png")
        pyplot.close()
    '''

def predictTestSet(test_set, model_name):
    """Does prediction with the test set and evaluates the model"""

    dataset_name = "Test set"
    color = 'red'
    if "humvar" in test_set:
        dataset_name = "HumVar"
        color = "orange"
    elif "humdiv" in test_set:
        dataset_name = "HumDiv"
        color = "black"

    test_set = pd.read_csv(test_set, sep="\t", low_memory=False)
    test_output = test_set["output"].values
    test_features = test_set[test_set.columns.drop(["labels","output"])].values

    loaded_model = pickle.load(open(model_name, 'rb'))

    pred = loaded_model.predict(test_features)

    mean_fpr = np.linspace(0, 1, 100)
    fpr, tpr, t = roc_curve(test_output, pred)
    tprs = [np.interp(mean_fpr, fpr, tpr)]
    plotAUC(tprs, mean_fpr, dataset_name, color, merge=True)

    evaluateModel(pred, test_output, show_results=True)


def extractImportantFeaturesDatasets():
    """Creates datasets based on the computed most important features"""

    exrtactImportantFeatures(dataset=script_path + "/data/train_set2.tsv",
                             important_features=script_path + "/results/important_features.txt",
                             output_file=script_path + "/results/important_features_train_set.tsv")

    exrtactImportantFeatures(dataset=script_path + "/data/test_set2.tsv",
                             important_features=script_path + "/results/important_features.txt",
                             output_file=script_path + "/results/important_features_test_set.tsv")

    exrtactImportantFeatures(dataset=script_path + "/data/humvar_dataset.tsv",
                             important_features=script_path + "/results/important_features.txt",
                             output_file=script_path + "/results/humvar_important_features.tsv")

    exrtactImportantFeatures(dataset=script_path + "/data/humdiv_dataset.tsv",
                             important_features=script_path + "/results/important_features.txt",
                             output_file=script_path + "/results/humdiv_important_features.tsv")


def exrtactImportantFeatures(dataset, important_features, output_file):

    important_features = pd.read_csv(important_features, sep="\t")
    dataset = pd.read_csv(dataset, sep="\t")

    features = important_features["Feature name"]#.append(pd.Series('labels'))#.append("labels")

    dataset[features].to_csv(output_file, sep="\t",
                                            index=False)


def runPreprocessTrain():
    """Gives options for preprocessing, training and evaluation of the classifier """

    '''
    # for ploting the AUROC 
    trainRandomForest(train_set=DATADIR + "train_set2.tsv",
                      model_name=script_path + "/random_forest_model.sav")

    predictTestSet(DATADIR + "test_set2.tsv", script_path + "/random_forest_model.sav")
    predictTestSet(DATADIR + "humvar_dataset_imputed.tsv", script_path + "/random_forest_model.sav")
    predictTestSet(DATADIR + "humdiv_dataset_imputed.tsv", script_path + "/random_forest_model.sav")
    '''

    ans = input("Run Random Forest for HumVar dataset? y/n")
    if ans == "y" or ans == "Y":
        answer = input("Execute preprocessing and train-test splitting?")
        if answer == "y" or answer == "Y":
            '''
            humvar_deleterious_dataset = encodeCategoricalFeatures(features_path + "humvar_deleterious_dataset.tsv")
            humvar_deleterious_dataset["output"] = 1
            humvar_deleterious_dataset.to_csv(script_path + "/data/humvar_deleterious_dataset.tsv", sep="\t", index=False)

            humvar_neutral_dataset = encodeCategoricalFeatures(features_path + "humvar_neutral_dataset.tsv")
            humvar_neutral_dataset["output"] = 0
            humvar_neutral_dataset.to_csv(script_path + "/data/humvar_neutral_dataset.tsv", sep="\t", index=False)
            
            preprocessing.mergeSets(DATADIR + "humvar_deleterious_dataset.tsv",
                      DATADIR + "humvar_neutral_dataset.tsv",
                      DATADIR + "humvar_dataset.tsv")

            merged_dataset_imputed = preprocessing.findBestK(DATADIR + "humvar_dataset.tsv", default_k=3)# k = 3
            merged_dataset_imputed.to_csv(DATADIR + "humvar_dataset_imputed.tsv", sep="\t", index=True)
            '''
        else:
            print("Will not run preprocessing!")

        answer = input("Do you want to evaluate the model? y/n")
        if answer == "y" or answer == "Y":
            predictTestSet(DATADIR + "humvar_dataset_imputed.tsv", script_path + "/random_forest_model.sav")
        else:
            print("Will not evaluate the model!")

    ans = input("Run Random Forest for HumDiv dataset? y/n")
    if ans == "y" or ans == "Y":
        answer = input("Execute preprocessing and train-test splitting?")
        if answer == "y" or answer == "Y":
            '''
            humdiv_deleterious_dataset = encodeCategoricalFeatures(features_path + "humdiv_deleterious_dataset.tsv")
            humdiv_deleterious_dataset["output"] = 1
            humdiv_deleterious_dataset.to_csv(script_path + "/data/humdiv_deleterious_dataset.tsv", sep="\t", index=False)

            humdiv_neutral_dataset = encodeCategoricalFeatures(features_path + "humdiv_neutral_dataset.tsv")
            humdiv_neutral_dataset["output"] = 0
            humdiv_neutral_dataset.to_csv(script_path + "/data/humdiv_neutral_dataset.tsv", sep="\t", index=False)

            preprocessing.mergeSets(DATADIR + "humdiv_deleterious_dataset.tsv",
                      DATADIR + "humdiv_neutral_dataset.tsv",
                      DATADIR + "humdiv_dataset.tsv")
            
            merged_dataset_imputed = preprocessing.findBestK(DATADIR + "humdiv_dataset.tsv", default_k=7)# k = 7
            merged_dataset_imputed.to_csv(DATADIR + "humdiv_dataset_imputed.tsv", sep="\t", index=True)
            '''
        else:
            print("Will not run preprocessing!")

        answer = input("Do you want to evaluate the model? y/n")
        if answer == "y" or answer == "Y":
            predictTestSet(DATADIR + "humdiv_dataset_imputed.tsv", script_path + "/random_forest_model.sav")
        else:
            print("Will not evaluate the model!")

    ans = input("Run Random Forest for our dataset? y/n")
    if ans == "y" or ans == "Y":
        answer = input("Execute preprocessing and train-test splitting?")
        if answer == "y" or answer == "Y":
            '''
            pathogenic_dataset = encodeCategoricalFeatures(features_path + "pathogenic_dataset.tsv")
            pathogenic_dataset["output"] = 1
            pathogenic_dataset.to_csv(DATADIR+"path_dataset.tsv", sep="\t", index=False)

            non_pathogenic_dataset = encodeCategoricalFeatures(features_path + "non_pathogenic_dataset.tsv")
            non_pathogenic_dataset["output"] = 0
            non_pathogenic_dataset.to_csv(DATADIR+"non_path_dataset.tsv", sep="\t", index=False)
            
            preprocessing.mergeSets(DATADIR + "path_dataset.tsv",
                      DATADIR + "non_path_dataset.tsv",
                      DATADIR + "merged_dataset.tsv")
                      
            merged_dataset_imputed = preprocessing.findBestK(DATADIR + "merged_dataset.tsv", default_k=1)  # k = 1
            merged_dataset_imputed.to_csv(DATADIR + "merged_dataset_imputed.tsv", sep="\t", index=True)

            preprocessing.splitTrainTestSets(DATADIR + "merged_dataset_imputed.tsv",
                                DATADIR + "train_set2.tsv",
                                DATADIR + "test_set2.tsv")
            '''
        else:
            print("Will not run preprocessing!")

        answer = input("Do you want to train the classifier? y/n")
        if answer == "y" or answer == "Y":
            trainRandomForest(train_set=DATADIR + "train_set2.tsv",
                              model_name=script_path + "/random_forest_model.sav")
        else:
            print("Will not train the classifier!")

        answer = input("Do you want to evaluate the model? y/n")
        if answer == "y" or answer == "Y":
            predictTestSet(DATADIR + "test_set2.tsv", script_path + "/random_forest_model.sav")
        else:
            print("Will not evaluate the model!")


if __name__ == "__main__":
    runPreprocessTrain()
