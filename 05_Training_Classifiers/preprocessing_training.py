import pandas as pd
import numpy as np
import os
import errno
import random

from sklearn.feature_extraction import FeatureHasher
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import train_test_split, KFold, learning_curve
from sklearn.metrics import f1_score, mean_squared_error,accuracy_score, mean_absolute_error,\
	precision_score, recall_score, roc_auc_score, fbeta_score, roc_curve
from sklearn import svm

from matplotlib import pyplot
import warnings
warnings.filterwarnings('always')  # "error", "ignore", "always", "default", "module" or "once"

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
projectPath = os.path.abspath(os.path.join(script_path, os.pardir))
featuresPath = projectPath+"/04_Features_computation"

DATADIR = script_path + "/data/"

try:
	os.makedirs(DATADIR)
except OSError as e:
	if e.errno != errno.EEXIST: raise


def encodeCategoricalFeatures(dataset_file):
	""" Creates a list of dictionaries from the labels-values of every column that belongs to the
	categorical features. The dictionary is hashed with the FeatureHasher of sklearn and 2 new columns
	with hashed/encoded values is produced. These columns replace the single column of the feature. """

	print("Encoding categorical features...")
	dataset = pd.read_csv(dataset_file, sep="\t", low_memory=False)

	features_to_be_hashed = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T",
							 "V","W","Y","U","O"]
	#list=	["nA","nC","nD","nE","nF","nG","nH","nI","nK","nL","nM","nN","nP",
	#						 "nQ","nR","nS","nT","nV","nW","nY","nU","nO","nX","struct_interactions"]

	for feat in features_to_be_hashed:
		# find the correct n_features value. the value must be multiplied with 2
		# and give result >= #distinct feature values
		# e.x. 3 distinct values 2*n_features>=3 -> n_features=2
		# e.x. 10 distinct values 2*n_features>=10 -> n_features=4
		n_feat = 1
		while True:
			if 2 ** n_feat >= len(dataset[feat].value_counts()):
				break
			else:
				n_feat += 1

		n_feat=2
		feats_list = [feat+"_"+str(i+1) for i in range(n_feat)]

		hasher = FeatureHasher(n_features=n_feat, input_type="dict")
		column = dataset[["labels",feat]].set_index('labels')#.T.to_dict("records")

		# covert each index-column value to a dictionary and put them all in a list
		column_dict_list = list()
		for ind, c in zip(column.index, column.values):
			column_dict = dict({ind: int(c[0])})
			column_dict_list.append(column_dict)

		column_encoded = hasher.fit_transform(column_dict_list).toarray()
		dataset = dataset.drop(feat, axis=1)
		dataset[feats_list] = column_encoded

	print("Encoding done!")
	return dataset


def applyScaleKnn(dataset, neighbors, indices):
	""" """

	new_dataset = dataset.set_index("labels").astype(object).replace('None', np.nan)
	new_dataset = StandardScaler().fit_transform(new_dataset.to_numpy())

	imputer = KNNImputer(n_neighbors=neighbors)

	print("performing knn imputation on dataset!")
	imputed = imputer.fit_transform(new_dataset)
	cols = [x for x in dataset.columns if x != 'labels']
	imputed = pd.DataFrame(imputed, columns=cols, index=indices)
	print("knn imputation completed!")
	return imputed


def diminishSet(full_data):
	""" Randomly selects the 1% of rows of dataset and deletes the value from a randomly selected column"""

	full_data = pd.DataFrame(full_data)
	rowIndices = np.random.choice(full_data.shape[0], int(full_data.shape[0]/100), replace=False)

	nan_replacement_dict={}
	for i in rowIndices:
		col = random.randint(0, full_data.shape[1]-1)
		if not pd.isnull(full_data.loc[i][col]): # only if the value is known
		# keep old values to a dictionary with row,col as key and old values as values.
			nan_replacement_dict[str(i) + "," + str(col)] = full_data.loc[i][col]
			full_data.loc[i][col] = np.nan

	return full_data, nan_replacement_dict


def applyOptimizedKnn(data):
	"""
	Sets labels ase indices, replaces "None" with np.nan, scales with StandardScaler and
	performs Knn missing values Imputation.

	After using this function the outcome is that k=7 results on the minimum RMSE
	"""
	# assign labels as indexes and replace "None" with np.nan

	print("Applying Standard Scaling and Knn imputation!")

	new_dataset = data.set_index("labels").astype(object).replace('None', np.nan)
	scaled_dataset = StandardScaler().fit_transform(new_dataset.to_numpy())

	scaled_dataset, nan_replacement_dict = diminishSet(scaled_dataset)

	errors_dict = dict()
	for k in range(1, 21, 2):
		imputer = KNNImputer(n_neighbors=k)
		imputed = imputer.fit_transform(scaled_dataset)
		cols = [x for x in data.columns if x != 'labels']
		df_imputed = pd.DataFrame(imputed, columns=cols)  #, index=new_dataset.index

		print("knn imputation applied with K = " + str(k))

		predicted = list()
		for key, val in nan_replacement_dict.items():
			row, col = key.split(",")[0], key.split(",")[1]
			predicted.append(df_imputed.loc[int(row)][int(col)])

		mae = mean_absolute_error(predicted, list(nan_replacement_dict.values()))
		print({"K" : k, "Mean absolute error": mae})
		errors_dict[k] = mae

	minimum_error_k = min(errors_dict.keys(), key=(lambda k: errors_dict[k]))

	print('Minimum MAE with K:', minimum_error_k, "MAE:", errors_dict[minimum_error_k])
	print("New dataset created!")

	df_imputed = applyScaleKnn(data, neighbors=7, indices=data["labels"])

	return df_imputed


def diminishNonPathogenicSet(pathogenic_dataset, non_pathogenic_dataset):
	"""
	We keep the 19.674 records of the non-pathogenic dataset in order to have same
	lengths in both datasets
	"""

	print("Splitting non pathogenic dataset to train and test sets.")

	pathogenic_dataset = pd.read_csv(pathogenic_dataset, sep="\t", low_memory=False)
	non_pathogenic_dataset = pd.read_csv(non_pathogenic_dataset, sep="\t", low_memory=False)

	rowIndices = np.random.choice(len(non_pathogenic_dataset), len(pathogenic_dataset), replace=False)
	# randomly select 19.674 indices in order to get their records

	# keep the rest non patogenic rows for the test set
	all_indices = [i for i in range(len(non_pathogenic_dataset))]

	for ind in all_indices:
		if ind in rowIndices:
			all_indices.remove(ind) # keep only indices not existing in the rowIndices

	# create a list of lists with all of the rest records to put them into a file
	rest_records = list()
	for i in all_indices:
		values_list = non_pathogenic_dataset.loc[i].tolist()
		record = ""
		for i in range(len(values_list)):
			if i == len(values_list)-1:
				record += str(values_list[i])
			else:
				record += str(values_list[i]) + "\t"

		rest_records.append(record.split("\t"))

	df = pd.DataFrame(rest_records, columns=non_pathogenic_dataset.columns)
	df["output"] = 0
	df.to_csv(DATADIR + "/non_pathogenic_test_set.tsv", sep="\t", index=False)

	# create a list of lists with all the records to put them into a file
	frame_lists = list()
	for i in rowIndices:
		values_list = non_pathogenic_dataset.loc[i].tolist()
		record = ""
		for i in range(len(values_list)):
			if i == len(values_list)-1:
				record += str(values_list[i])
			else:
				record += str(values_list[i]) + "\t"

		frame_lists.append(record.split("\t"))

	df = pd.DataFrame(frame_lists, columns= non_pathogenic_dataset.columns)
	df.to_csv(DATADIR+"/non_pathogenic_train_set.tsv", sep="\t", index=False)

	print("Splitting done!")

	return DATADIR+"/non_pathogenic_train_set.tsv", DATADIR + "/non_pathogenic_test_set.tsv"


def mergeSets(pathogenic_dataset, non_pathogenic_dataset):
	"""
	We keep the 2/3 of both pathogenic and non-pathogenic datasets as the train set
	and the 1/3 as the test set.
	"""

	print("Merging two datasets into merged_dataset.tsv")

	pathogenic_dataset = pd.read_csv(pathogenic_dataset, sep="\t", low_memory=False, index_col=0)
	pathogenic_dataset["output"] = 1
	#print(pathogenic_dataset)

	non_pathogenic_dataset = pd.read_csv(non_pathogenic_dataset, sep="\t", low_memory=False,  index_col=0)
	non_pathogenic_dataset["output"] = 0
	#print(non_pathogenic_dataset)
	df_complete = pd.concat([pathogenic_dataset, non_pathogenic_dataset])
	df_complete.to_csv(DATADIR+"/merged_dataset.tsv",sep="\t", index=True)


def mergeTestSets(x_test_set, y_test_set, non_pathogenic_test_set):
	""" Merges the test set of the merged file with the rest of the records of the non pathogenic dataset
	which are stored in non_pathogenic_test_set.tsv """

	print("Merging test set with the rest non pathogenic records")

	non_pathogenic_test_set = pd.read_csv(non_pathogenic_test_set, sep="\t", low_memory=False).drop(["labels"],axis=1)

	x = non_pathogenic_test_set.drop(["output"], axis=1)#[["DNAbindingMut","deltaG"]].values
	y = non_pathogenic_test_set["output"].values

	x_test_set = np.concatenate((x_test_set, x), axis=0)
	y_test_set = np.concatenate((y_test_set, y), axis=0)

	return x_test_set, y_test_set


def testFunc(merged_dataset, non_pathogenic_test_set):

	from scipy import stats

	merged_dataset = pd.read_csv(merged_dataset, sep="\t", low_memory=False).drop(["labels"], axis=1)
	columns = merged_dataset.columns

	for c in columns:
		print(stats.chisquare(merged_dataset[c], merged_dataset["output"]))
		#print(stats.ttest_ind(merged_dataset[c], merged_dataset["output"]))


	labels = merged_dataset["output"].values
	features = merged_dataset[columns.drop(["output"])].values

	kf = KFold(n_splits=5, shuffle=True, random_state=42)  # Use for 5-Fold classification

	rf = RandomForestClassifier(n_estimators=100, criterion="gini", random_state=42,
								max_features="auto", bootstrap=True, n_jobs=-1)


	X_train, X_test, y_train, y_test  = train_test_split(features, labels, test_size=0.45)

	X_test, y_test = mergeTestSets(X_test, y_test, non_pathogenic_test_set)



def trainRandomForest(merged_dataset, non_pathogenic_test_set):
	""" """

	print("Training random forest classifier")

	merged_dataset = pd.read_csv(merged_dataset, sep="\t", low_memory=False).drop(["labels"], axis=1)

	columns = merged_dataset.columns
	labels = merged_dataset["output"].values

	features = merged_dataset[columns.drop(["output"])].values


	kf = KFold(n_splits=5, shuffle=True, random_state=42) # Use for 5-Fold classification

	rf = RandomForestClassifier(n_estimators=100, criterion="gini", random_state=42,
								max_features="auto", bootstrap=True, n_jobs=-1)

	classifier = svm.SVC(kernel='rbf', decision_function_shape='ovo', cache_size=20000)

	mae_list = list()
	precision_list = list()
	recall_list = list()
	f1_score_list = list()
	area_under_curve_list = list()
	accuracy_list = list()
	f2_score_list = list()


	for train_index, test_index in kf.split(features):

		X_train, X_test, y_train, y_test = features[train_index], features[test_index], \
										   labels[train_index], labels[test_index]

		'''
		print('Training Features Shape:', X_train.shape)
		print('Training Labels Shape:', y_train.shape)
		print('Testing Features Shape:', X_test.shape)
		print('Testing Labels Shape:', y_test.shape)
		'''

		X_test, y_test = mergeTestSets(X_test, y_test, non_pathogenic_test_set)

		rf.fit(X_train,y_train)
		prediction = rf.predict(X_test)

		#classifier.fit(X_train,y_train)
		#prediction = classifier.predict(X_test)

		mae, precision, recall, f1, area_under_curve, accuracy, f2 = countMetrics(prediction, y_test)

		mae_list.append(mae)
		precision_list.append(precision)
		recall_list.append(recall)
		f1_score_list.append(f1)
		area_under_curve_list.append(area_under_curve)
		accuracy_list.append(accuracy)
		f2_score_list.append(f2)

		'''	
		lr_fpr, lr_tpr, _ = roc_curve(y_test, prediction.round())
		# plot the roc curve for the model
		pyplot.plot(lr_fpr, lr_tpr, marker='.', label='RandomForest')
		# axis labels
		pyplot.xlabel('False Positive Rate')
		pyplot.ylabel('True Positive Rate')
		# show the legend
		pyplot.legend()
		# show the plot
		pyplot.show()

		'''

		pyplot.scatter(prediction, y_test)
		pyplot.show()

	print("Average Mean absolute error (MAE): %f" % np.mean(mae_list))
	print("Average Precision Score : %f" % np.mean(precision_list))
	print("Average Recall Score : %f" % np.mean(recall_list))
	print("Average f1 Score : %f " %  np.mean(f1_score_list))
	print("Average AUC :  %f " % np.mean(area_under_curve_list))
	print("Average Accuracy :  %f " %  np.mean(accuracy_list))
	print("Average F2 score :  %f" %  np.mean(f2_score_list))



def countMetrics(prediction, y_test):

	rmse = mean_squared_error(y_test, list(map(int, prediction)))
	print("Root mean squared error (RMSE): %f" % rmse)
	mae = mean_absolute_error(y_test, list(map(int, prediction)))
	print("Mean absolute error (MAE): %f" % mae)
	precision = precision_score(y_test, list(map(int, prediction)), average='weighted', labels=np.unique(prediction))
	print("Precision Score : ", precision)
	recall = recall_score(y_test, list(map(int, prediction)), average='weighted', labels=np.unique(prediction))
	print("Recall Score : ", recall)
	f1 = f1_score(y_test, list(map(int, prediction)), average='weighted', labels=np.unique(prediction))
	print("f1 Score : ", f1)
	area_under_curve = roc_auc_score(y_test, list(map(int, prediction)))
	print("AUC : ", area_under_curve)
	accuracy = accuracy_score(y_test, list(map(int, prediction)))
	print("Accuracy : ", accuracy)
	f2 = fbeta_score(y_test, list(map(int, prediction)), average='weighted', beta=0.5)
	print("F2 score : ", f2)
	# print(precision_recall_fscore_support(y_test, prediction.round(), average='weighted'))

	return mae, precision, recall, f1, area_under_curve, accuracy, f2


def learningCurve(merged_dataset):

	merged_dataset = pd.read_csv(merged_dataset, sep="\t", low_memory=False).drop(["labels"], axis=1)

	columns = merged_dataset.columns
	labels = merged_dataset["output"].values
	features = merged_dataset[columns.drop(["output"])].values

	train_sizes, train_scores, test_scores = learning_curve(RandomForestClassifier(),
                                                        features,
                                                        labels,
                                                        # Number of folds in cross-validation
                                                        cv=5,
                                                        # Evaluation metric
                                                        scoring='accuracy',
                                                        # Use all computer cores
                                                        n_jobs=-1,
                                                        # 50 different sizes of the training set
                                                        train_sizes=np.linspace(0.01, 1.0, 50))

	# Create means and standard deviations of training set scores
	train_mean = np.mean(train_scores, axis=1)
	train_std = np.std(train_scores, axis=1)

	# Create means and standard deviations of test set scores
	test_mean = np.mean(test_scores, axis=1)
	test_std = np.std(test_scores, axis=1)

	# Draw lines
	pyplot.subplots(1, figsize=(10, 10))
	pyplot.plot(train_sizes, train_mean, '--', color="#111111", label="Training score")
	pyplot.plot(train_sizes, test_mean, color="#111111", label="Cross-validation score")

	pyplot.fill_between(train_sizes, train_mean - train_std, train_mean + train_std, color="#DDDDDD")
	pyplot.fill_between(train_sizes, test_mean - test_std, test_mean + test_std, color="#DDDDDD")

	pyplot.title("Learning Curve")
	pyplot.xlabel("Training Set Size"), pyplot.ylabel("Accuracy Score"), pyplot.legend(loc="best")
	pyplot.tight_layout()
	pyplot.show()


if __name__=="__main__":
	'''
	pathogenic_dataset_file = featuresPath +"/pathogenic_dataset.tsv"
	non_pathogenic_dataset_file = featuresPath + "/non_pathogenic_dataset.tsv"

	pathogenic_dataset = encodeCategoricalFeatures(pathogenic_dataset_file)
	non_pathogenic_dataset = encodeCategoricalFeatures(non_pathogenic_dataset_file)

	pathogenic_imputed = applyOptimizedKnn(pathogenic_dataset)
	pathogenic_imputed.to_csv(script_path + "/data/pathogenic_dataset.tsv", sep="\t")

	non_pathogenic_imputed = applyOptimizedKnn(non_pathogenic_dataset)
	non_pathogenic_imputed.to_csv(script_path + "/data/non_pathogenic_dataset.tsv", sep="\t")

	non_pathogenic_train_set, non_pathogenic_test_set = diminishNonPathogenicSet(DATADIR+"/pathogenic_dataset.tsv",
														DATADIR+"/non_pathogenic_dataset.tsv")

	mergeSets(DATADIR+"/pathogenic_dataset.tsv", DATADIR+"/non_pathogenic_train_set.tsv")

	'''

	#trainRandomForest(DATADIR + "/merged_dataset.tsv", DATADIR + "/non_pathogenic_test_set.tsv")
	testFunc(DATADIR + "/merged_dataset.tsv", DATADIR + "/non_pathogenic_test_set.tsv")
