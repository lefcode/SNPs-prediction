import pandas as pd
import os
import plotly.graph_objs as go

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.svm import SVC
import itertools as it
from sklearn.metrics import make_scorer, accuracy_score, roc_auc_score

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
DATADIR = script_path + "/data/"

'''
def HyperParameterTuning(train_set):
    """Hyper parameters tuning using graphics card"""
    from tune_sklearn import TuneGridSearchCV
    # Set training and validation sets

    train_set = pd.read_csv(train_set, sep="\t", low_memory=False)

    train_output = train_set["output"].values
    train_features = train_set[train_set.columns.drop(["labels", "output"])].values

    X_train, X_test, y_train, y_test = train_test_split(train_features, train_output, test_size=0.20)

    parameters = {
        'n_estimators': [int(x) for x in range(100, 2000, 400)],
        'max_features': ['log2', 'sqrt'],
        'criterion': ["gini","entropy"],
    }

    tune_search = TuneGridSearchCV(
        RandomForestClassifier(),# n_jobs=-1
        parameters,
    )

    tune_search.fit(X_train, y_train)

    print(tune_search.best_params)
    print(tune_search.best_estimator_)

    pred = tune_search.predict(X_test)
    accuracy = np.count_nonzero(np.array(pred) == np.array(y_test)) / len(pred)
    print("Tune Accuracy:", accuracy)
'''

def tuneRandomForest(train_set):
    """Performs hyper-parameter tuning for Random Forest based on specific metrics
    Results indicate best combination of parameter values is:
    n_estimators = 2000, max_features = sqrt, criterion = gini
    """

    auc_score = make_scorer(roc_auc_score)
    acc = make_scorer(accuracy_score)

    train_set = pd.read_csv(train_set, sep="\t", low_memory=False)

    train_output = train_set["output"].values
    train_features = train_set[train_set.columns.drop(["labels", "output"])].values

    #X_train, X_test, y_train, y_test = train_test_split(train_features, train_output, test_size=0.20)

    # define parameters to be optimized
    parameters = {
        'n_estimators': [int(x) for x in range(200, 3000, 300)],
        'max_features': ['log2', 'sqrt', "auto"],
        'criterion': ["gini", "entropy"],
    }
    #plotGrid(parameters, script_path + "/results/GridSearchPlot.png")

    scores = ['precision', 'recall', 'f1', auc_score, acc] # compute efficiency based on scores
    for score in scores:
        print("# Tuning hyper-parameters for %s" % score)

        tune_search = GridSearchCV(
                RandomForestClassifier(n_jobs=-1),
                parameters,
                scoring=score
        )
        #tune_search.fit(X_train, y_train)
        tune_search.fit(train_features, train_output)
        print(tune_search.best_params_)

        means = tune_search.cv_results_['mean_test_score']
        stds = tune_search.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, tune_search.cv_results_['params']):
            print("%0.3f (+/-%0.03f) for %r" % (mean, std * 2, params))

        #y_true, y_pred = y_test, tune_search.predict(X_test)
       # print(classification_report(y_true, y_pred))
        #print()
    

def plotGrid(parameters, plot_name):
    """Creates a 3d plot of the grid for GridSearch"""

    allNames = sorted(parameters)
    combinations = it.product(*(parameters[Name] for Name in allNames))

    combs = [list(c) for c in list(combinations) ]
    combs_df = pd.DataFrame(combs, columns=["criterion", "max_features", "n_estimators"])

    trace = go.Scatter3d(
        x=combs_df['criterion'],
        y=combs_df['max_features'],
        z=combs_df['n_estimators'],
        mode='markers',
        marker=dict(
            size= 5,
            color='green',
            opacity=0.99,
            colorscale='Viridis',
            showscale=False
        )
    )
    data = [trace]
    layout = go.Layout(
        margin=dict(
            l=30,
            r=30,
            b=30,
            t=30
        ),
        height=700,
        width=960,
        showlegend=True,
        scene=dict(
            xaxis=dict(
                title='criterion',
                nticks=10,
            ),
            yaxis=dict(
                title='max_features',
                nticks=10,
            ),
            zaxis=dict(
                title='n_estimators',
                nticks=10,
            ),
            camera=dict(
                eye=dict(
                    y=2.089757339892154,
                    x=-0.5464711077183096,
                    z=0.14759264478960377,
                )
            ),
        ),
    )

    fig = go.Figure(data=data, layout=layout)
    fig.show()
    fig.write_image(plot_name)
    fig.write_html(script_path + "/results/GridSearchPlot.html")


def HyperSVM(train_set):

    from sklearn.metrics import make_scorer, accuracy_score, roc_auc_score
    auc_score = make_scorer(roc_auc_score)
    acc = make_scorer(accuracy_score)

    train_set = pd.read_csv(train_set, sep="\t", low_memory=False)

    train_output = train_set["output"].values
    train_features = train_set[train_set.columns.drop(["labels", "output"])].values

    X_train, X_test, y_train, y_test = train_test_split(train_features, train_output, test_size=0.20)

    # define parameters to be optimized
    parameters = {
        'kernel': ["poly", "rbf"],
        'C': [0.1, 1, 10, 100],
        'gamma': [0.01, 0.001, 0.0001]
    }

    scores = ['precision', 'recall', acc]  # compute efficiency based on scores
    for score in scores:
        print("# Tuning hyper-parameters for %s" % score)

        tune_search = GridSearchCV(
            SVC(), parameters, scoring=score, n_jobs=-1
        )

        tune_search.fit(X_train, y_train)
        print(tune_search.best_params_)

        means = tune_search.cv_results_['mean_test_score']
        stds = tune_search.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, tune_search.cv_results_['params']):
            print("%0.3f (+/-%0.03f) for %r" % (mean, std * 2, params))

        y_true, y_pred = y_test, tune_search.predict(X_test)
        print(classification_report(y_true, y_pred))
        print()


if __name__ == "__main__":
    # HyperParameterTuning(DATADIR + "train_set2.tsv")
    tuneRandomForest(DATADIR + "train_set2.tsv")
 #   HyperSVM(script_path + "/results/important_features_train_set.tsv")



