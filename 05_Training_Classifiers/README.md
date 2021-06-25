# Training and evaluation

In this part of the project uses the created datasets that emerged from the previous step to predict SNPs' pathogenity. <br>

Firstly, preprocessing takes place to
1. Encode categorical features
2. Impute missing values and scale data
3. Segregate a balanced training set and a completely distinct testing set to evaluate the predictions


Secondly, Hyper-parameter tuning on the parameters of Random Forest Classifier takes place to optimize the forest's ability to predict.
The resulting parameters are already defined so we use them to the next step.


Thirdly, 5-fold cross-validation training takes place on the Random Forest Classifier with the training set.
Also, we visualize the forest, find the most important features for the classification based on Gini Information Gain, plot the AUROC plot and save the model to be reusable to predict new data. 

Lastly, the classifier is evaluated based on several metrics and the results are stored to the results/ folder.

In order to execute the training procedure one must follow these steps:

```
cd 05_Training_Classifier/
python3 training.py
```

During the script's execution questions are asked to the user for preprocessing, training and testing
