import pandas as pd
import os


def findDifference(old_set, new_set, type_dat):

    print(type_dat)
    script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
    df_old = pd.read_csv(script_path + old_set, sep="\t", low_memory=False)
    df_new = pd.read_csv(script_path + new_set, sep="\t", low_memory=False)

    old_labels = df_old["labels"]
    new_labels = df_new["labels"]

    for old in old_labels.values:
        if old not in new_labels.values:
            print(old)

    for new in new_labels.values:
        if new not in old_labels.values:
            print(new)

findDifference("/old_features/humvar_deleterious_dataset.tsv", "/humvar_deleterious_dataset.tsv", "var del")
findDifference("/old_features/humdiv_deleterious_dataset.tsv", "/humdiv_deleterious_dataset.tsv", "div del")
findDifference("/old_features/humvar_neutral_dataset.tsv", "/humvar_neutral_dataset.tsv", "var neut")
findDifference("/old_features/humdiv_neutral_dataset.tsv", "/humdiv_neutral_dataset.tsv", "div neut")