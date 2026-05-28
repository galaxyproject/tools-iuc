import argparse
import copy
import itertools
import warnings
import os

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from imblearn.over_sampling import SMOTE, RandomOverSampler
from imblearn.pipeline import Pipeline
from imblearn.under_sampling import NearMiss, RandomUnderSampler
from lightgbm import LGBMClassifier
from matplotlib.colors import ListedColormap
from sklearn import datasets
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.datasets import load_iris
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.feature_selection import (
    SelectFromModel,
    SelectKBest,
    VarianceThreshold,
    f_classif,
)
from sklearn.metrics import (
    auc,
    matthews_corrcoef,
    precision_recall_fscore_support,
    roc_auc_score,
    roc_curve,
    precision_recall_curve,
    f1_score,
    precision_score,
    recall_score,
)  
    
from sklearn.model_selection import (
    GridSearchCV,
    KFold,
    StratifiedKFold,
    cross_val_predict,
    cross_val_score,
    train_test_split,
)
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import LabelBinarizer, LabelEncoder, label_binarize, StandardScaler
from sklearn.svm import SVC

#from tabpfn import TabPFNClassifier
# from tabpfn_extensions.post_hoc_ensembles.sklearn_interface import AutoTabPFNClassifier
from xgboost import XGBClassifier
from itertools import combinations

# Suppress FutureWarnings
warnings.filterwarnings("ignore", category=FutureWarning)



def split_classes(X, y):
    return {
        (c1, c2): (X[(y == c1) | (y == c2)], y[(y == c1) | (y == c2)])
        for c1, c2 in itertools.combinations(np.unique(y), 2)
    }


def ovo_and_ova_multiclass_auc(X, y, base_clf, p_grid, random_state, model_name):
    results = {}
    plot_data = []
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    class_names = le.classes_

    # Stratified K-Folds
    inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)

    ####################
    # One-vs-Rest Classification
    ####################
    print("Performing One vs Rest classification")

    #checking grid search enabled or not 
    if p_grid is not None:
        ovr_clf = GridSearchCV(
            estimator=OneVsRestClassifier(base_clf),
            param_grid=p_grid,
            cv=inner_cv,
            scoring="roc_auc_ovr",
        )
    else:
        ovr_clf = OneVsRestClassifier(base_clf)
   
    
    y_score = cross_val_predict(ovr_clf, X, y_encoded, cv=outer_cv, method="predict_proba") 
    y_pred = np.argmax(y_score, axis=1) 
    
    # Per-class metrics for OvR 
    per_class_precision = [] 
    per_class_recall = [] 
    per_class_f1 = [] 
    per_class_mcc = []

    for idx, cls in enumerate(class_names): 
        y_bin = (y_encoded == idx).astype(int) 
        cls_score = y_score[:, idx] 
        
        # Ensure minority class is positive 
        if np.sum(y_bin) > np.sum(1 - y_bin): 
            y_bin = 1 - y_bin 
            cls_score = 1 - cls_score 
            
        y_pred_bin = (y_pred == idx).astype(int) 
        precision, recall, f1, _ = precision_recall_fscore_support(y_bin, y_pred_bin, average="binary") 
        mcc = matthews_corrcoef(y_bin, y_pred_bin) 
        prec_curve, rec_curve, _ = precision_recall_curve(y_bin, cls_score) 
        pr_auc_val = auc(rec_curve, prec_curve) 
        roc_auc_val = roc_auc_score(y_bin, cls_score) 
        
        results[f"{cls} vs Rest - Precision"] = precision 
        results[f"{cls} vs Rest - Recall"] = recall 
        results[f"{cls} vs Rest - F1"] = f1 
        results[f"{cls} vs Rest - MCC"] = mcc 
        results[f"{cls} vs Rest - PR AUC"] = pr_auc_val 
        results[f"{cls} vs Rest - ROC AUC"] = roc_auc_val 
        
        per_class_precision.append(precision) 
        per_class_recall.append(recall) 
        per_class_f1.append(f1) 
        per_class_mcc.append(mcc) 
    
    # Macro metrics OvR 

    macro_ovr_auc = np.mean([results[f"{cls} vs Rest - ROC AUC"] for cls in class_names]) 
    macro_ovr_precision = np.mean(per_class_precision) 
    macro_ovr_recall = np.mean(per_class_recall) 
    macro_ovr_f1 = np.mean(per_class_f1) 
    macro_ovr_mcc = np.mean(per_class_mcc)
    macro_ovr_pr_auc = np.mean([results[f"{cls} vs Rest - PR AUC"] for cls in class_names])

    
    results["OvR Macro ROC AUC"] = macro_ovr_auc
    results["OvR Macro Precision"] = macro_ovr_precision
    results["OvR Macro Recall"] = macro_ovr_recall
    results["OvR Macro F1"] = macro_ovr_f1 
    results["OvR Macro MCC"] =  macro_ovr_mcc
    results["OvR Macro PR AUC"] = macro_ovr_pr_auc 

    '''
    print(f"Macro ROC AUC (OvR): {macro_ovr_auc:.4f}")
    print(f"Macro Precision (OvR): {macro_ovr_precision:.4f}")
    print(f"Macro Recall (OvR): {macro_ovr_recall:.4f}")
    print(f"Macro F1 (OvR): {macro_ovr_f1:.4f}")
    print(f"Macro MCC (OvR): {macro_ovr_mcc:.4f}")
    print(f"Macro PR AUC (OvR): {macro_ovr_pr_auc:.4f}")  '''

    

    #avoiding  meaningless computation as OvO metrics won’t make sense with TabPFN
    if model_name == "tabpfn":
        print("Skipping One-vs-One metrics for TabPFN")
    else:
        
        ####################
        # One-vs-One Classification
        ####################
        print("Performing One vs One classification")
    
        ovo_auc = {} 
        ovo_precision = {} 
        ovo_recall = {} 
        ovo_f1 = {} 
        ovo_mcc = {} 
        
        for c1, c2 in combinations(range(len(class_names)), 2): 
            mask = np.isin(y_encoded, [c1, c2]) 
            X_pair, y_pair = X[mask], y_encoded[mask] 
    
            # checking grid search enabled or not
            if p_grid is not None:
                ovo_clf = GridSearchCV(
                    estimator=base_clf,
                    param_grid={k.replace("estimator__", ""): v for k, v in p_grid.items()},
                    cv=inner_cv,
                    scoring="roc_auc"
                )
            else:
                ovo_clf = base_clf
                
            y_score_pair = cross_val_predict(ovo_clf, X_pair, y_pair, cv=outer_cv, method="predict_proba") 
            
            # Identify minority 
            
            vals, counts = np.unique(y_pair, return_counts=True) 
            minority = vals[np.argmin(counts)] 
            minority_idx = np.where([c1, c2] == minority)[0][0] 
            
            y_bin = (y_pair == minority).astype(int) 
            y_score_cls = y_score_pair[:, minority_idx] 
            
            
            
            # Ensure minority positive 
            
            if np.sum(y_bin) > np.sum(1 - y_bin): 
                y_bin = 1 - y_bin 
                y_score_cls = 1 - y_score_cls 
                
            y_pred_bin = (np.argmax(y_score_pair, axis=1) == minority_idx).astype(int)
                                                                                 
            precision, recall, f1, _ = precision_recall_fscore_support(y_bin, y_pred_bin, average="binary") 
            mcc = matthews_corrcoef(y_bin, y_pred_bin) 
            prec_curve, rec_curve, _ = precision_recall_curve(y_bin, y_score_cls) 
            pr_auc_val = auc(rec_curve, prec_curve) 
            roc_auc_val = roc_auc_score(y_bin, y_score_cls)
            
            pair_name = f"{le.inverse_transform([c1])[0]} vs {le.inverse_transform([c2])[0]}" 
            
            results[f"{pair_name} - Precision"] = precision 
            results[f"{pair_name} - Recall"] = recall 
            results[f"{pair_name} - F1"] = f1 
            results[f"{pair_name} - MCC"] = mcc 
            results[f"{pair_name} - PR AUC"] = pr_auc_val 
            results[f"{pair_name} - ROC AUC"] = roc_auc_val 
            
            ovo_auc[(c1, c2)] = roc_auc_val 
            ovo_precision[(c1, c2)] = precision 
            ovo_recall[(c1, c2)] = recall 
            ovo_f1[(c1, c2)] = f1 
            ovo_mcc[(c1, c2)] = mcc 
            
            # for plotting 
            plot_data.append({
                "class_a": le.inverse_transform([c1])[0],
                "class_b": le.inverse_transform([c2])[0],
                "pair_name": pair_name,
                "y_true": y_bin.copy(),
                "y_prob": y_score_cls.copy(),
                "roc_auc": roc_auc_val,
                "pr_auc": pr_auc_val
            })
            
        # Macro metrics OvO 
        macro_ovo_auc = np.mean(list(ovo_auc.values()))
        macro_ovo_precision = np.mean(list(ovo_precision.values()))
        macro_ovo_recall = np.mean(list(ovo_recall.values())) 
        macro_ovo_f1 = np.mean(list(ovo_f1.values())) 
        macro_ovo_mcc = np.mean(list(ovo_mcc.values())) 
        macro_ovo_pr_auc = np.mean([results[k] for k in results if "vs" in k and "PR AUC" in k]) 
    
        results["OvO Macro ROC AUC"] =  macro_ovo_auc
        results["OvO Macro Precision"] = macro_ovo_precision
        results["OvO Macro Recall"] = macro_ovo_recall
        results["OvO Macro F1"] = macro_ovo_f1
        results["OvO Macro MCC"] = macro_ovo_mcc
        results["OvO Macro PR AUC"] =  macro_ovo_pr_auc
    
    
        ''' 
        print(f"Macro ROC AUC (OvO): {macro_ovo_auc:.4f}")
        print(f"Macro Precision (OvO): {macro_ovo_precision:.4f}")
        print(f"Macro Recall (OvO): {macro_ovo_recall:.4f}")
        print(f"Macro F1 (OvO): {macro_ovo_f1:.4f}")
        print(f"Macro MCC (OvO): {macro_ovo_mcc:.4f}")
        print(f"Macro PR AUC (OvO): {macro_ovo_pr_auc:.4f}") '''
    
    return results, plot_data


def repeat_clf(n_seeds, ks, X, y, label, model, sampling_strategy, use_grid=False):

    print("features(ks): ", ks)
    print("seeds: ", n_seeds)

    # Define sampling strategies
    sampling_strategies = {
        "No Sampling": None,
        "Random OverSampling": RandomOverSampler(random_state=42),
        "SMOTE": SMOTE(random_state=42),
        "Random UnderSampling": RandomUnderSampler(random_state=42),
        "NearMiss (v1)": NearMiss(version=1),
        "NearMiss (v2)": NearMiss(version=2),
        "NearMiss (v3)": NearMiss(version=3),
    }

    # If the selected strategy is not in the dictionary, use "No Sampling"
    sampler = sampling_strategies.get(sampling_strategy, None)

    seed_results = {}

    for seed in range(n_seeds):

        ks_results = {}
        for k in ks:

            print(f"CV for seed {seed} and {k} features")

            # Create a Random Forest Classifier
            rf = RandomForestClassifier(random_state=seed)

            # Create a SelectFromModel using the Random Forest Classifier
            selector = SelectFromModel(rf, max_features=k)

            if model == "rf":
                ml_model = rf
                ml_model_grid = {
                    "estimator__classification__n_estimators":[100, 300, 500],  # Number of trees in the forest
                    "estimator__classification__max_depth": [None, 10, 20, 30],  # tree depth
                    "estimator__classification__max_features": ["sqrt", "log2"],  # Feature selection strategy
                    "estimator__classification__criterion": ["entropy"],  # Split criterion
                    "estimator__classification__min_samples_leaf": [1, 2, 4],  # Minimum samples per leaf
                }
            elif model == "xgb":
                ml_model = XGBClassifier(
                    use_label_encoder=False, eval_metric="logloss", random_state=seed
                )
                ml_model_grid = {
                    "estimator__classification__n_estimators": [100, 300, 500], 
                    "estimator__classification__gamma": [0, 0.1, 0.3], # min loss reduction
                    "estimator__classification__max_depth": [3, 5, 7], 
                    "estimator__classification__learning_rate": [0.01, 0.05, 0.1], # step size
                }
            elif model == "etc":
                ml_model = ExtraTreesClassifier(random_state=seed)
                ml_model_grid = {
                    "estimator__classification__n_estimators": [100, 300, 500],
                    "estimator__classification__max_depth": [None, 10, 20],       # tree depth
                    "estimator__classification__max_features": ["sqrt", "log2"],  # features per split
                    "estimator__classification__min_samples_leaf": [1, 2, 4],     # min leaf samples
                    
                }
            elif model == "lgbm":
                ml_model = LGBMClassifier(random_state=seed, verbose=-1)
                ml_model_grid = {
                    "estimator__classification__n_estimators": [100, 300, 500],  
                    "estimator__classification__learning_rate": [0.01, 0.05, 0.1],
                    "estimator__classification__num_leaves": [31, 63, 127],      # leaves per tree
                            
                }
            elif model == "tabpfn":
                from tabpfn import TabPFNClassifier
                ml_model = TabPFNClassifier(
                    device="cpu",  
                    n_estimators=32  # default
                )
                ml_model_grid = None  # TabPFN does not use GridSearch

            # If there is a sampler, include it in the pipeline
            steps = []
            if sampler:
                steps.append(("sampling", sampler))
            steps.append(("feature_selection", selector))
            steps.append(("classification", ml_model))

            # Create a pipeline with feature selection, sampling, and classification
            pipeline = Pipeline(steps=steps)

            ###########################

            # Run the classification with the sampling strategy
            if use_grid:
                results, plot_data = ovo_and_ova_multiclass_auc(
                    X, y, pipeline, ml_model_grid, random_state=seed ,  model_name=model
                )
            else:
                results , plot_data = ovo_and_ova_multiclass_auc(
                    X, y, pipeline, None, random_state=seed,  model_name=model
                )
                       

            # print(results)

            ks_results[k] = {
                "results": results,
                "plot_data": plot_data,
                "Label": label,
                "Model": model,
                "Sampling_Strategy": sampling_strategy,
            }


        seed_results[seed] = copy.copy(ks_results)

    return seed_results


def store_results(seed_results, output):

    # Flatten the nested dictionary into a DataFrame
    '''df = pd.DataFrame(
        {
            (outer_key, inner_key): values
            for outer_key, inner_dict in seed_results.items()
            for inner_key, values in inner_dict.items()
        }
    ).T

    # '''
    
    final_results = []
    metrics = ["ROC AUC", "Precision", "Recall", "F1", "MCC", "PR AUC"]
    
    for seed, ks_results in seed_results.items():
        for k, result_info in ks_results.items():
            result = result_info["results"]
            model = result_info["Model"]
            sampling_strategy = result_info["Sampling_Strategy"]
            label=result_info["Label"]
                        
            # Determine Class and Type
           
            groups = set()
            for key in result.keys():
                if "Macro" in key:
                    groups.add(("Macro", "OvR" if "OvR" in key else "OvO"))
                elif "vs Rest" in key:
                    groups.add((key.split(" vs Rest")[0], "OvR"))
                else:
                    groups.add((key.split(" - ")[0], "OvO"))

            # assign metric values according to class and type
            for class_name, type_name in groups:

                metric_values = {}
            
                for metric in metrics:
                    metric_key = None
            
                    for key in result.keys():
                        if class_name == "Macro":
                            if metric in key and "Macro" in key and type_name in key:
                                metric_key = key
                                break
                        elif type_name == "OvR":
                            if metric in key and f"{class_name} vs Rest" in key:
                                metric_key = key
                                break
                        else:  # OvO
                            if metric in key and "vs" in key and class_name in key:
                                metric_key = key
                                break
            
                    metric_values[metric] = result[metric_key] if metric_key else np.nan

                final_results.append({
                    "Seed": seed,
                    "Features (k)": k,
                    "Label": label,
                    "Model": model,
                    "Sampling_Strategy": sampling_strategy,
                    "Class": class_name,
                    "Type": type_name,
                    **metric_values
                })
                
    df = pd.DataFrame(final_results)

    '''#Set multi-level index names for clarity
    df.set_index(["Seed", "Features (k)", "Label", "Model", "Sampling_Strategy"], inplace=True)

    df.index.names = ["Seed", "Features (k)","Label","Model","Sampling_Strategy"]
    # Display the DataFrame
    df = df.reset_index()'''
    print(df.head())
    print(df.shape)

    df.to_csv(output, index=False)

    print(df)


def run_classification(X, y, ks, n_seeds,output, label,model, sampling_strategy,use_grid=False):

    '''# Ensure ks does not exceed the number of columns in X
    max_features = len(X.columns)
    ks = [k for k in ks if k <= max_features]
    if max_features not in ks:
        ks.append(max_features)'''

    seed_results = repeat_clf(n_seeds, ks, X, y, label,model, sampling_strategy, use_grid=use_grid)
    store_results(seed_results, output)
    
    return seed_results




        
def plot_pairwise_diagnostics(plot_data, diagnostic_plot):

    n_pairs = len(plot_data)

    fig, axes = plt.subplots(n_pairs, 3, figsize=(18, 4*n_pairs))

    if n_pairs == 1:
        axes = np.array([axes])

    for ax_row, item in zip(axes, plot_data):
        class_a = item["class_a"]
        class_b = item["class_b"]
        y_true = item["y_true"]
        y_prob = item["y_prob"]
        pair_name = item["pair_name"]

        ax_roc, ax_pr, ax_hist = ax_row

        # ROC
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        roc_auc = roc_auc_score(y_true, y_prob)

        ax_roc.plot(fpr, tpr, lw=2, label=f"AUC = {roc_auc:.2f}")
        ax_roc.plot([0, 1], [0, 1], 'k--', alpha=0.5)
        ax_roc.legend()
        ax_roc.set_title(f" {class_a} vs {class_b}\nROC Curve")
        ax_roc.set_xlabel("False Positive Rate")
        ax_roc.set_ylabel("True Positive Rate")

        # PR
        # PR curve
        prec, rec, _ = precision_recall_curve(y_true, y_prob)

        ax_pr.plot(
            rec,
            prec,
            lw=2,
            label="PR Curve"
        )

        # baseline
        baseline = np.mean(y_true)
        ax_pr.plot([0, 1], [baseline, baseline], 'k--', alpha=0.5)

        ax_pr.set_xlabel("Recall")
        ax_pr.set_ylabel("Precision")
        ax_pr.set_title(f"{class_a} vs {class_b}\nPrecision–Recall Curve")
        ax_pr.legend()

        # Histogram
        ax_hist.hist(
            [y_prob[y_true == 0], y_prob[y_true == 1]],
            bins=20,
            alpha=0.7,
            label=[f"Negative ({class_a})", f"Positive ({class_b})"]
            
        )

        ax_hist.set_title(f"{class_a} vs {class_b}\nPredicted Probabilities")
        ax_hist.set_xlabel("Predicted Probability")
        ax_hist.set_ylabel("Count")
        ax_hist.legend()

    fig.tight_layout()
    fig.savefig(diagnostic_plot, format="png")

def plot_model_performance_by_features(result_file, plot_per_feature):

    df = pd.read_csv(result_file, sep=",", engine="python")
    df.columns = df.columns.str.strip()
    
    # non macro classes 
    if "Class" in df.columns:
        df = df[df["Class"] != "Macro"]

    df["Features (k)"] = pd.to_numeric(df["Features (k)"], errors="coerce")

    metrics = ["ROC AUC", "PR AUC"]
    feature_sets = sorted(df["Features (k)"].dropna().unique())

    fig, axes = plt.subplots(
        len(feature_sets), 2,
        figsize=(14, 5 * len(feature_sets)),
        sharey=True
    )

    if len(feature_sets) == 1:
        axes = [axes]

    for i, feature in enumerate(feature_sets):

        for j, t in enumerate(["OvR", "OvO"]):

            ax = axes[i][j] if len(feature_sets) > 1 else axes[j]

            df_plot = df[(df["Features (k)"] == feature) & (df["Type"] == t)]

            
            df_melt = df_plot.melt(
                id_vars=["Class"],
                value_vars=metrics,
                var_name="Metric",
                value_name="Score"
            )

            sns.barplot(
                data=df_melt,
                x="Class",
                y="Score",
                hue="Metric",  
                edgecolor="black",
                errorbar="sd",
                ax=ax
            )

            for container in ax.containers:
                for bar in container:
                    h = bar.get_height()
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        h / 2,
                        f"{h:.2f}",
                        ha="center",
                        va="center",
                        fontsize=8
                    )

            ax.set_title(f"{feature} features - {t}")
            ax.set_xlabel("Class")
            ax.set_ylim(0, 1)
            ax.tick_params(axis="x", rotation=30)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", title="Metric")

    fig.suptitle("Model Performance Across Feature Sets", fontsize=14)

    fig.tight_layout(rect=[0, 0, 0.9, 0.95])

    fig.savefig(plot_per_feature, format="png", bbox_inches="tight")
    plt.close(fig)

            
           
      

def main():
    parser = argparse.ArgumentParser(description="Run Classification Model")

    parser.add_argument("--X", type=str, required=True, help="path to X")
    parser.add_argument("--y", type=str, required=True, help="path to y")
    parser.add_argument("--ks", type=str, required=True, help="list of values of k")
    parser.add_argument("--n_seeds", type=int, default=2, help="number of seeds")
    parser.add_argument("--label", type=str, required=True, help="add label for clarity")
    parser.add_argument("--model", type=str, required=True, help="choose model :['rf', 'XGB', 'ETC', 'lgbm', 'TabPFN']")
    parser.add_argument("--sampling_strategy", type=str, required=True, help="choose sampling strategy: ['No Sampling','Random OverSampling','SMOTE','Random UnderSampling','NearMiss (v1)','NearMiss (v2)','NearMiss (v3)']")
    parser.add_argument("--grid_search", action="store_true", help="grid search")
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--diagnostic_plot", type=str, required=True)
    parser.add_argument("--plot_per_feature", type=str, required=True)
    args = parser.parse_args()
    
    
    
    # reading str file paths
    X = pd.read_csv(args.X, sep="\t",index_col=0)
    y = pd.read_csv(args.y, sep="\t",index_col=0)
    
    # keep only numeric columns for sampling strategies
    X = X.select_dtypes(include=[np.number])
   
    # ks str value to int list 
    ks = [int(x.strip()) for x in args.ks.split(",")]
    
    # finding target column 
    possible_targets = ["target", "Sample_Condition"]

    for col in possible_targets:
        if col in y.columns:
            y = y[col]
            break
    else:
        raise ValueError(f"No valid target found. Tried: {possible_targets}")

    # flattening y into 1D array
    y = y.values.ravel()
      
    
    seed_results = run_classification(X, y, ks, args.n_seeds,args.output, args.label,args.model, args.sampling_strategy, args.grid_search)
    #plot_bar(result_path,args.ks)
    #plot_all_modes(result_path)
    
    #diagnostic plot 
    first_seed = list(seed_results.keys())[0]
    first_k = list(seed_results[first_seed].keys())[0]
    plot_data = seed_results[first_seed][first_k]["plot_data"]
    plot_pairwise_diagnostics(plot_data, args.diagnostic_plot)
   

    # performance plot per feature 
    plot_model_performance_by_features(args.output, args.plot_per_feature)
    
    


    
if __name__ == "__main__":
    main()

 
