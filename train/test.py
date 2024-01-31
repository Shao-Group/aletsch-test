import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import make_scorer, precision_score, recall_score, precision_recall_curve, accuracy_score, roc_auc_score, auc
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score, cross_validate, StratifiedKFold
import ast
import matplotlib.pyplot as plt
import seaborn as sns
from joblib import load
import argparse

version="aletsch-39"

# Constants / Configurations
DATA_CONFIG_en = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID2": 'encode10_ensembl_other_chrs-star',
    #"dataID2": 'encode10_refseq_other_chrs-star',
    "test_sample_size": 10,
    "default_threshold": 0.5,
    "single_cell": 0,
    "ref_size": 223356,#ensembl
    #"ref_size":170097 #refseq
}

DATA_CONFIG_sc = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID2": 'smartseq3_ensembl_human_chr1_92-star',
    #"dataID2": 'smartseq3_ensembl_human_other_chrs_92-star',
    #"dataID2": 'smartseq3_refseq_human_chr1_92-star',
    #"dataID2": 'smartseq3_refseq_human_other_chrs_92-star',
    "test_sample_size": 92,
    "default_threshold": 0.5,
    "single_cell": 1,
    "ref_size": 223356,#ensembl
    #"ref_size":170097 #refseq
}

DATA_CONFIG_mouse = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID2": 'smartseq3_ensembl_mouse-star',
    #"dataID2": 'smartseq3_refseq_mouse-star',
    "test_sample_size": 369,
    "default_threshold": 0.5,
    "single_cell": 1,
    "ref_size": 120608, #ensembl
    #"ref_size": 128199 #refseq
}

DATA_CONFIG_xpress = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID2": 'smartseq3_ensembl_xpress_run2_other_chrs-star',
    #"dataID2": 'smartseq3_refseq_xpress_run2_other_chrs-star',
    #"dataID2": 'smartseq3_ensembl_xpress_run2_top500_other_chrs-star',
    "test_sample_size": 1066,
    "default_threshold": 0.5,
    "single_cell": 1,
    "ref_size": 223356,#ensembl
    #"ref_size":170097 #refseq
}

DATA_CONFIG_xx230 = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID2": 'PRJNA575230_ensembl_other_chrs-star',
    #"dataID2": 'PRJNA575230_refseq_other_chrs-star',
    "test_sample_size": 73,
    "default_threshold": 0.5,
    "single_cell": 0,
    "ref_size": 223356,#ensembl
    #"ref_size":170097 #refseq
}

DATA_CONFIG_pcsim = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID2": 'psiclass_sim_chr2-hisat2',
    "test_sample_size": 25,
    "default_threshold": 0.5,
    "single_cell": 0,
    "ref_size": 16917
}

DATA_CONFIG_xx891 = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID2": 'PRJNA489891_ensembl_other_chrs-star',
    #"dataID2": 'PRJNA489891_ensembl_chr1-star',
    #"dataID2": 'PRJNA489891_refseq_chr1-star',
    #"dataID2": 'PRJNA489891_refseq_other_chrs-star',
    "train_sample_size": 12,
    "test_sample_size": 12,
    "default_threshold": 0.5,
    "single_cell": 0,
    "ref_size": 223356,#ensembl
    #"ref_size":170097 #refseq
}

DATA_CONFIG_mouse790 = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID2": 'PRJEB18790_ensembl-star',
    #"dataID2": 'PRJEB18790_refseq-star',
    "test_sample_size": 44,
    "default_threshold": 0.5,
    "single_cell": 0,
    "ref_size": 120608, #ensembl
    #"ref_size": 128199 #refseq
}

DATA_CONFIG_poly = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID2": 'polyester_test7_ensembl-star',
    "test_sample_size": 30,
    "default_threshold": 0.5,
    "single_cell": 0,
    "ref_size": 101453
}

configurations = {
    "encode10": DATA_CONFIG_en,
    "sc-human": DATA_CONFIG_sc,
    "xx230": DATA_CONFIG_xx230,
    "sc-xpress": DATA_CONFIG_xpress,
    "sc-mouse": DATA_CONFIG_mouse,
    "mouse790": DATA_CONFIG_mouse790,
    "poly": DATA_CONFIG_poly,
    "xx891": DATA_CONFIG_xx891
}

DATA_CONFIG = None


def load_data(dataID, sample_size, balance):
    features_to_normalize = ['cov', 'abundance', 'count1', 'count2',
            'start_loss1', 'end_loss1', 'start_loss2', 'end_loss2', 
            'start_loss3', 'end_loss3','start_merged_loss', 'end_merged_loss',
            'seq_min_cnt','seq_min_abd', 'seq_max_cnt','seq_max_abd'
            ]    
    frames = []

    # Load all data frames
    for i in range(-1, sample_size):
        df = pd.read_csv(f"{DATA_CONFIG['base_path']}/{version}/{dataID}/{i}.stats.csv", dtype={2: str})
        df['meta_only'] = (i == -1)
        df['sample_id'] = i
        #df['single_cell'] = DATA_CONFIG['single_cell']
        frames.append(df)

    combined_df = pd.concat(frames, ignore_index=True)

    max_cnt = combined_df['count2'].max()
    combined_df[features_to_normalize] = combined_df[features_to_normalize]/max_cnt
    print("max_cnt: ", max_cnt)

    # Downsample the larger group
    if balance:
        pos_cases = combined_df[combined_df['Matched'] == 1]
        neg_cases = combined_df[combined_df['Matched'] == 0]
        if len(pos_cases) < len(neg_cases):
            neg_cases = neg_cases.sample(len(pos_cases), random_state=42)
            combined_df = pd.concat([pos_cases, neg_cases], ignore_index=True)

    combined_df['sample_size'] = sample_size
    combined_df['meta_tid'] = combined_df['meta_tid'].apply(lambda x: f"{dataID}.{x}")

    print("Load and Normalize: ", dataID, combined_df.shape)
    return combined_df

def load_all_data(data_samples, balance):
    all_frames = []
    for dataID, sample_size in data_samples.items():
        data = load_data(dataID, sample_size, balance)
        all_frames.append(data)

    return pd.concat(all_frames, ignore_index=True)

def preprocess_data(data):
    X = data[['cov', 'cov2', 'abundance', 'confidence', 'count1', 'count2',
              'num_exons', 'gr_vertices', 'gr_edges', 'v', 'e', 
              'junc_ratio', 'max_mid_exon_len',
              'start_loss1', 'end_loss1', 'start_loss2', 'end_loss2', 
              'start_loss3', 'end_loss3','start_merged_loss', 'end_merged_loss',
              'introns', 'intron_ratio', 'start_introns', 'end_introns',
              'start_intron_ratio', 'end_intron_ratio','uni_junc', 
              'seq_min_wt','seq_min_cnt','seq_min_abd','seq_min_ratio', 
              'seq_max_wt','seq_max_cnt','seq_max_abd','seq_max_ratio', 
              'meta_only', 'sample_size',
              'start_cnt', 'start_weight', 'start_abd',
              'end_cnt', 'end_weight', 'end_abd',
              'gr_reads', 'gr_subgraph',
              "unbridge_start_coming_count", "unbridge_start_coming_ratio",
              "unbridge_end_leaving_count", "unbridge_end_leaving_ratio",
              #"single_cell"
              ]]
    y = data['Matched']
    return X, y

def compute_meta_data(X_test, y_prob, test_data):
    meta_data = test_data[['meta_tid', 'cov', 'Matched', 'cov2']].copy()
    meta_data['y_prob'] = y_prob
    meta_data['cov'] = meta_data['cov']*(DATA_CONFIG['test_sample_size']+1)
    meta_data['cov2'] = np.exp(meta_data['cov2'])-1

    # Group by meta_tid and compute the mean for y_prob and cov
    meta_data = meta_data.groupby('meta_tid').agg({'y_prob': 'mean', 'cov': 'mean', 'cov2': 'sum', 'Matched': 'first'}).reset_index()

    return meta_data

def read_transmeta_data():
    true_positives = pd.read_csv(f"{DATA_CONFIG['base_path']}/../results/transmeta/{DATA_CONFIG['dataID2']}/transmeta.mat", header=None).iloc[:, 0].tolist()
    precision = pd.read_csv(f"{DATA_CONFIG['base_path']}/../results/transmeta/{DATA_CONFIG['dataID2']}/transmeta.pre", header=None).iloc[:, 0].tolist()
    return true_positives, precision

def read_st2merge_data():
    true_positives = pd.read_csv(f"{DATA_CONFIG['base_path']}/../results/st2merge/{DATA_CONFIG['dataID2']}/st2merge.mat", header=None).iloc[:, 0].tolist()
    precision = pd.read_csv(f"{DATA_CONFIG['base_path']}/../results/st2merge/{DATA_CONFIG['dataID2']}/st2merge.pre", header=None).iloc[:, 0].tolist()
    return true_positives, precision

def read_sc2taco_data():
    true_positives = pd.read_csv(f"{DATA_CONFIG['base_path']}/../results/sc2taco/{DATA_CONFIG['dataID2']}/sc2taco.mat", header=None).iloc[:, 0].tolist()
    precision = pd.read_csv(f"{DATA_CONFIG['base_path']}/../results/sc2taco/{DATA_CONFIG['dataID2']}/sc2taco.pre", header=None).iloc[:, 0].tolist()
    return true_positives, precision

def read_psiclass_data():
    mat_file = f"{DATA_CONFIG['base_path']}/../results/psiclass/{DATA_CONFIG['dataID2']}/psiclass.mat"
    pre_file = f"{DATA_CONFIG['base_path']}/../results/psiclass/{DATA_CONFIG['dataID2']}/psiclass.pre"

    true_positives, precision = [], []

    if os.path.exists(mat_file):
        true_positives = pd.read_csv(mat_file, header=None).iloc[:, 0].tolist()
        precision = pd.read_csv(pre_file, header=None).iloc[:, 0].tolist()

    return true_positives, precision


def calculate_metrics_at_thresholds(meta_data_sorted, tool, true_positives_all, precision_all, default_index):
    print(f"\n----------{tool}----------")
    results1 = []

    true_positives = [true_positives_all[0], true_positives_all[default_index]]
    precision = [precision_all[0], precision_all[default_index]]

    for i, tp_target in enumerate(true_positives):
        cumsum_matched = meta_data_sorted['Matched'].cumsum()
        threshold_index = (cumsum_matched - tp_target).abs().idxmin()
        threshold = meta_data_sorted.loc[threshold_index, 'y_prob']
        relevant_rows = meta_data_sorted[meta_data_sorted['y_prob'] >= threshold]
        tp = relevant_rows['Matched'].sum()
        fp = len(relevant_rows) - tp
        precision_y_prob = tp / (tp + fp) if (tp + fp) != 0 else 0

        recall = tp_target/DATA_CONFIG['ref_size']
        f1score=2 * (precision[i]/100 * recall) / (precision[i]/100 + recall)

        results1.append({
            f'true_positives_{tool}': tp_target,
            f'precision_{tool}(%)': precision[i],
            'true_positives_Aletsch': tp,
            'precision_Aletsch(%)': round(precision_y_prob*100, 1),
            '+(abs,%)': round(precision_y_prob*100, 1)-precision[i],
            f'Fscore_{tool}(%)': round(f1score*100, 1)
        })

    results1_df = pd.DataFrame(results1)
    print(results1_df)

    results2 = []
    for j, precision_target in enumerate(precision):
        precision_target = precision_target/100
        for index in reversed(range(len(meta_data_sorted))):
            tp = meta_data_sorted.iloc[:index + 1]['Matched'].sum()
            fp = (index + 1) - tp
            precision_current = tp / (tp + fp) if (tp + fp) != 0 else 0
            if precision_current >= precision_target:
                break

        recall = true_positives[j] / DATA_CONFIG['ref_size']  
        f1score = 2 * (precision_target * recall) / (precision_target + recall)

        results2.append({
            f'precision_{tool}': precision_target,
            f'true_positives_{tool}':true_positives[j],
            'precision_Aletsch': round(precision_current * 100, 1),
            'true_positives_Aletsch': tp,
            '+(rlt,%)': round((tp-true_positives[j])/true_positives[j]*100, 1),
            f'Fscore_{tool}(%)': round(f1score*100, 1)
        })

    results2_df = pd.DataFrame(results2)
    print(results2_df)
    return results1_df, results2_df

def calculate_for_new_cov(meta_df, covlist):
    results = []

    for cov in covlist:
        subset_df = meta_df[meta_df['new_cov'] >= cov]
        tp = subset_df['Matched'].sum()
        total = len(subset_df)
        precision = tp / total if total > 0 else 0

        recall = tp/DATA_CONFIG['ref_size']
        f1score=2 * (precision * recall) / (precision + recall)
        #print(precision, tp, recall)

        results.append({
            "Cov": cov,
            "Prob": subset_df['y_prob'].min(), 
            "True Positives": tp, 
            "Precision": round(precision*100, 1),
            "F1 score": f1score
            })

    results_df = pd.DataFrame(results)
    print(results_df)

    results_df[['True Positives']].to_csv(f"{DATA_CONFIG['base_path']}/../results/{version}/{DATA_CONFIG['dataID2']}/{version}.mat", index=False, header=False)
    results_df[['Precision']].to_csv(f"{DATA_CONFIG['base_path']}/../results/{version}/{DATA_CONFIG['dataID2']}/{version}.pre", index=False, header=False)
    return results_df

def calculate_for_default_prob(meta_df):
    results = []
    problist=[0.05 * i for i in range(0, 20)]

    for prob in problist:
        subset_df = meta_df[meta_df['y_prob'] >= prob]
        tp = subset_df['Matched'].sum()
        total = len(subset_df)
        precision = tp / total if total > 0 else 0

        recall = tp/DATA_CONFIG['ref_size']
        f1score=2 * (precision * recall) / (precision + recall)
        #print(precision, tp, recall)

        results.append({
            "Prob": prob,
            "True Positives": tp, 
            "Precision": round(precision*100, 1),
            "F1 score": f1score
            })

    results_df = pd.DataFrame(results)
    print(results_df)

    results_df[['True Positives']].to_csv(f"{DATA_CONFIG['base_path']}/../results/{version}/{DATA_CONFIG['dataID2']}/{version}.default.mat", index=False, header=False)
    results_df[['Precision']].to_csv(f"{DATA_CONFIG['base_path']}/../results/{version}/{DATA_CONFIG['dataID2']}/{version}.default.pre", index=False, header=False)
    return results_df

def calculate_prob_thresholds(meta_df_sorted):
    tp = meta_df_sorted['Matched'].sum()
    last_recorded_tp = tp 
    results = []

    # Iterate through descreasing sorted meta_df
    for index in range(len(meta_df_sorted) - 1, -1, -1):
        row = meta_df_sorted.iloc[index]
        if tp == last_recorded_tp:
            total_predictions = index + 1
            precision = tp / total_predictions if total_predictions > 0 else 0

            results.append({
                'True Positives': tp,
                'Precision': round(precision*100, 1)
            })
            last_recorded_tp -= 10

        if row['Matched'] == 1:
            tp -= 1
    
    results_df = pd.DataFrame(results)

    results_df[['True Positives']].to_csv(f"{DATA_CONFIG['base_path']}/../results/{version}/{DATA_CONFIG['dataID2']}/{version}.mat", index=False, header=False)
    results_df[['Precision']].to_csv(f"{DATA_CONFIG['base_path']}/../results/{version}/{DATA_CONFIG['dataID2']}/{version}.pre", index=False, header=False)
    return results_df


def calculate_individual_metrics(group, prediction_col):
    TP = ((group[prediction_col] == 1) & (group['Matched'] == 1)).sum()
    FP = ((group[prediction_col] == 1) & (group['Matched'] == 0)).sum()
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    return pd.Series({'TP': TP, 'Precision': precision})

def calculate_pr_auc(precision_percent, true_positives, total_positives, max_tp, min_tp):
    """
    Calculate the area under the Precision-Recall curve.

    :param precision_percent: List of precision values in percentages.
    :param true_positives: List of true positive counts.
    :param total_positives: Total number of actual positives in the dataset.
    :param max_tp: Maximum true positive count for the range of interest.
    :param min_tp: Minimum true positive count for the range of interest.
    :return: Area under the Precision-Recall curve for the specified range.
    """

    # Convert precision from percentage to proportion
    precision = np.array(precision_percent) / 100

    # Calculate recall
    recall = np.array(true_positives) / total_positives

    # Calculate max and min recall for the specified TP range
    max_recall = max_tp / total_positives
    min_recall = min_tp / total_positives

    # Filter based on recall range
    filter_mask = (recall >= min_recall) & (recall <= max_recall)
    filtered_precision = precision[filter_mask]
    filtered_recall = recall[filter_mask]

    # Sort the filtered data by recall (ascending)
    sorted_indices = np.argsort(filtered_recall)
    sorted_recall = filtered_recall[sorted_indices]
    sorted_precision = filtered_precision[sorted_indices]

    print(sorted_precision, sorted_recall)

    # Print filtered precision and recall
    #print("Filtered Precision:", filtered_precision)
    #print("Filtered Recall:", filtered_recall)

    # Calculate the area under the PR curve
    pr_auc = auc(sorted_recall, sorted_precision)
    print("Area under the PR curve:", pr_auc)

    return pr_auc

def main():
    parser = argparse.ArgumentParser(description="Test Script")
    parser.add_argument('-i', type=str, required=True, choices=configurations.keys(), help="Input test dataset")
    args = parser.parse_args()

    if args.i not in configurations:
        print(f"Error: '{args.i}' is not a valid choice.")
        sys.exit(1)
    
    global DATA_CONFIG
    DATA_CONFIG = configurations[args.i]

    #model = load(f"{DATA_CONFIG['base_path']}/models/{version}.(50).{DATA_CONFIG['dataID1']}.joblib")
    model_version=f"{version}.ensembl.allchrs.(en+xx230+sc-human192+sc-xpress500&1066).joblib"
    print(model_version)
    model = load(f"{DATA_CONFIG['base_path']}/models/{model_version}")

    # Test
    test_data = load_data(DATA_CONFIG['dataID2'], DATA_CONFIG['test_sample_size'], False)
    X_test, y_test = preprocess_data(test_data)
    print(y_test.value_counts())

    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]

    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_prob)
    print(f'Precision: {precision:.3f}')
    print(f'Recall: {recall:.3f}')
    print(f'ROC curve (area = {roc_auc:.3f})')

    individual_data = test_data[['sample_id', 'Matched']].copy()
    individual_data['y_prob'] = y_prob
    individual_data['individual_y_pred'] = (individual_data['y_prob'] >= 0.4).astype(int)
    individual_metrics = individual_data.groupby('sample_id').apply(calculate_individual_metrics, 'individual_y_pred')
    print(individual_metrics)

    meta_df = compute_meta_data(X_test, y_prob, test_data)
    roc_auc_meta = roc_auc_score(meta_df['Matched'], meta_df['y_prob'])
    print(f'Meta ROC curve (area = {roc_auc_meta:.3f})')

    precision_meta_array, recall_meta_array, _ = precision_recall_curve(meta_df['Matched'], meta_df['y_prob'])


    calculate_for_default_prob(meta_df)

    meta_df_sorted = meta_df.sort_values(by='y_prob', ascending=False).reset_index(drop=True)

    #sorted_cov = sorted(meta_df['cov2'], ascending=False)
    #meta_df_sorted['new_cov'] = sorted_cov

    #nor_cov = meta_df_sorted[meta_df_sorted['y_prob'] > 0.4].head(1)['new_cov'].iloc[0]
    #print("normalized cov:", nor_cov)

    #max_new_cov = meta_df_sorted['new_cov'].max()
    #indices_subset1 = meta_df_sorted[meta_df_sorted['new_cov'] >= nor_cov].index
    #indices_subset2 = meta_df_sorted[meta_df_sorted['new_cov'] < nor_cov].index
    #df_subset1 = meta_df_sorted.loc[indices_subset1, 'new_cov']
    #meta_df_sorted.loc[indices_subset1, 'new_cov'] = (df_subset1 - nor_cov) / (max_new_cov - nor_cov) * (max_new_cov - 1) + 1
    #df_subset2 = meta_df_sorted.loc[indices_subset2, 'new_cov']
    #meta_df_sorted.loc[indices_subset2, 'new_cov'] = df_subset2 / nor_cov

    pd.set_option('display.max_rows', 30)
    pd.set_option('display.max_columns', 50)
    pd.set_option('display.width', 1000)

    #covlist = [0, 0.5, 1, 2, 5, 10, 20, 30, 50, 80, 100, 150, 200, 250, 300, 400, 500]
    #calculate_for_new_cov(meta_df_sorted, covlist)

    #interest region
    max_tp = 0
    min_tp = DATA_CONFIG['ref_size']

    tm_true_positives, tm_precision = read_transmeta_data()
    calculate_metrics_at_thresholds(meta_df_sorted, 'transmeta', tm_true_positives, tm_precision, 3)
    tm_max_tp=max(tm_true_positives)
    tm_min_tp=min(tm_true_positives)
    calculate_pr_auc(tm_precision, tm_true_positives, DATA_CONFIG['ref_size'], tm_max_tp, tm_min_tp)
    max_tp = max(max_tp, tm_max_tp)
    min_tp = min(min_tp, tm_min_tp)
    print("Aletsch rpAUC:")
    calculate_pr_auc(precision_meta_array*100, recall_meta_array*np.sum(meta_df['Matched']), DATA_CONFIG['ref_size'], tm_max_tp, tm_min_tp)

    psi_true_positives, psi_precision = read_psiclass_data()
    if(psi_true_positives):
        calculate_metrics_at_thresholds(meta_df_sorted, 'psiclass', psi_true_positives, psi_precision, 2)
        psi_max_tp=max(psi_true_positives)
        psi_min_tp=min(psi_true_positives)
        calculate_pr_auc(psi_precision, psi_true_positives, DATA_CONFIG['ref_size'], psi_max_tp, psi_min_tp)
        max_tp = max(max_tp, psi_max_tp)
        min_tp = min(min_tp, psi_min_tp)
        print("Aletsch rpAUC:")
        calculate_pr_auc(precision_meta_array*100, recall_meta_array*np.sum(meta_df['Matched']), DATA_CONFIG['ref_size'], psi_max_tp, psi_min_tp)

    st2_true_positives, st2_precision = read_st2merge_data()
    calculate_metrics_at_thresholds(meta_df_sorted, 'st2merge', st2_true_positives, st2_precision, 0)
    st2_max_tp=max(st2_true_positives)
    st2_min_tp=min(st2_true_positives)
    calculate_pr_auc(st2_precision, st2_true_positives, DATA_CONFIG['ref_size'], st2_max_tp, st2_min_tp)
    max_tp = max(max_tp, st2_max_tp)
    min_tp = min(min_tp, st2_min_tp)
    print("Aletsch rpAUC:")
    calculate_pr_auc(precision_meta_array*100, recall_meta_array*np.sum(meta_df['Matched']), DATA_CONFIG['ref_size'], st2_max_tp, st2_min_tp)

    sc2_true_positives, sc2_precision = read_sc2taco_data()
    calculate_metrics_at_thresholds(meta_df_sorted, 'sc2taco', sc2_true_positives, sc2_precision, 1)
    sc2_max_tp=max(sc2_true_positives)
    sc2_min_tp=min(sc2_true_positives)
    calculate_pr_auc(sc2_precision, sc2_true_positives, DATA_CONFIG['ref_size'], sc2_max_tp, sc2_min_tp)
    max_tp = max(max_tp, sc2_max_tp)
    min_tp = min(min_tp, sc2_min_tp)
    print("Aletsch rpAUC:")
    calculate_pr_auc(precision_meta_array*100, recall_meta_array*np.sum(meta_df['Matched']), DATA_CONFIG['ref_size'], sc2_max_tp, sc2_min_tp)

    # pAUC
    print("Aletsch pAUC:")
    calculate_pr_auc(precision_meta_array*100, recall_meta_array*np.sum(meta_df['Matched']), DATA_CONFIG['ref_size'], max_tp, 0)


    #generate points to plot Aletsch's curve
    meta_df_sorted = meta_df_sorted[meta_df_sorted['y_prob'] >= 0.2]
    calculate_prob_thresholds(meta_df_sorted)

if __name__ == "__main__":
    main()





