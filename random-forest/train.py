import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import make_scorer, precision_score, recall_score, precision_recall_curve, accuracy_score, roc_auc_score
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score, cross_validate, StratifiedKFold
import ast
import matplotlib.pyplot as plt
import seaborn as sns
from joblib import dump

# Constants / Configurations
DATA_CONFIG_en = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID1": 'encode10_ensembl_chr1-star',
    "dataID2": 'encode10_ensembl_other_chrs-star',
    "train_sample_size": 10,
    "test_sample_size": 10,
    "default_threshold": 0.5,
    "single_cell": 0
}

DATA_CONFIG_sc = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID1": 'smartseq3_ensembl_human_chr1-star',
    "dataID2": 'smartseq3_ensembl_human_other_chrs-star',
    #"dataID1": 'smartseq3_refseq_human_chr1-star',
    #"dataID2": 'smartseq3_refseq_human_other_chrs-star',
    "train_sample_size": 192,
    "test_sample_size": 192,
    "default_threshold": 0.5,
    "single_cell": 1
}

DATA_CONFIG_mouse = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID1": 'smartseq3_ensembl_mouse_chr1-star',
    "dataID2": 'smartseq3_ensembl_mouse_other_chrs-star',
    "train_sample_size": 369,
    "test_sample_size": 369,
    "test_chr": 'other_chrs',
    "default_threshold": 0.5,
    "single_cell": 1
}

DATA_CONFIG_xpress = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID1": 'smartseq3_ensembl_xpress_run2_top500_chr1-star',
    "dataID2": 'smartseq3_ensembl_xpress_run2_top500_other_chrs-star',
    "train_sample_size": 500,
    "test_sample_size": 500,
    "default_threshold": 0.5,
    "single_cell": 1
}

DATA_CONFIG_xx230 = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID1": 'PRJNA575230_ensembl_chr1-star',
    "dataID2": 'PRJNA575230_ensembl_other_chrs-star',
    "train_sample_size": 73,
    "test_sample_size": 73,
    "default_threshold": 0.5,
    "single_cell": 0
}

DATA_CONFIG_poly = {
    "base_path": "/data/qzs23/projects/aletsch-test/train",
    "dataID1": 'polyester_train5_ensembl-star',
    "train_sample_size": 30,
    "default_threshold": 0.5,
    "single_cell": 0
}

data_samples = {
    'smartseq3_ensembl_human-star': 192, 
    'smartseq3_ensembl_human_other_chrs-star': 192,
    #'smartseq3_ensembl_human_chr1_92-star': 92, 
    #'smartseq3_ensembl_human_other_chrs_92-star': 92,
    #'smartseq3_ensembl_human_chr1_100-star': 100, 
    #'smartseq3_ensembl_human_other_chrs_100-star': 100,
    'smartseq3_ensembl_xpress_run2_top500_chr1-star': 500,
    'smartseq3_ensembl_xpress_run2_top500_other_chrs-star': 500,
    'smartseq3_ensembl_xpress_run2_chr1-star': 1066,
    'smartseq3_ensembl_xpress_run2_other_chrs-star': 1066,
    'PRJNA575230_ensembl_chr1-star': 73,
    'PRJNA575230_ensembl_other_chrs-star': 73,
    'encode10_ensembl_chr1-star': 10,
    'encode10_ensembl_other_chrs-star': 10
    #'smartseq3_refseq_human_chr1_100-star': 100, 
    #'smartseq3_refseq_human_other_chrs-star': 192,
    #'smartseq3_refseq_xpress_run2_chr1-star': 1066,
    #'smartseq3_refseq_xpress_run2_other_chrs-star': 1066,
    #'smartseq3_refseq_xpress_run2_top500_chr1-star': 500,
    #'smartseq3_refseq_xpress_run2_top500_other_chrs-star': 500,
    #'PRJNA575230_refseq_chr1-star': 73,
    #'PRJNA575230_refseq_other_chrs-star': 73,
    #'encode10_refseq_chr1-star': 10,
    #'encode10_refseq_other_chrs-star': 10
}

DATA_SC = {
    'smartseq3_ensembl_human_chr1-star': 1, 
    #'smartseq3_ensembl_human_other_chrs-star': 1,
    #'smartseq3_ensembl_human_chr1_92-star': 92, 
    #'smartseq3_ensembl_human_other_chrs_92-star': 92,
    #'smartseq3_ensembl_human_chr1_100-star': 1, 
    #'smartseq3_ensembl_human_other_chrs_100-star': 100,
    #'smartseq3_ensembl_xpress_run2_top500_chr1-star': 1,
    #'smartseq3_ensembl_xpress_run2_top500_other_chrs-star': 1,
    'PRJNA575230_ensembl_chr1-star': 0,
    'PRJNA575230_ensembl_other_chrs-star': 0,
    'encode10_ensembl_chr1-star': 0,
    'encode10_ensembl_other_chrs-star': 0
}
DATA_CONFIG = DATA_CONFIG_poly
version="aletsch-39"

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
        #df['single_cell'] = DATA_SC[dataID]
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
    meta_data = test_data[['meta_tid', 'cov', 'Matched']].copy()
    meta_data['y_prob'] = y_prob

    # Group by meta_tid and compute the mean for y_prob and cov
    meta_data = meta_data.groupby('meta_tid').agg({'y_prob': 'mean', 'cov': 'mean', 'Matched': 'first'}).reset_index()

    return meta_data

def main():
    # Train
    #train_data = load_data(DATA_CONFIG['dataID1'], DATA_CONFIG['train_sample_size'], True)
    train_data = load_all_data(data_samples, True)
    X_train, y_train = preprocess_data(train_data)
    print(X_train.shape)
    print(y_train.value_counts())

    model = RandomForestClassifier(n_estimators=100, max_depth=20, random_state=42)
    model.fit(X_train, y_train)
    y_train_pred = model.predict(X_train)
    y_train_prob = model.predict_proba(X_train)[:, 1]

    precision = precision_score(y_train, y_train_pred)
    recall = recall_score(y_train, y_train_pred)
    roc_auc = roc_auc_score(y_train, y_train_prob)
    print(f'Train Precision: {precision:.3f}')
    print(f'Train Recall: {recall:.3f}')
    print(f'Train ROC curve (area = {roc_auc:.3f})')

    meta_train_df = compute_meta_data(X_train, y_train_prob, train_data)

    roc_auc = roc_auc_score(meta_train_df['Matched'], meta_train_df['y_prob'])
    print(f'Train Meta ROC curve (area = {roc_auc:.3f})')

    #dump(model, f"{DATA_CONFIG['base_path']}/models/{version}.(50).{DATA_CONFIG['dataID1']}.joblib")
    dump(model, f"{DATA_CONFIG['base_path']}/models/{version}.refseq.chr1.(en+xx230+sc-human100+sc-xpress500&1066).joblib")


if __name__ == "__main__":
    main()





