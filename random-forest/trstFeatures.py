import csv
import sys

def compare_meta_tmap(version,dataID):
    tmap_path = f"../results/{version}/{dataID}/gffmul.meta.gtf.tmap"
    stat_path = f"./{version}/{dataID}/-1.stats.csv"
    path_feature_path = f"../results/{version}/{dataID}/meta.trstFeature.csv"

    trst_hit_miss = {}

    with open(tmap_path, 'r') as tmap:
        next(tmap)  # skip the first line
        for line in tmap:
            fields = line.strip().split('\t')
            # Extract relevant fields
            _, _, class_code, _, qry_id, *_ = fields
            trst_hit_miss[qry_id] = class_code

    with open(stat_path, 'w', newline='') as stat:
        csv_writer = csv.writer(stat)
        csv_writer.writerow(["tid", "meta_tid", "chr", "cov", "cov2", "abundance", "confidence", "count1", "count2", "num_exons", "gr_vertices", "gr_edges", "gr_reads", "gr_subgraph", "v", "e", "junc_ratio", "max_mid_exon_len", "start_loss1","start_loss2", "start_loss3","end_loss1", "end_loss2","end_loss3", "start_merged_loss", "end_merged_loss", "introns", "intron_ratio","start_introns", "start_intron_ratio","end_introns", "end_intron_ratio", "uni_junc", "seq_min_wt", "seq_min_cnt", "seq_min_abd", "seq_min_ratio","seq_max_wt", "seq_max_cnt", "seq_max_abd", "seq_max_ratio","start_cnt","start_weight","start_abd","end_cnt","end_weight","end_abd","unbridge_start_coming_count", "unbridge_start_coming_ratio", "unbridge_end_leaving_count", "unbridge_end_leaving_ratio", "class_code", "Matched"])

        with open(path_feature_path, 'r') as path_feature:
            for line in path_feature:
                fields = line.strip().split('\t')
                tid, *remaining_fields = fields

                if tid not in trst_hit_miss:
                    print(f"Warning: {tid} is missing")
                    continue

                matched = 1 if trst_hit_miss[tid] == "=" else 0
                stat.write(f"{tid},{','.join(remaining_fields)},{trst_hit_miss[tid]},{matched}\n")

def compare_tmaps(version, dataID, sid):
    tmap_path = f"../results/{version}/{dataID}/gtf/{sid}.mul.{sid}.gtf.tmap"
    stat_path = f"./{version}/{dataID}/{sid}.stats.csv"
    path_feature_path = f"../results/{version}/{dataID}/gtf/{sid}.trstFeature.csv"

    trst_hit_miss = {}

    with open(tmap_path, 'r') as tmap:
        next(tmap)  # skip the first line
        for line in tmap:
            fields = line.strip().split('\t')
            # Extract relevant fields
            _, _, class_code, _, qry_id, *_ = fields
            trst_hit_miss[qry_id] = class_code

    with open(stat_path, 'w', newline='') as stat:
        csv_writer = csv.writer(stat)
        csv_writer.writerow(["tid", "meta_tid", "chr", "cov", "cov2", "abundance", "confidence", "count1", "count2", "num_exons", "gr_vertices", "gr_edges", "gr_reads", "gr_subgraph", "v", "e", "junc_ratio", "max_mid_exon_len", "start_loss1","start_loss2", "start_loss3","end_loss1", "end_loss2","end_loss3", "start_merged_loss", "end_merged_loss", "introns", "intron_ratio","start_introns", "start_intron_ratio","end_introns", "end_intron_ratio", "uni_junc", "seq_min_wt", "seq_min_cnt", "seq_min_abd", "seq_min_ratio", "seq_max_wt", "seq_max_cnt", "seq_max_abd", "seq_max_ratio","start_cnt","start_weight","start_abd","end_cnt","end_weight","end_abd","unbridge_start_coming_count", "unbridge_start_coming_ratio", "unbridge_end_leaving_count", "unbridge_end_leaving_ratio", "class_code", "Matched"])

        with open(path_feature_path, 'r') as path_feature:
            for line in path_feature:
                fields = line.strip().split('\t')  
                tid, *remaining_fields = fields

                if tid not in trst_hit_miss:
                    print(f"Warning: {tid} is missing")
                    continue

                matched = 1 if trst_hit_miss[tid] == "=" else 0
                csv_writer.writerow([tid] + remaining_fields + [trst_hit_miss[tid], matched])




def main():
    if len(sys.argv) != 4:
        print("Usage: python test.py version dataset sample_size")
        sys.exit(1)

    version = sys.argv[1]
    dataID = sys.argv[2]
    sample_size = int(sys.argv[3])  

    print(f"Version: {version}")
    print(f"Dataset: {dataID}")
    print(f"Sample Size: {sample_size}")

    compare_meta_tmap(version,dataID)
    for sid in range(sample_size):
        compare_tmaps(version,dataID, sid)

if __name__ == "__main__":
    main()

