#!/bin/bash
covlist='0 0.5 1 2 5 10 20 30 50 80 100 150 200 250 300 400 500'

# Common Configuration
dir=`pwd`
stringtie2=$dir/../programs/stringtie
gffcompare=$dir/../programs/gffcompare
gtfcuff=$dir/../programs/gtfcuff
results=$dir/st2merge
aligner="star"

# Check for required tools
missing_tools=()
[ ! -f "$stringtie2" ] && missing_tools+=("stringtie2")
[ ! -f "$gffcompare" ] && missing_tools+=("gffcompare")
[ ! -f "$gtfcuff" ] && missing_tools+=("gtfcuff")

if [ ${#missing_tools[@]} -ne 0 ]; then
    echo "Missing required tools: ${missing_tools[*]}. Exiting."
    exit 1
fi

encode10_datasets=("encode10_ensembl_chr1" "encode10_ensembl_other_chrs" "encode10_ensembl" "encode10_refseq_chr1" "encode10_refseq_other_chrs" "encode10_refseq")
PRJNA575230_datasets=("PRJNA575230_ensembl_chr1" "PRJNA575230_ensembl_other_chrs" "PRJNA575230_ensembl" "PRJNA575230_refseq_chr1" "PRJNA575230_refseq_other_chrs" "PRJNA575230_refseq")
smartseq3_human_datasets=("smartseq3_ensembl_human_chr1" "smartseq3_ensembl_human_other_chrs" "smartseq3_ensembl_human" "smartseq3_ensembl_human_chr1_92" "smartseq3_ensembl_human_other_chrs_92" "smartseq3_ensembl_human_chr1_100" "smartseq3_ensembl_human_other_chrs_100" "smartseq3_refseq_human_chr1" "smartseq3_refseq_human_other_chrs" "smartseq3_refseq_human" "smartseq3_refseq_human_chr1_92" "smartseq3_refseq_human_other_chrs_92" "smartseq3_refseq_human_chr1_100" "smartseq3_refseq_human_other_chrs_100")
smartseq3_xpress_datasets=("smartseq3_ensembl_xpress_run2_chr1" "smartseq3_ensembl_xpress_run2_other_chrs" "smartseq3_ensembl_xpress_run2" "smartseq3_ensembl_xpress_run2_top500_chr1" "smartseq3_ensembl_xpress_run2_top500_other_chrs" "smartseq3_refseq_xpress_run2_chr1" "smartseq3_refseq_xpress_run2_other_chrs" "smartseq3_refseq_xpress_run2" "smartseq3_refseq_xpress_run2_top500_chr1" "smartseq3_refseq_xpress_run2_top500_other_chrs")

smartseq3_mouse_datasets=("smartseq3_ensembl_mouse" "smartseq3_refseq_mouse")
psiclass_sim_chr2_datasets=("psiclass_sim_chr2")
PRJNA489891_datasets=("PRJNA489891_ensembl_chr1" "PRJNA489891_ensembl_other_chrs" "PRJNA489891_ensembl" "PRJNA489891_refseq_chr1" "PRJNA489891_refseq_other_chrs" "PRJNA489891_refseq")
PRJEB18790_datasets=("PRJEB18790_ensembl" "PRJEB18790_refseq")
polyester_datasets=("polyester_test1_ensembl" "polyester_test3_ensembl" "polyester_test5_ensembl" "polyester_test7_ensembl" "polyester_test9_ensembl" "polyester_train5_ensembl")

# Function to display dataset choices
display_dataset_choices() {
    local datasets=("$@")
    for i in "${!datasets[@]}"; do
        echo "  $i: ${datasets[$i]}"
    done
    echo "  all: submit all choices in this group"
}

# Check for required input parameters
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 {encode10|xx230|sc-human|sc-mouse|sc-xpress} dataset_index"
    echo "For example: $0 xx230 3"
    echo "Dataset Choices:"
    echo "  For encode10(encode10):"
    display_dataset_choices "${encode10_datasets[@]}"
    echo "  For xx230(PRJNA575230):"
    display_dataset_choices "${PRJNA575230_datasets[@]}"
    echo "  For sc-human(smartseq3-human):"
    display_dataset_choices "${smartseq3_human_datasets[@]}"
    echo "  For sc-mouse(smartseq3-mouse):"
    display_dataset_choices "${smartseq3_mouse_datasets[@]}"
    echo "  For sc-xpress(smartseq3-xpress):"
    display_dataset_choices "${smartseq3_xpress_datasets[@]}"
    echo "  For pc-sim(psiclass simulated chr2):"
    display_dataset_choices "${psiclass_sim_chr2_datasets[@]}"
    echo "  For xx891(PRJNA489891):"
    display_dataset_choices "${PRJNA489891_datasets[@]}"
    echo "  For mouse790(PRJEB18790):"
    display_dataset_choices "${PRJEB18790_datasets[@]}"  
    echo "  For poly(polyester simulated):"
    display_dataset_choices "${polyester_datasets[@]}"  
    exit 1
fi
process_dataset() {
    dataset=$1

    case $dataset in
        encode10_ensembl*)
            strand="--rf"
            ref=/data/shared/data/ensembl/release-107/Homo_sapiens.GRCh38.107.gtf
            ;;
        encode10_refseq*)
            strand="--rf"
            ref=/data/shared/data/refseq/Homo_sapiens_GRCh38.p14/genomic.gtf
            ;;
        PRJNA575230_ensembl*)
            strand="--rf"
            ref=/data/shared/data/ensembl/release-107/Homo_sapiens.GRCh38.107.gtf
            ;;
        PRJNA575230_refseq*)
            strand="--rf"
            ref=/data/shared/data/refseq/Homo_sapiens_GRCh38.p14/genomic.gtf
            ;;
        smartseq3_ensembl_human*)
            strand=
            ref=/data/shared/data/ensembl/release-107/Homo_sapiens.GRCh38.107.gtf
            ;;
        smartseq3_refseq_human*)
            strand=
            ref=/data/shared/data/refseq/Homo_sapiens_GRCh38.p14/genomic.gtf
            ;;
        smartseq3_ensembl_mouse)
            strand=
            ref=/data/shared/data/mouse/Mus_musculus.GRCm39.110.gtf            
            ;;
        smartseq3_refseq_mouse)
            strand=
            ref=/data/shared/data/refseq/Mus_musculus_GRCm39/genomic.gtf           
            ;;
        smartseq3_ensembl_xpress*)
            strand="--fr"
            ref=/data/shared/data/ensembl/release-107/Homo_sapiens.GRCh38.107.gtf
            ;;
        smartseq3_refseq_xpress*)
            strand="--fr"
            ref=/data/shared/data/refseq/Homo_sapiens_GRCh38.p14/genomic.gtf
            ;;
        psiclass_sim_chr2*)
            strand=
            ref=/data/qzs23/data/psiclass_sim_chr2/Homo_sapiens.GRCh38.107.chr2.gtf
            aligner="hisat2"
            ;;
        PRJNA489891_ensembl*)
            strand=unstranded
            ref=/data/shared/data/ensembl/release-107/Homo_sapiens.GRCh38.107.gtf
            ;;
        PRJNA489891_refseq*)
            strand=unstranded
            ref=/data/shared/data/refseq/Homo_sapiens_GRCh38.p14/genomic.gtf
            ;;
        PRJEB18790_ensembl*)
            strand=unstranded
            ref=/data/shared/data/mouse/Mus_musculus.GRCm39.110.gtf
            ;;
        PRJEB18790_refseq*)
            strand=unstranded
            ref=/data/shared/data/refseq/Mus_musculus_GRCm39/genomic.gtf
            ;; 
        polyester_*)
            strand=
            ref=/data/qzs23/data/polyester_simulated2/annotations/Homo_sapiens.GRCh38.107.trst-min3.half2.gtf
            ;;       
        *)
            echo "Unknown dataset: $dataset"
            exit 1
            ;;
    esac

    echo "Processing dataset: $dataset"

    cur="$results/$dataset-$aligner"
    list="$cur/st.gtf.list"

    script=$cur/$dataset-$aligner.sh
    cd $cur
    rm -f st2merge.mat
    rm -f st2merge.pre

    #find $cur -maxdepth 1 -type f ! -name "st.gtf.list" -exec rm -f {} +
    rm -rf $cur/tmp_*
    echo "cd $cur" > $script

    temp_file="$cur/st.gtf.tmp"
    rm -f $temp_file
    while IFS= read -r gtf_path; do
        if [ $(wc -l < "$gtf_path") -gt 2 ]; then
            echo "$gtf_path" >> "$temp_file"
        fi
    done < "$list"
    mv "$temp_file" "$list"

    path_count=$(wc -l < "$list")
    if [ "$path_count" -gt 1000 ]; then
        head -n 1000 "$list" > "${list}1"
        tail -n +1001 "$list" > "${list}2"
        echo "Large list(#$path_count) split into ${list}1 and ${list}2"
    fi

    for cov in $covlist;
    do
        stats=$cur/$cov.mul.stats

        if [[ -f "${list}1" && -f "${list}2" ]]; then
            echo "$cur/st2merge-list1-$cov.gtf" > $list-$cov
            echo "$cur/st2merge-list2-$cov.gtf" >> $list-$cov
            echo "{ /usr/bin/time -v $stringtie2 --merge ${list}1 -o $cur/st2merge-list1-$cov.gtf -c $cov -p 10 ; } 2> st2merge-list1-$cov.time" >> $script
            echo "{ /usr/bin/time -v $stringtie2 --merge ${list}2 -o $cur/st2merge-list2-$cov.gtf -c $cov -p 10 ; } 2> st2merge-list2-$cov.time" >> $script
            echo "{ /usr/bin/time -v $stringtie2 --merge $list-$cov -o $cur/st2merge-$cov.gtf -c $cov -p 10 ; } 2> st2merge-$cov.time" >> $script
        else
            echo "{ /usr/bin/time -v $stringtie2 --merge $list -o $cur/st2merge-$cov.gtf -c $cov -p 10 ; } 2> st2merge-$cov.time" >> $script
        fi

        echo "$gffcompare -M -N -r $ref -o $stats st2merge-$cov.gtf" >> $script
        echo "less $stats | awk '{print \$4}' | sed -n '18p' >> st2merge.mat" >> $script
        echo "less $stats | awk '{print \$6}' | sed -n '14p' >> st2merge.pre" >> $script
    done    
    chmod +x $script

    nohup bash $script > "${cur}/nohup.log" 2>&1 &
    echo "Submit job: StringTie2-merge for $dataset"
}

submit_all_jobs() {
    local datasets=("$@")
    for dataset in "${datasets[@]}"; do
        process_dataset "$dataset"
    done
}

group=$1
index=$2

if [ "$index" = "all" ]; then
    echo "Submitting jobs for all datasets..."
    case $group in 
        encode10)
            submit_all_jobs "${encode10_datasets[@]}"
            ;;
        xx230)
            submit_all_jobs "${PRJNA575230_datasets[@]}"
            ;;
        sc-human)
            submit_all_jobs "${smartseq3_human_datasets[@]}"
            ;;
        sc-mouse)
            submit_all_jobs "${smartseq3_mouse_datasets[@]}"
            ;;
        sc-xpress)
            submit_all_jobs "${smartseq3_xpress_datasets[@]}"
            ;;
        pc-sim)
            submit_all_jobs "${psiclass_sim_chr2_datasets[@]}"
            ;;
        xx891)
            submit_all_jobs "${PRJNA489891_datasets[@]}"
            ;;
        mouse790)
            submit_all_jobs "${PRJEB18790_datasets[@]}"
            ;;
        poly)
            submit_all_jobs "${polyester_datasets[@]}"
            ;;
        *)
            exit 1
    esac
    exit 0
fi

case $group in
    encode10)
        if [[ $index -ge 0 && $index -lt ${#encode10_datasets[@]} ]]; then
            dataset=${encode10_datasets[$index]}
            process_dataset "$dataset"
        else
            echo "Invalid index for encode10. Please choose between 0 and $((${#encode10_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    xx230)
        if [[ $index -ge 0 && $index -lt ${#PRJNA575230_datasets[@]} ]]; then
            dataset=${PRJNA575230_datasets[$index]}
            process_dataset "$dataset" 
        else
            echo "Invalid index for PRJNA575230. Please choose between 0 and $((${#PRJNA575230_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    sc-human)
        if [[ $index -ge 0 && $index -lt ${#smartseq3_human_datasets[@]} ]]; then
            dataset=${smartseq3_human_datasets[$index]}
            process_dataset "$dataset" 
        else
            echo "Invalid index for smartseq3. Please choose between 0 and $((${#smartseq3_human_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    sc-mouse)
        if [[ $index -ge 0 && $index -lt ${#smartseq3_mouse_datasets[@]} ]]; then
            dataset=${smartseq3_mouse_datasets[$index]}
            process_dataset "$dataset" 
        else
            echo "Invalid index for smartseq3. Please choose between 0 and $((${#smartseq3_mouse_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    sc-xpress)
        if [[ $index -ge 0 && $index -lt ${#smartseq3_xpress_datasets[@]} ]]; then
            dataset=${smartseq3_xpress_datasets[$index]}
            process_dataset "$dataset" 
        else
            echo "Invalid index for smartseq3. Please choose between 0 and $((${#smartseq3_xpress_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    pc-sim)
        if [[ $index -ge 0 && $index -lt ${#psiclass_sim_chr2_datasets[@]} ]]; then
            dataset=${psiclass_sim_chr2_datasets[$index]}
            process_dataset "$dataset" 
        else
            echo "Invalid index for psiclass simulated dataset. Please choose between 0 and $((${#psiclass_sim_chr2_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    xx891)
        if [[ $index -ge 0 && $index -lt ${#PRJNA489891_datasets[@]} ]]; then
            dataset=${PRJNA489891_datasets[$index]}
            process_dataset "$dataset"         
        else
            echo "Invalid index for PRJNA489891. Please choose between 0 and $((${#PRJNA489891_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    mouse790)
        if [[ $index -ge 0 && $index -lt ${#PRJEB18790_datasets[@]} ]]; then
            dataset=${PRJEB18790_datasets[$index]}
            process_dataset "$dataset"         
        else
            echo "Invalid index for PRJEB18790. Please choose between 0 and $((${#PRJEB18790_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    poly)
        if [[ $index -ge 0 && $index -lt ${#polyester_datasets[@]} ]]; then
            dataset=${polyester_datasets[$index]}
            process_dataset "$dataset"         
        else
            echo "Invalid index for polyester. Please choose between 0 and $((${#polyester_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    *)
        echo "Invalid dataset group. Please choose 'encode10' or 'PRJNA575230' or 'smartseq3'."
        exit 1
        ;;
esac
