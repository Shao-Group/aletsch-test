#!/bin/bash
version=aletsch-39
aligner=star
mkdir -p ./$version

encode10_datasets=("encode10_ensembl_chr1" "encode10_ensembl_other_chrs" "encode10_ensembl" "encode10_refseq_chr1" "encode10_refseq_other_chrs" "encode10_refseq")
PRJNA575230_datasets=("PRJNA575230_ensembl_chr1" "PRJNA575230_ensembl_other_chrs" "PRJNA575230_ensembl" "PRJNA575230_refseq_chr1" "PRJNA575230_refseq_other_chrs" "PRJNA575230_refseq")
smartseq3_human_datasets=("smartseq3_ensembl_human_chr1_92" "smartseq3_ensembl_human_other_chrs_92" "smartseq3_ensembl_human" "smartseq3_refseq_human_chr1_100" "smartseq3_refseq_human_other_chrs_100" "smartseq3_refseq_human")
smartseq3_xpress_datasets=("smartseq3_ensembl_xpress_run2_top500_chr1" "smartseq3_ensembl_xpress_run2_top500_other_chrs" "smartseq3_ensembl_xpress_run2" "smartseq3_refseq_xpress_run2_top500_chr1" "smartseq3_refseq_xpress_run2_top500_other_chrs" "smartseq3_refseq_xpress_run2")

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

group=$1
index=$2

case $group in
    encode10)
        if [[ $index -ge 0 && $index -lt ${#encode10_datasets[@]} ]]; then
            dataset=${encode10_datasets[$index]}
            sample_size=10
        else
            echo "Invalid index for encode10. Please choose between 0 and $((${#encode10_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    xx230)
        if [[ $index -ge 0 && $index -lt ${#PRJNA575230_datasets[@]} ]]; then
            dataset=${PRJNA575230_datasets[$index]}
            sample_size=73
        else
            echo "Invalid index for PRJNA575230. Please choose between 0 and $((${#PRJNA575230_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    sc-human)
        if [[ $index -ge 0 && $index -lt ${#smartseq3_human_datasets[@]} ]]; then
            dataset=${smartseq3_human_datasets[$index]}
            sample_size=192
        else
            echo "Invalid index for smartseq3. Please choose between 0 and $((${#smartseq3_human_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    sc-mouse)
        if [[ $index -ge 0 && $index -lt ${#smartseq3_mouse_datasets[@]} ]]; then
            dataset=${smartseq3_mouse_datasets[$index]}
            sample_size=369
        else
            echo "Invalid index for smartseq3. Please choose between 0 and $((${#smartseq3_mouse_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    sc-xpress)
        if [[ $index -ge 0 && $index -lt ${#smartseq3_xpress_datasets[@]} ]]; then
            dataset=${smartseq3_xpress_datasets[$index]}
            sample_size=500
        else
            echo "Invalid index for smartseq3. Please choose between 0 and $((${#smartseq3_xpress_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    pc-sim)
        if [[ $index -ge 0 && $index -lt ${#psiclass_sim_chr2_datasets[@]} ]]; then
            dataset=${psiclass_sim_chr2_datasets[$index]}
            sample_size=25
        else
            echo "Invalid index for psiclass simulated dataset. Please choose between 0 and $((${#psiclass_sim_chr2_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    xx891)
        if [[ $index -ge 0 && $index -lt ${#PRJNA489891_datasets[@]} ]]; then
            dataset=${PRJNA489891_datasets[$index]}
            sample_size=12
        else
            echo "Invalid index for PRJNA489891. Please choose between 0 and $((${#PRJNA489891_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    mouse790)
        if [[ $index -ge 0 && $index -lt ${#PRJEB18790_datasets[@]} ]]; then
            dataset=${PRJEB18790_datasets[$index]}
            sample_size=44
        else
            echo "Invalid index for PRJEB18790. Please choose between 0 and $((${#PRJEB18790_datasets[@]} - 1))."
            exit 1
        fi
        ;;
    poly)
        if [[ $index -ge 0 && $index -lt ${#polyester_datasets[@]} ]]; then
            dataset=${polyester_datasets[$index]}
            sample_size=30
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

mkdir -p ./$version/$dataset-$aligner
python trstFeatures.py $version $dataset-$aligner $sample_size
