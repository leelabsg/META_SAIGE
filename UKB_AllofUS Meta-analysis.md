## Example: Running Meta-Analysis between UKBB and AllofUS

## Description
To perform integrated association analysis across UKBB and AllofUS, there is a difficulty in that both platforms are not compatible for jointly using individual-level genotypes. This page describes example scripts running a meta-analysis between UKBB and AllofUS using META-SAIGE in the [AllofUS Workbench system](https://workbench.researchallofus.org/). 

In this example, we are going to talk about how to run META-SAIGE for two phenotypes: T2D(250.2) and colorectal cancer(153). 

## Workflow Overview

![image](https://github.com/user-attachments/assets/f68437ba-1987-4a7e-86e1-49e51f857245)

## Summary statistics generation from allofUS workbench
### input files description
1. UKB GWAS summary information
   
   Genotypes required for generation of the following information can be acquired from [UKBiobank RAP](https://ukbiobank.dnanexus.com/). 
   - Single-variant Level GWAS
   - Sparse LD Matrix Generation
3. AllofUS summary information
   
   For general help in using workbench and dsub commands, please refer to snippets and sample workspace provided by AllofUS workbench.
   Following is the detailed process on each step of SAIGE and META-SAIGE.


### Helper file for running batch jobs in AllofUS dsub.
```
%%writefile ~/aou_dsub.bash

#!/bin/bash

# This shell function passes reasonable defaults for several dsub parameters, while
# allowing the caller to override any of them. It creates a nice folder structure within
# the workspace bucket for dsub log files.

# --[ Parameters ]--
# any valid dsub parameter flag

#--[ Returns ]--
# the job id of the job created by dsub

#--[ Details ]--
# The first five parameters below should always be those values when running on AoU RWB.

# Feel free to change the values for --user, --regions, --logging, and --image if you like.

# Note that we insert some job data into the logging path.
# https://github.com/DataBiosphere/dsub/blob/main/docs/logging.md#inserting-job-data

function aou_dsub () {

  # Get a shorter username to leave more characters for the job name.
  local DSUB_USER_NAME="$(echo "${OWNER_EMAIL}" | cut -d@ -f1)"

  # For AoU RWB projects network name is "network".
  #local AOU_NETWORK=network
  #local AOU_SUBNETWORK=subnetwork

  dsub \
      --provider google-batch \
      --user-project "${GOOGLE_PROJECT}"\
      --project "${GOOGLE_PROJECT}"\
      --image 'ubuntu:latest' \
      --network "global/networks/network" \
      --subnetwork "regions/us-central1/subnetworks/subnetwork" \
      --service-account "$(gcloud config get-value account)" \
      --use-private-address \
      --user "${DSUB_USER_NAME}" \
      --regions us-central1 \
      --logging "${WORKSPACE_BUCKET}/dsub/logs/{job-name}/{user-id}/$(date +'%Y%m%d/%H%M%S')/{job-id}-{task-id}-{task-attempt}.log" \
      "$@"
}
```
### Parameter description for dsub commands
- `name` : job name
- `provider` : gcp provider (Now google-batch)
- `image` : docker image identifier for gcr.io, and docker hub
- `logging` : path to logging
- `mount` : Designate bucket name for mount in the system
- `boot-disk-size` : Disk size of boot
- `disk-size ` : Disk size per machine
- `min-ram': Minimum size of memory(GB)
- `min-cores': Minimum number of cores
- `task`: delimitered table file path including task parameters, environments. 
- `command` : string command or path to bash script file for running batch jobs.

### GWAS summary generation Step for AllofUS
```
# Example command for SAIGE in AllofUS workbench dsub
# Step 1: Fitting the Null Model
bucket = os.getenv('WORKSPACE_BUCKET')
gwas_dir = 'SAIGE_GENE'
USER_NAME = os.getenv('OWNER_EMAIL').split('@')[0].replace('.','-')
%env USER_NAME={USER_NAME}
JOB_NAME=f'saige-step1-{USER_NAME}'
%env JOB_NAME={JOB_NAME}
%env PHEN={phen_name}

num_traits = 2
traits = ['250.2','153']

# Parameter data frame for running batch jobs
params_df = pd.DataFrame(data={
    '--input-recursive INPUT_DIR': [f"{bucket}/{gwas_dir}/step1_input/"]*num_traits,
    '--input-recursive GRM_DIR': [f"{bucket}/{gwas_dir}/step0_output/"]*num_traits,
    '--output-recursive OUT_DIR': [f"{bucket}/{gwas_dir}/step1_output/"]*num_traits,
    '--env traitType':['binary' for x in traits],
    '--env invNormalize':[True if x.__contains__('f.') else False for x in traits],    
    '--env PHEN': traits
})

PARAMETER_FILENAME = f'{JOB_NAME}_params.tsv'
%env PARAMETER_FILENAME={PARAMETER_FILENAME}

params_df.to_csv(PARAMETER_FILENAME, sep='\t', index=False)

# Pruned genotypes were prepared and used from called genotypes
job_output = !source ~/aou_dsub.bash; aou_dsub \
  --name "${JOB_NAME}_${PHEN}" \
  --provider google-batch \
  --image "wzhou88/saige:1.3.0" \
  --logging "${WORKSPACE_BUCKET}/dsub_logs/saige_step1" \
  --boot-disk-size 50 \
  --disk-size 128 \
  --min-ram 32 \
  --min-cores 8 \
  --tasks "${PARAMETER_FILENAME}" \
  --command 'step1_fitNULLGLMM.R \
                --bedFile=$INPUT_DIR/pruned_arrays_eur.bed \
                --bimFile=$INPUT_DIR/pruned_arrays_eur_new.bim \
                --famFile=$INPUT_DIR/pruned_arrays_eur.fam \
                --phenoFile=$INPUT_DIR/eur_basic_traits.tsv  \
                --invNormalize=$invNormalize \
                --phenoCol=${PHEN} \
                --covarColList=SEX,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
                --qCovarColList=SEX \
                --sampleIDColinphenoFile=IID \
                --traitType=$traitType \
                --isCateVarianceRatio=TRUE \
                --outputPrefix=$OUT_DIR/${PHEN}_step1_output \
                --IsOverwriteVarianceRatioFile=TRUE \
                --LOCO=FALSE \
                --nThreads=8'
job_output

# Step 2: Single variant association testing

USER_NAME = os.getenv('OWNER_EMAIL').split('@')[0].replace('.','-')
%env USER_NAME={USER_NAME}
JOB_NAME=f'saige-step2-{USER_NAME}'
%env JOB_NAME={JOB_NAME}
bucket = os.getenv('WORKSPACE_BUCKET')
gwas_dir = 'SAIGE_GENE'


exome_bgen_dir = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome_v7.1/bgen"

params_df = pd.DataFrame(data={
    '--input-recursive STEP1_OUT': [f"{bucket}/{gwas_dir}/step1_output" for _ in range(22*num_traits)],
    '--input bgenFile': [f"{exome_bgen_dir}/exome.chr{x+1}.bgen" for x in range(22)]*num_traits,
    '--input bgenFileIndex': [f"{exome_bgen_dir}/exome.chr{x+1}.bgen.bgi" for x in range(22)]*num_traits,
    '--input sampleFile': [f"{exome_bgen_dir}/exome.chr{x+1}.sample" for x in range(22)]*num_traits,
    '--input groupFile': [f"{bucket}/{gwas_dir}/step2_input/UKB470k_chr_{x+1}_groupfile.loftee.edit.txt" for x in range(22)]*num_traits,
    '--env PHEN': [x for x in traits for _ in range(22)],
    '--env CHROM': [x+1 for x in range(22)]*num_traits,
    '--output-recursive OUT_DIR': [f"{bucket}/{gwas_dir}/step2_output/" for _ in range(22*num_traits)]
})


PARAMETER_FILENAME = f'{JOB_NAME}_params.tsv'
%env PARAMETER_FILENAME={PARAMETER_FILENAME}

params_df.to_csv(PARAMETER_FILENAME, sep='\t', index=False)
job_id = !source ~/aou_dsub.bash; aou_dsub \
  --name "${JOB_NAME}" \
  --provider google-batch \
  --image "gcr.io/pheweb/saige:1.3.0" \
  --logging "${WORKSPACE_BUCKET}/dsub_logs/saige_step2" \
  --boot-disk-size 100 \
  --disk-size 200 \
  --min-ram 32 \
  --tasks "${PARAMETER_FILENAME}" \
  --command 'step2_SPAtests.R \
                --bgenFile=${bgenFile} \
                --bgenFileIndex=${bgenFileIndex} \
                --minMAF=0 \
                --AlleleOrder=ref-first \
                --is_output_moreDetails=TRUE \
                --is_overwrite_output=TRUE \
                --chrom=chr${CHROM} \
                --GMMATmodelFile=$STEP1_OUT/${PHEN}_step1_output.rda \
                --varianceRatioFile=$STEP1_OUT/${PHEN}_step1_output.varianceRatio.txt \
                --LOCO=FALSE \
                --SAIGEOutputFile=$OUT_DIR/${PHEN}_chr${CHROM}_step2_output_single.txt \
                '
print("\n".join(job_id))
job_id = job_id[1].split(" ")[-1]
%env JOB_ID={job_id}
```
  
### Sparse LD Matrix Generation

```
cmd_line = 'step3_LDmat.R \
    --bgenFile=${bgenFile} \
    --bgenFileIndex=${bgenFileIndex} \
    --sample_include_inLDMat_File=$ids \
    --AlleleOrder=ref-first \
    --chrom=$CHROM \
    --SAIGEOutputFile=$OUT_DIR/${maf}_${anno}_chr${CHR} \
    --groupFile=$groupFile/UKB470k_chr_${CHR}_groupfile.txt \
    --annotation_in_groupTest=${annotation_in_groupTest} \
    --is_overwrite_output=TRUE \
    --maxMAF_in_groupTest=${maf} \
    '
with open("Cmd_step3.sh", "w") as text_file:
    text_file.write(cmd_line)

USER_NAME = os.getenv('OWNER_EMAIL').split('@')[0].replace('.','-')
%env USER_NAME={USER_NAME}
JOB_NAME=f'saige-step3-{USER_NAME}'
%env JOB_NAME={JOB_NAME}
bucket = os.getenv('WORKSPACE_BUCKET')
gwas_dir = 'SAIGE_GENE'

exome_plink_dir = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome/plink_bed"

exome_bgen_dir = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome_v7.1/bgen"
anno = ['missense_lof_synonymous']
chr = [x+1 for x in range(22)]
mafs=([0.01)*22
i=0
print(f"{bucket}/{gwas_dir}/step3_output/{mafs[i]}_{anno[i]}_chr{str(chr[i])}/")
params_df = pd.DataFrame(data={
    '--input ids': [f"{bucket}/{gwas_dir}/step1_input/ehr_IDs.txt" for _ in range(22)],
    '--input bgenFile': [f"{exome_bgen_dir}/exome.chr{x+1}.bgen" for x in range(22)],
    '--input bgenFileIndex': [f"{exome_bgen_dir}/exome.chr{x+1}.bgen.bgi" for x in range(22)],
    '--input-recursive groupFile': [f"{bucket}/{gwas_dir}/step2_input/" for _ in range(22)],
    '--env maf': mafs,
    '--env CHR': chr,
    '--env CHROM': ['chr'+str(x+1) for x in range(22)],
    '--output-recursive OUT_DIR': [f"{bucket}/{gwas_dir}/step3_output" for i in range(22)],
    '--env annotation_in_groupTest': ['missense;lof;synonymous']*22,
    '--env anno': anno
})

PARAMETER_FILENAME = f'{JOB_NAME}_params.tsv'
%env PARAMETER_FILENAME={PARAMETER_FILENAME}

params_df.to_csv(PARAMETER_FILENAME, sep='\t', index=False)

job_id = !source ~/aou_dsub.bash; aou_dsub \
  --name "${JOB_NAME}" \
  --provider google-batch \
  --regions us-central1 \
  --image "wzhou88/saige:1.3.0" \
  --logging "${WORKSPACE_BUCKET}/dsub_logs/saige_step3" \
  --boot-disk-size 50 \
  --disk-size 100 \
  --min-ram 10 \
  --tasks "${PARAMETER_FILENAME}" \
  --script "Cmd_step3.sh"
print("\n".join(job_id))
job_id = job_id[1].split(" ")[-1]
%env JOB_ID={job_id}

```
<br>

Step3 generates a sparse LD matrix for each gene using the specified gene_file_prefix. This prefix corresponds to the gene-based LD matrix files and matches the prefix of the marker_info.txt file. Additionally, Step3 produces the `marker_info.txt` file, which contains variant information within the LD matrix and serves as an input for Meta-SAIGE. Examples are provided in the `extdata/test_input` directory.

## Running Meta-Analysis using META-SAIGE
This section describes running META-SAIGE. To begin with, we are ready with sparse LD generated files and single variant level summary for two cohorts.
By providing groupfiles, one can generate all masks along with 
```
cmd_line = '''
INPUT_DIR1=${BUCKET}/SAIGE_GENE/imported/step3_docker/step3_docker/WES470k_${maf}_${anno}_chr${i}/
INPUT_DIR2=${BUCKET}/SAIGE_GENE/step3_output/${maf}_${anno}_chr${i}/
INPUT_GWAS1=${BUCKET}/SAIGE_GENE/imported/step2_phenome/step2/WES470k_Whites_${Phecode}_chr${i}
INPUT_GWAS2=${BUCKET}/SAIGE_GENE/step2_output/Pheno_${Phecode}_chr${i}_step2_output_single.txt

chrom=$i
cd /app/
Rscript ./inst/scripts/RV_meta_GC.R \
    --num_cohorts 4 \
    --chr ${chrom} \
    --col_co 10 \
    --info_file_path /app/tmp_chr${chrom}_mkr_info1.txt $INPUT_DIR2/${maf}_${anno}_chr${chrom}_loftee.marker_info.txt \
    --gene_file_prefix $INPUT_DIR1/WES470k_chr${chrom}_ $INPUT_DIR2/${maf}_${anno}_chr${chrom}_loftee_ \
    --gwas_path /app/ukb_dr_chr${i}_GWAS1.txt $INPUT_GWAS2 \
    --groupFile=$groupFile/UKB470k_chr_${CHR}_groupfile.txt \
    --trait_type binary \
    --verbose FALSE \
    --output_prefix $OUT_DIR/META_SAIGE_${phen}_${maf}_${anno}_chr${chrom}
'''
with open("Cmd_meta_saige_0.3.0_MA.sh", "w") as text_file:
    text_file.write(cmd_line)

USER_NAME = os.getenv('OWNER_EMAIL').split('@')[0].replace('.','-')
%env USER_NAME={USER_NAME}
#bucket = os.getenv('WORKSPACE_BUCKET')
bucket="gs://fc-secure-edebf927-ae0c-43f6-87a4-3bae475b5e9f"
gwas_dir = 'SAIGE_GENE'
maf=(['0.01']*22+['0.001']*22+['0.0001']*22)*6
anno=(['lof']*66+['missense_lof']*66+['missense_lof_synonymous']*66)*2
phen = ['T2D']*198+['colorectal']*198
phen_ukb = ['t2d']*198+['colca']*198
phecodes = ['250.2']*198+['153']*198
chrom=[i for i in range(1,23)]*3*3*2
tot_len = len(maf)
params_df0 = pd.DataFrame(data={
    '--env maf': maf,
    '--env phen': phen,
    '--env phen_ukb': phen_ukb,
    '--output-recursive OUT_DIR': [f"{bucket}/{gwas_dir}/step4_output_recode/" for _ in range(tot_len)],
    '--env Phecode': phecodes,
    '--env anno': anno,
    '--env i': chrom,    
})

print(params_df0.shape)
params_df0.to_csv('saige-step4-seokho92_ma.tsv',sep = '\t', index = False)
%env PARAMETER_FILENAME=saige-step4-seokho92_ma.tsv
%env JOB_NAME=Meta_MA
job_id = !source ~/aou_dsub.bash; aou_dsub \
  --name "${JOB_NAME}" \
  --provider google-batch \
  --image "gcr.io/pheweb/meta-saige:0.3.0" \
  --logging "gs://fc-secure-edebf927-ae0c-43f6-87a4-3bae475b5e9f/dsub_logs/saige_step4" \
  --tasks saige-step4-seokho92_ma.tsv \
  --mount BUCKET="gs://fc-secure-edebf927-ae0c-43f6-87a4-3bae475b5e9f" \
  --script "Cmd_meta_saige_0.3.0_MA.sh" \
  --min-ram 5 \
  --disk-size 20
job_id2 = job_id[1].split(" ")[-1]
%env JOB_ID={job_id2}
job_id
```
### Parameter Description
- `n.cohorts` : number of cohorts
- `chr` : chromosome number
- `gwas_path` : path to the GWAS summary. Need to specify GWAS summary file from each and every cohort delimited by white-space (`' '`)
- `info_path` : path to the marker_info.txt file generated from SAIGE 'step3_LDmat.R'. Need to specify marker_info.txt file from each and every cohort delimited by white-space (`' '`)
- `gene_file_prefix` : prefix to the LD matrix separated by genes (also generated from SAIGE `step3_LDmat.R`) usually same as marker_info.txt file's prefix
- `col_co` : ultra-rare variant collapsing cut-off. (default is 10)
- `output_path` : path for the meta-analysis resutls
- `verbose`(optional): verbose mode. TRUE or FALSE
- `ancestry`(optional) :  Ancestry indicator (ex. 1 1 1 2). Need to specify ancestry indicator from each and every cohort delimited by white-space (`' '`). Optional input for multi-ancestry meta-analysis.
- `trait_type` : trait type (binary or continuous)
- `groupfile`(optional) : path to the group file for gene-based analysis (Same format as SAIGE-GENE+)
- `annotation`(optional) : functional annotation for the variants of interest. ex. c('lof', 'missense_lof')
- `mafcutoff`(optional) : Maximum MAF for group-based analysis ex. c(0.01 0.001 0.0001)
