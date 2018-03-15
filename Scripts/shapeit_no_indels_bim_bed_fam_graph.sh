#!/bin/bash

#SBATCH -J shapeit

#SBATCH -A MRC-BSU-SL2

#SBATCH --nodes 1

#SBATCH --ntasks 16

#SBATCH --time 2:00:00

#SBATCH --mail-type FAIL

#SBATCH -p  mrc-bsu-sand

#SBATCH --output=shape.graph.err

. /etc/profile.d/modules.sh # Leave this line (enables the module command)

module purge                # Removes all modules still loaded

module load default-impi    # REQUIRED - loads the basic environment


export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets

#JOBID=$SLURM_JOB_ID
#echo -e "JobID: $JOBID

echo "Time: `date`

echo "Running on master node: `hostname`


export swP3='/scratch/wallace/1000GP_Phase3'
export mseas='/mrc-bsu/scratch/ev250/ase-il2ra/shapeit'
export mseao='/mrc-bsu/scratch/ev250/ase-il2ra/objects'

cd $mseas
## only use EUR ancestry from ref panel and exclude fake sample from bed, bim , fam files
##echo "EUR" > group.list
##echo "fake" > exclude.inds 

# $swlb/shapeit -check \
#     --input-bed $mseao/genoB37no_I_D_hetsample.bed $mseao/genoB37no_I_D_hetsample.bim$mseao/genoB37no_I_D_hetsample.bim $mseao/genoB37no_I_D_hetsample.fam  \
#     -M $swP3/genetic_map_chr10_combined_b37.txt \
#     --input-ref $swP3/1000GP_Phase3_chr10.hap.gz $swP3/1000GP_Phase3_chr10.legend.gz $swP3/1000GP_Phase3.sample \
#     --exclude-ind exclude.inds \
#     --output-log chr10.aligments3


$swlb/shapeit --input-bed $mseao/genoB37no_I_D_hetsample.bed $mseao/genoB37no_I_D_hetsample.bim $mseao/genoB37no_I_D_hetsample.fam  \
    -M $swP3/genetic_map_chr10_combined_b37.txt \
    --input-ref $swP3/1000GP_Phase3_chr10.hap.gz $swP3/1000GP_Phase3_chr10.legend.gz $swP3/1000GP_Phase3.sample \
    --exclude-snp $mseas/chr10.aligments3.snp.strand.exclude \
    --exclude-ind exclude.inds \
    --include-grp group.list \
    --output-graph il2ra.bbfam.graph \
    -T 16

## phasing uncertainty with selected SNPs (B37 position)
echo "hap.uncertainty 6053374 6094697 6124257" > hap.snps

$swlb/shapeit -convert \
	      --input-graph il2ra.bbfam.graph \
	      --input-sets hap.snps
