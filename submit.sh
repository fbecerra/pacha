#!/bin/bash

if [ $# -lt 2 ]; then
    echo Not enough arguments! Exiting...
    echo Usage: ./submit.sh script base ncores
    exit
fi

tmp=submit_tmp.sh
homes=/scratch/02563/fbecerra/paha

script=$1
base=$2
ncores=$3

rm -f $tmp

echo '#!/bin/bash' > $tmp

echo '#SBATCH -A TG-AST130020' >> $tmp
echo "#SBATCH -n $ncores" >> $tmp
echo "#SBATCH -J $base" >> $tmp
echo '#SBATCH -p normal' >> $tmp
echo '#SBATCH -t 24:00:00' >> $tmp
echo "#SBATCH -o $homes/output_files/output_%j.out" >> $tmp

echo "ibrun tacc_affinity python $script > $homes/output_files/outlog_$base.out" >> $tmp

sbatch $tmp

rm -f $tmp
