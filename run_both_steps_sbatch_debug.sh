kmin=2
kmax=15
klist='[ 4 8 ]'
echo $klist

matrix='filedoesnotexist'

mkdir "Results/$matrix"

cp "Count_Matrices/$matrix.txt"  "Results/$matrix/"

#write sbatch files

echo "writing first run file"

sed -e 's/MATRIX/'$matrix'/g' -e 's/KMIN/'$kmin'/' -e 's/KMAX/'$kmax'/' <Templates/run-step1-template.sh > run-step1-"$matrix".sh

echo "writing second run file"

sed -e 's/MATRIX/'$matrix'/g' -e "s/KLIST/$klist/" <Templates/run-step2-template.sh > run-step2-"$matrix".sh

#run step 1

jobid=$(sbatch run-step1-"$matrix".sh | cut -d " " -f 4)

#run step 2 only after step 1 is complete

sbatch --dependency=afterok:$jobid run-step2-"$matrix".sh
