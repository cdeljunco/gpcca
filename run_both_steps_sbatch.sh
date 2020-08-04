kmin=2
kmax=15
klist='[ 8 ]'
echo $klist

#these lines are for splitting up the vars below by telling bash that the separator is a comma
OLDIFS=$IFS 
IFS=','
#bash implementation of zip(list1, list2) in python
for vars in 0.1,1.53e-01 0.5,9.86e-02 1,5.23e-02 5,2.96e-02 10,1.24e-02 50,2.70e-03 1000,1.00e-04
do set -- $vars;
ratio=$1
lag=$2

for nsims in 2 3 4 6 8 11 15 20 27 37 49 66 88 118 159 213 285 382 512
do

let y=1024/nsims

x=$((y<100 ? y : 100))

for ((i=10;i<x;i++));

do 

matrix='minimal-'$ratio'-'$i'-nsims'$nsims'-nclusters16-tc-lag'$lag'-interval1.00e-01'

#echo $matrix
mkdir "../G-PCCA-Results/$matrix"

cp "Count_Matrices/$matrix.txt"  "../G-PCCA-Results/$matrix/"

#write sbatch files

#echo "writing first run file"

sed -e 's/MATRIX/'$matrix'/g' -e 's/KMIN/'$kmin'/' -e 's/KMAX/'$kmax'/' <Templates/run-step1-template.sh > run-step1-"$matrix".sh

#echo "writing second run file"

sed -e 's/MATRIX/'$matrix'/g' -e "s/KLIST/$klist/" <Templates/run-step2-template.sh > run-step2-"$matrix".sh

#run step 1

jobid=$(sbatch run-step1-"$matrix".sh | cut -d " " -f 4)

#run step 2 only after step 1 is complete

sbatch --dependency=afterok:$jobid run-step2-"$matrix".sh

done
done
done

IFS=$OLDIFS
