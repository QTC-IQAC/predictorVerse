#!/bin/bash
#SBATCH --job-name=RFAA
#SBATCH -e %j.err
#SBATCH -o %j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --gres=gpu:1


cwd=$(pwd)


inputs_dir=$cwd/RFAA/inputs
outputs_dir=$cwd/RFAA/outputs


cd /home/ramon/progs/RoseTTAFold-All-Atom  # path to the directory where RFAA is installed

#### Main program ####

start_date=$(date)


for file in $inputs_dir/*.yaml; do
    # Change the relative path for absolute ones
	sed -i "s|\./|$cwd/|g" $file
	
	
	# Split input file into file name and folder
	filename=$(basename "$file")
	foldername=$(dirname "$file")
	
	python -m rf2aa.run_inference --config-name $filename --config-dir "$dir"$foldername
	
	
done


final_date=$(date)

cd $cwd

echo "START DATE: $start_date"
echo "FINAL DATE: $final_date"

