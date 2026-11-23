# Write the partitions of the algorithms in partition/ (MCFISP-exe) for the thesis
# data set in positions of the original genomes, without reducing them (see
# partition_fixed.py and PARTICAO_CYC.md), so that cyc can build the adjacency
# graph with fixed adjacencies (and fixed characters) of the thesis.
#
# For each input file and algorithm, the joined outputs of MCFISP-exe
# out_folder/PART/<ALGO>/N<n>/<mod>/<file name>.txt are reused when they exist
# (they are the same files written by run_part_tese.sh); otherwise MCFISP-exe is
# run, as in run_part_tese.sh, and the PART file is written. Then:
#   out_folder/FIX/<ALGO>/N<n>/<mod>/<file name>.txt   3 lines per instance:
#       breakpoints of the origin, breakpoints of the target (indices of
#       intergenic regions of the input genomes) and the correspondence of genes;
#   out_folder/COST/<ALGO>/N<n>/<mod>/<file name>.txt  the costs, as in
#       run_part_tese.sh;
#   out_folder/PHASE2/Approx/N<n>/<mod>/<file name>.txt only for Approx, the number
#       of times the second stage of the approximation was used in each instance, as
#       in run_part_tese.sh. MCFISP-exe is run again for Approx when the PART file
#       exists but this one does not, since it is only written by MCFISP-exe.
# The instances of FIX are in the same order as the instances of the input file
# in input_folder, which is still the input of cyc. Input files whose results
# already exist are skipped, so the script can be stopped and started again.
#
# The variables INPUT_FOLDER, OUT_FOLDER, N_INST and ALGOS override the defaults below.

scripts=$(cd "$(dirname "$0")" && pwd)
cd "$scripts/../partition"

input_folder=${INPUT_FOLDER:-$scripts/../db}
out_folder=${OUT_FOLDER:-$scripts/../TESE_OUT/ANS_TESE_PART}

sizes="200"
mods="SDCJ SREV SREV_TRANS"
# Partition algorithms (option -a of MCFISP-exe), see run_part_tese.sh
algos=${ALGOS:-"Greedy SOARComb Approx"}
algo_params="-s" # the genomes of the data set are signed
n_inst=${N_INST:-100} # instances per input file

built=false
for n in $sizes
do
    for mod in $mods
    do
        for file in $input_folder/N${n}/${mod}/N${n}_${mod}[0-9]*_I*_L*_F*.txt
        do
            [ -f "$file" ] || continue # no file with this size/operation
            name=$(basename $file .txt)
            for algo in $algos
            do
                sub=${algo}/N${n}/${mod}
                out=$out_folder/ALL/${sub}/${name}
                part=$out_folder/PART/${sub}/${name}.txt
                cost=$out_folder/COST/${sub}/${name}.txt
                fix=$out_folder/FIX/${sub}/${name}.txt
                phase2=$out_folder/PHASE2/${sub}/${name}.txt
                need_phase2=false
                [ $algo = Approx ] && [ ! -f $phase2 ] && need_phase2=true
                if [ -f $fix ] && [ -f $cost ] && ! $need_phase2
                then
                    continue
                fi
                mkdir -p $out_folder/PART/${sub} $out_folder/COST/${sub} $out_folder/FIX/${sub}
                [ $algo = Approx ] && mkdir -p $out_folder/PHASE2/${sub}
                if [ ! -f $part ] || $need_phase2
                then
                    if ! $built
                    then
                        # build MCFISP-exe (copied to the partition folder)
                        stack install --local-bin-path . > /dev/null || exit 1
                        built=true
                    fi
                    rm -rf $out
                    mkdir -p $out
                    echo "$algo $name"
                    ./MCFISP-exe -a $algo $algo_params -i $file -o $out/inst
                    outs=$(ls $out/inst_[0-9]* 2> /dev/null | sort)
                    if [ $(echo "$outs" | grep -c .) -ne $n_inst ]
                    then
                        echo "incomplete results for $algo $name, keeping $out"
                        continue
                    fi
                    if [ $algo = Approx ] && [ $(grep -c . $out/inst_phase2 2> /dev/null || echo 0) -ne $n_inst ]
                    then
                        echo "missing or incomplete $out/inst_phase2 for $algo $name, keeping $out"
                        continue
                    fi
                    cat $outs > $part.tmp && mv $part.tmp $part
                    [ $algo = Approx ] && mv $out/inst_phase2 $phase2
                    rm -r $out
                fi
                if [ $(grep -c "^# Time" $part) -ne $n_inst ]
                then
                    echo "$part does not have $n_inst instances, skipping"
                    continue
                fi
                python3 $scripts/partition_fixed.py $fix.tmp $cost.tmp $part $file || {
                    echo "conversion failed for $algo $name"
                    rm -f $fix.tmp $cost.tmp
                    continue
                }
                mv $cost.tmp $cost
                mv $fix.tmp $fix
            done
        done
    done
done

# remove the empty folders left in the temporary folder
find $out_folder/ALL -type d -empty -delete 2>/dev/null

cd - > /dev/null
