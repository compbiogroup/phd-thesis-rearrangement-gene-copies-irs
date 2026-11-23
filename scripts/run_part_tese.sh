# Run the partition algorithms in partition/ (MCFISP-exe) on the thesis data set
# (db/, all the files of N200 and N500), and turn each
# partition into a reduced instance (each block becomes a single gene, see
# reduce_partition.py) that can be given to dec/dist of cyc.
#
# MCFISP-exe writes the result of each instance of an input file to the
# temporary folder out_folder/ALL/<ALGO>/N<n>/<mod>/<file name without .txt>/.
# When it finishes, three files are written for the input file:
#   out_folder/PART/<ALGO>/N<n>/<mod>/<file name>.txt  the outputs of MCFISP-exe
#       joined (per instance: the partition of the origin and of the target, with
#       the blocks separated by |, the compatible blocks, and the time);
#   out_folder/COST/<ALGO>/N<n>/<mod>/<file name>.txt  one line per instance:
#       breakpoints_origin breakpoints_target cost blocks_origin blocks_target
#       exclusive_origin exclusive_target time   (cost = number of breakpoints of
#       the two partitions, the cost of the thesis);
#   out_folder/RED/<ALGO>/N<n>/<mod>/<file name>.txt   the reduced instances, in
#       the format of the thesis data set (4 lines per instance, genomes without
#       caps), with the same folder structure as input_folder.
# For Approx there is also
#   out_folder/PHASE2/Approx/N<n>/<mod>/<file name>.txt one line per instance: the
#       number of times the second stage of the approximation (the one that breaks a
#       block D) was used, or NA if the instance failed (MCFISP-exe writes it to
#       <output prefix>_phase2).
# and the temporary folder is deleted. Input files whose results already exist
# are skipped, so the script can be stopped and started again.
#
# The variables INPUT_FOLDER, OUT_FOLDER, N_INST and ALGOS override the defaults below.
#
# To run the heuristics of cyc on the reduced instances of an algorithm ALGO, use
# the option PART of run_dec_tese.sh and run_dist_tese.sh, which take
# out_folder/RED/<ALGO> as the input folder and write to TESE_OUT/ANS_TESE_DEC_<ALGO>
# and TESE_OUT/ANS_TESE_DIST_<ALGO>, e.g.
#   PART=SOARComb ./run_dist_tese.sh map
# The reduced genomes have a single copy of each gene, so a single mapping is enough
# (see scripts/PARTICAO_CYC.md).

scripts=$(cd "$(dirname "$0")" && pwd)
cd "$scripts/../partition"

# build MCFISP-exe (copied to the partition folder)
stack install --local-bin-path . > /dev/null || exit 1

input_folder=${INPUT_FOLDER:-$scripts/../db}
out_folder=${OUT_FOLDER:-$scripts/../TESE_OUT/ANS_TESE_PART}

sizes="200 500"
mods="SDCJ SREV SREV_TRANS TRANS"
# Partition algorithms (option -a of MCFISP-exe): Greedy, GreedySin, SOAR,
# SOARComb, Approx and FPT. FPT is exponential, so it is not in the default list
# (neither is Approx, which used to take about 8 minutes per instance of N200).
algos=${ALGOS:-"Greedy GreedySin SOAR SOARComb"}
algo_params="-s" # the genomes of the data set are signed
n_inst=${N_INST:-100} # instances per input file

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
                red=$out_folder/RED/${sub}/${name}.txt
                phase2=$out_folder/PHASE2/${sub}/${name}.txt
                if [ -f $part ] && [ -f $cost ] && [ -f $red ] && { [ $algo != Approx ] || [ -f $phase2 ]; }
                then
                    continue
                fi
                rm -rf $out
                mkdir -p $out $out_folder/PART/${sub} $out_folder/COST/${sub} $out_folder/RED/${sub}
                [ $algo = Approx ] && mkdir -p $out_folder/PHASE2/${sub}
                echo "$algo $name"
                ./MCFISP-exe -a $algo $algo_params -i $file -o $out/inst

                # join only complete results (written to .tmp and renamed, so an
                # interrupted join is not taken as complete)
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
                if grep -q "^# Error\|^# No solution" $outs
                then
                    echo "some instances of $algo $name have no partition (NA in $cost)"
                fi
                python3 $scripts/reduce_partition.py $red.tmp $cost.tmp $outs || {
                    echo "reduction failed for $algo $name, keeping $out"
                    rm -f $red.tmp $cost.tmp
                    continue
                }
                cat $outs > $part.tmp && mv $part.tmp $part
                mv $cost.tmp $cost
                mv $red.tmp $red
                [ $algo = Approx ] && mv $out/inst_phase2 $phase2
                rm -r $out
            done
        done
    done
done

# remove the empty folders left in the temporary folder
find $out_folder/ALL -type d -empty -delete 2>/dev/null

cd - > /dev/null
