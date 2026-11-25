# Run the FPT of partition/ (MCFISP-exe -a FPT) on the first instances of each file
# of the subset of the thesis data set used by run_dec_tese.sh (db/, see
# db/README.md), together with heuristics for comparison on the same instances.
# Without the variables below, it uses the first 10 instances of each of the 144
# files of N200 with op in 100 200, I in 0 50, L in 10 100 200 and F in 0 40 100
# (the four operations, as in the smaller data set of the thesis).
#
# The FPT starts with the best partition of Greedy and SOARComb, stops after
# TIME_LIMIT seconds per instance (default 600) and then writes the best partition
# found, marked as "# Partial solution". MCFISP-exe runs the instances of a file in
# parallel, so each input file takes at most about TIME_LIMIT seconds with the FPT.
#
# For each input file and algorithm, MCFISP-exe writes the result of each instance
# to the temporary folder out_folder/ALL/<ALGO>/N<n>/<mod>/<file name without .txt>/.
# When it finishes, the following files are written (one line per instance in the
# order of the input file):
#   out_folder/PART/<ALGO>/N<n>/<mod>/<file name>.txt  the outputs of MCFISP-exe
#       joined;
#   out_folder/COST/<ALGO>/N<n>/<mod>/<file name>.txt  breakpoints_origin
#       breakpoints_target cost blocks_origin blocks_target exclusive_origin
#       exclusive_target time (as in run_part_tese.sh), or NA. They are computed
#       by partition_fixed.py, which, unlike reduce_partition.py, accepts circular
#       chromosomes (the SDCJ files);
#   out_folder/STATUS/FPT/N<n>/<mod>/<file name>.txt  only for FPT, Exact if the
#       search finished, Partial if it reached the time limit, or NA.
# and the temporary folder is deleted. Input files whose results already exist
# are skipped, so the script can be stopped and started again. The table of the
# costs is made with
#   python3 table_part_tese.py ../TESE_OUT/ANS_TESE_FPT
#
# Each call of MCFISP-exe also has an external time limit (enough for its instances
# with TIME_LIMIT seconds each), since in rare cases it stopped without finishing
# some instances. The results of such a file are incomplete and kept in the
# temporary folder, and the file is run again by the next execution of the script.
#
# The variables INPUT_FOLDER, OUT_FOLDER, SIZES, N_INST (instances taken from each
# file), ALGOS, MODS and TIME_LIMIT override the defaults below. OP_VALS, I_VALS,
# L_VALS and F_VALS select the input files as in run_dec_tese.sh (an empty value
# keeps all the values of the parameter).

scripts=$(cd "$(dirname "$0")" && pwd)
cd "$scripts/../partition"

# build MCFISP-exe (copied to the partition folder)
stack install --local-bin-path . > /dev/null || exit 1

input_folder=${INPUT_FOLDER:-$scripts/../db}
out_folder=${OUT_FOLDER:-$scripts/../TESE_OUT/ANS_TESE_FPT}

sizes=${SIZES:-"200"}
mods=${MODS:-"SDCJ SREV SREV_TRANS TRANS"}
# Partition algorithms (option -a of MCFISP-exe). The heuristics come first, so
# that their results are available while the FPT runs.
algos=${ALGOS:-"Greedy SOARComb Approx FPT"}
algo_params="-s" # the genomes of the data set are signed
n_inst=${N_INST:-10} # instances taken from the start of each input file
time_limit=${TIME_LIMIT:-600} # seconds per instance (only for the FPT)
# external limit of a call of MCFISP-exe, which runs nproc instances at a time
call_limit=$(( (n_inst / $(nproc) + 2) * time_limit + 600 ))

# subset of the data set (an empty list keeps all the values), the same in
# run_dec_tese.sh
op_vals=${OP_VALS-"100 200"}
i_vals=${I_VALS-"0 50"}
l_vals=${L_VALS-"10 100 200"}
f_vals=${F_VALS-"0 40 100"}

# Whether value is in the list (an empty list has every value)
selected () {
    local value=$1 list=$2
    [ -z "$list" ] || [[ " $list " == *" $value "* ]]
}

for n in $sizes
do
    for mod in $mods
    do
        for file in $input_folder/N${n}/${mod}/N${n}_${mod}[0-9]*_I*_L*_F*.txt
        do
            [ -f "$file" ] || continue # no file with this size/operation
            name=$(basename $file .txt)
            # values of op, I, L and F of the file
            [[ $name =~ ^N${n}_${mod}([0-9]+)_I([0-9]+)_L([0-9]+)_F([0-9]+)$ ]] ||
                continue
            selected ${BASH_REMATCH[1]} "$op_vals" &&
                selected ${BASH_REMATCH[2]} "$i_vals" &&
                selected ${BASH_REMATCH[3]} "$l_vals" &&
                selected ${BASH_REMATCH[4]} "$f_vals" || continue
            for algo in $algos
            do
                sub=${algo}/N${n}/${mod}
                out=$out_folder/ALL/${sub}/${name}
                part=$out_folder/PART/${sub}/${name}.txt
                cost=$out_folder/COST/${sub}/${name}.txt
                status=$out_folder/STATUS/${sub}/${name}.txt
                if [ -f $part ] && [ -f $cost ] && { [ $algo != FPT ] || [ -f $status ]; }
                then
                    continue
                fi
                rm -rf $out
                mkdir -p $out $out_folder/PART/${sub} $out_folder/COST/${sub}
                [ $algo = FPT ] && mkdir -p $out_folder/STATUS/${sub}
                # first n_inst instances (4 lines each) of the input file
                head -n $((4 * n_inst)) $file > $out/input.txt
                echo "$algo $name"
                timeout $call_limit ./MCFISP-exe -a $algo $algo_params -t $time_limit -i $out/input.txt -o $out/inst ||
                    echo "MCFISP-exe failed or reached the external time limit for $algo $name"

                # join only complete results (written to .tmp and renamed, so an
                # interrupted join is not taken as complete)
                outs=$(ls $out/inst_[0-9]* 2> /dev/null | sort)
                if [ $(echo "$outs" | grep -c .) -ne $n_inst ]
                then
                    echo "incomplete results for $algo $name, keeping $out (run the script again)"
                    continue
                fi
                if grep -q "^# Error\|^# No solution" $outs
                then
                    echo "some instances of $algo $name have no partition (NA in $cost)"
                fi
                cat $outs > $part.tmp
                python3 $scripts/partition_fixed.py $out/fix.txt $cost.tmp $part.tmp $out/input.txt || {
                    echo "cost computation failed for $algo $name, keeping $out"
                    rm -f $part.tmp $cost.tmp
                    continue
                }
                if [ $algo = FPT ]
                then
                    for o in $outs
                    do
                        if grep -q "^# Exact solution" $o; then echo Exact
                        elif grep -q "^# Partial solution" $o; then echo Partial
                        else echo NA
                        fi
                    done > $status.tmp
                    grep -q Partial $status.tmp &&
                        echo "$(grep -c Partial $status.tmp) instances of $name reached the time limit"
                fi
                mv $part.tmp $part
                mv $cost.tmp $cost
                [ $algo = FPT ] && mv $status.tmp $status
                rm -r $out
            done
        done
    done
done

# remove the empty folders left in the temporary folder
find $out_folder/ALL -type d -empty -delete 2>/dev/null

cd - > /dev/null
