# Run dec on a subset of the thesis data set (db/, see db/README.md).
# Files are input_folder/N<n>/<mod>/N<n>_<mod><op>_I<indels>_L<alphabet>_F<flex>.txt
# with 100 instances each (4 lines per instance, genomes without caps).
#
# dec writes the results of each input file to the temporary folder
# out_folder/ALL/<HEUR>_DEC_<algo>/N<n>/<mod>/<file name without .txt>/
# (one file per instance). When it finishes, they are joined (as in join_dcj.sh)
# into out_folder/DCJ/<HEUR>_DEC_<algo>/N<n>/<mod>/<file name>.txt (the -dcj files)
# and out_folder/DEC/... (the -best files), and the temporary folder is deleted.
# Input files whose joined results already exist are skipped, so the script can
# be stopped and started again.

# The variables INPUT_FOLDER, OUT_FOLDER, K and N_INST override the
# defaults below.
#
# Only the input files with the values of op (conservative operations), I
# (indels), L (alphabet) and F (flexibility) in OP_VALS, I_VALS, L_VALS and
# F_VALS (separated by spaces) are used. The defaults are the subset of the
# data set used in the thesis: op in 100 200, I in 0 50, L in 10 100 200 and
# F in 0 40 100 (36 files per operation). Set a variable to change its values,
# e.g. OP_VALS=100 I_VALS=0 ./run_dec_tese.sh, or to an empty value to keep all
# the values of its parameter, e.g. F_VALS= ./run_dec_tese.sh.
#
# With PART=<ALGO> (a partition algorithm of run_part_tese.sh, e.g. SOARComb),
# dec runs on the reduced instances of that algorithm
# (part_folder/RED/<ALGO>, part_folder is PART_FOLDER or TESE_OUT/ANS_TESE_PART)
# and the results go to TESE_OUT/ANS_TESE_DEC_<ALGO>. The reduced genomes have no replicas, so
# the cycle packing is determined and only rand with one iteration is used.
# Instances without a partition (lines NA) get the result NA.
#
# With FIX=<ALGO>, dec runs on the original instances with the adjacencies of
# the partitions of that algorithm fixed (option -f of dec with the files
# part_folder/FIX/<ALGO>/N<n>/<mod>/<file name>.txt of run_part_fix_tese.sh),
# and the results go to TESE_OUT/ANS_TESE_DEC_FIX_<ALGO>. With FIX_PAIRS=1 the
# correspondence of genes of the partitions is also fixed (option -p), the
# results go to TESE_OUT/ANS_TESE_DEC_FIX_<ALGO>_PAIRS and, since the packing is then
# determined, only rand with one iteration is used.

cd "$(dirname "$0")/../cyc"

make

input_folder=${INPUT_FOLDER:-../db}
out_folder=${OUT_FOLDER:-../TESE_OUT/ANS_TESE_DEC}

sizes="200"
mods="SDCJ SREV SREV_TRANS"
heurs="ga rand"
algos="bfs" # bfs, map or 1c2c3c
k=${K:-500} # iterations of dec
ga_params="--mut 20 -c 50"
n_inst=${N_INST:-100} # instances per input file

# subset of the data set (an empty list keeps all the values), the same in
# run_dist_tese.sh
op_vals=${OP_VALS-"100 200"}
i_vals=${I_VALS-"0 50"}
l_vals=${L_VALS-"10 100 200"}
f_vals=${F_VALS-"0 40 100"}

# Whether value is in the list (an empty list has every value)
selected () {
    local value=$1 list=$2
    [ -z "$list" ] || [[ " $list " == *" $value "* ]]
}

if [ -n "$PART" ]
then
    part_folder=${PART_FOLDER:-../TESE_OUT/ANS_TESE_PART} # out_folder of run_part_tese.sh
    input_folder=${INPUT_FOLDER:-$part_folder/RED/$PART}
    out_folder=${OUT_FOLDER:-../TESE_OUT/ANS_TESE_DEC_$PART}
    heurs="rand"
    k=1
fi

fix_params=""
if [ -n "$FIX" ]
then
    part_folder=${PART_FOLDER:-../TESE_OUT/ANS_TESE_PART} # out_folder of run_part_fix_tese.sh
    suffix=FIX_$FIX
    if [ -n "$FIX_PAIRS" ]
    then
        suffix=${suffix}_PAIRS
        fix_params="-p"
        heurs="rand"
        k=1
    fi
    out_folder=${OUT_FOLDER:-../TESE_OUT/ANS_TESE_DEC_$suffix}
fi

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
            for heur in $heurs
            do
                for algo in $algos
                do
                    sub=${heur^^}_DEC_${algo}/N${n}/${mod}
                    out=$out_folder/ALL/${sub}/${name}
                    dcj=$out_folder/DCJ/${sub}/${name}.txt
                    dec=$out_folder/DEC/${sub}/${name}.txt
                    if [ -f $dcj ] && [ -f $dec ]
                    then
                        continue
                    fi
                    rm -rf $out
                    mkdir -p $out $out_folder/DCJ/${sub} $out_folder/DEC/${sub}
                    params=""
                    if [ $heur == "ga" ]
                    then
                        params=$ga_params
                    fi
                    if [ -n "$FIX" ]
                    then
                        fix_file=$part_folder/FIX/$FIX/N${n}/${mod}/${name}.txt
                        if [ ! -f $fix_file ]
                        then
                            echo "missing $fix_file, skipping"
                            continue
                        fi
                        params="$params -f $fix_file $fix_params"
                    fi
                    echo "$heur $algo $name"
                    ./dec $heur $algo -e -i $file -o $out -k $k $params

                    # join only complete results (written to .tmp and renamed,
                    # so an interrupted join is not taken as complete)
                    if [ $(ls $out | grep -c -- "-best$") -ne $n_inst ] ||
                       [ $(ls $out | grep -c -- "-dcj$") -ne $n_inst ]
                    then
                        echo "incomplete results for $heur $algo $name, keeping $out"
                        continue
                    fi
                    cat $out/*-dcj > $dcj.tmp && mv $dcj.tmp $dcj
                    cat $out/*-best > $dec.tmp && mv $dec.tmp $dec
                    rm -r $out
                done
            done
        done
    done
done

# remove the empty folders left in the temporary folder
find $out_folder/ALL -type d -empty -delete 2>/dev/null

cd -
