# Run dist on a subset of the thesis data set (db/, the same input files as
# run_dec_tese.sh).
#
# usage: run_dist_tese.sh map|dec
#   map - map the genes with random mappings (k_map iterations per instance)
#         and keep the best one.
#   dec - map the genes with the decompositions found by run_dec_tese.sh
#         (dec_folder/DEC/<HEUR>_DEC_<algo>/N<n>/<mod>/<file name>.txt, the
#         joined -best files), k_dec iterations per instance.
#
# dist writes the results of each input file to the temporary folder
# out_folder/ALL/<SUB>/N<n>/<mod>/<file name without .txt>/ (one file per
# instance), where SUB is MAP_<dist_algo> or <HEUR>_DEC_<algo>_<dist_algo>.
# When it finishes, the -best files are joined into
# out_folder/DIST/<SUB>/N<n>/<mod>/<file name>.txt and the temporary folder is
# deleted. Input files whose joined results already exist are skipped, so the
# script can be stopped and started again.

option=$1
if [ "$option" != "map" ] && [ "$option" != "dec" ]
then
    echo "usage: $0 map|dec"
    exit 1
fi

# The variables INPUT_FOLDER, OUT_FOLDER, DEC_FOLDER, K_MAP and N_INST override the
# defaults below. OP_VALS, I_VALS, L_VALS and F_VALS select the input files, with
# the same defaults (the subset of the thesis) as in run_dec_tese.sh.
#
# With PART=<ALGO> (a partition algorithm of run_part_tese.sh, e.g. SOARComb),
# dist runs on the reduced instances of that algorithm
# (part_folder/RED/<ALGO>, part_folder is PART_FOLDER or TESE_OUT/ANS_TESE_PART),
# with the decompositions of run_dec_tese.sh with the same PART
# (TESE_OUT/ANS_TESE_DEC_<ALGO>), and the results go to
# TESE_OUT/ANS_TESE_DIST_<ALGO>. The
# reduced genomes have no replicas, so a single mapping is used. Instances
# without a partition (lines NA) get the result NA.
#
# With FIX=<ALGO> (and FIX_PAIRS=1), only the option dec can be used: dist
# runs on the original instances with the decompositions of run_dec_tese.sh with
# the same FIX (and FIX_PAIRS), and the results go to
# TESE_OUT/ANS_TESE_DIST_FIX_<ALGO> (or TESE_OUT/ANS_TESE_DIST_FIX_<ALGO>_PAIRS).

cd "$(dirname "$0")/../cyc"

make

input_folder=${INPUT_FOLDER:-../db}
dec_folder=${DEC_FOLDER:-../TESE_OUT/ANS_TESE_DEC} # out_folder of run_dec_tese.sh
out_folder=${OUT_FOLDER:-../TESE_OUT/ANS_TESE_DIST}

# must be the same values used in run_dec_tese.sh
sizes="200"
mods="SDCJ SREV SREV_TRANS"
heurs="ga rand"
algos="bfs"
dec_params="-e" # options -e, -z and -u used to run dec

# distance algorithm used for each type of conservative operation
declare -A dist_algos=(
    [SDCJ]=dcj
    [SREV]=reversal_indel25
    [SREV_TRANS]=reversal_transposition_indel6
    [TRANS]=none
)
k_map=${K_MAP:-1000}
k_dec=1
n_inst=${N_INST:-100} # instances per input file

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

if [ -n "$PART" ]
then
    part_folder=${PART_FOLDER:-../TESE_OUT/ANS_TESE_PART} # out_folder of run_part_tese.sh
    input_folder=${INPUT_FOLDER:-$part_folder/RED/$PART}
    dec_folder=${DEC_FOLDER:-../TESE_OUT/ANS_TESE_DEC_$PART}
    out_folder=${OUT_FOLDER:-../TESE_OUT/ANS_TESE_DIST_$PART}
    heurs="rand" # the same used by run_dec_tese.sh with PART
    k_map=${K_MAP:-1}
fi

if [ -n "$FIX" ]
then
    if [ $option == "map" ]
    then
        echo "FIX can only be used with the option dec"
        exit 1
    fi
    suffix=FIX_$FIX
    if [ -n "$FIX_PAIRS" ]
    then
        suffix=${suffix}_PAIRS
        heurs="rand" # the same used by run_dec_tese.sh with FIX_PAIRS
    fi
    dec_folder=${DEC_FOLDER:-../TESE_OUT/ANS_TESE_DEC_$suffix}
    out_folder=${OUT_FOLDER:-../TESE_OUT/ANS_TESE_DIST_$suffix}
fi

# Run dist on file with extra options params, saving the joined result in
# out_folder/DIST/sub/name.txt.
run () {
    local file=$1 sub=$2 params=$3
    local name=$(basename $file .txt)
    local out=$out_folder/ALL/${sub}/${name}
    local res=$out_folder/DIST/${sub}/${name}.txt
    if [ -f $res ]
    then
        return
    fi
    rm -rf $out
    mkdir -p $out $out_folder/DIST/${sub}
    echo "$sub $name"
    ./dist $dist_algo -F -i $file -o $out $dec_params $params

    # join only complete results (written to .tmp and renamed, so an
    # interrupted join is not taken as complete)
    if [ $(ls $out | grep -c -- "-best$") -ne $n_inst ]
    then
        echo "incomplete results for $sub $name, keeping $out"
        return
    fi
    cat $out/*-best > $res.tmp && mv $res.tmp $res
    rm -r $out
}

for n in $sizes
do
    for mod in $mods
    do
        dist_algo=${dist_algos[$mod]}
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
            if [ $option == "map" ]
            then
                run $file MAP_${dist_algo}/N${n}/${mod} "-k $k_map"
                continue
            fi
            for heur in $heurs
            do
                for algo in $algos
                do
                    dec_sub=${heur^^}_DEC_${algo}/N${n}/${mod}
                    dec=$dec_folder/DEC/${dec_sub}/${name}.txt
                    if [ ! -f $dec ]
                    then
                        echo "missing $dec, skipping"
                        continue
                    fi
                    run $file ${heur^^}_DEC_${algo}_${dist_algo}/N${n}/${mod} \
                        "-k $k_dec -D $dec"
                done
            done
        done
    done
done

# remove the empty folders left in the temporary folder
find $out_folder/ALL -type d -empty -delete 2>/dev/null

cd -
