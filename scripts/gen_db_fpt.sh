# Generate the small data set used to evaluate the FPT (db_fpt/ by default), with
# the generator of the thesis data set (DB-exe of partition/, as in gen_db.sh of
# the thesis repository) and the parameters proportional to the genome size n:
#   op (conservative operations)  n/4 and n/2
#   I (indels)                    0 and n/4 (I/2 insertions and I/2 deletions)
#   L (alphabet)                  n/2 and n
#   F (flexibility)               0 (equality relation, for which the FPT is exact)
# The files have the same names and folders as db/
# (N<n>/<mod>/N<n>_<mod><op>_I<I>_L<L>_F0.txt) and K pairs of genomes each.
#
# The variables OUT_FOLDER, SIZES (default "25 50"), MODS and K (default 100)
# override the defaults below. Existing files are not overwritten, since the
# generation is random.

scripts=$(cd "$(dirname "$0")" && pwd)
cd "$scripts/../partition"

# build DB-exe (copied to the partition folder)
stack install --local-bin-path . > /dev/null || exit 1

out_folder=${OUT_FOLDER:-$scripts/../db_fpt}
sizes=${SIZES:-"25 50"}
mods=${MODS:-"SDCJ SREV SREV_TRANS TRANS"}
k=${K:-100} # pairs of genomes per file

for n in $sizes
do
    for mod in $mods
    do
        # p is the percentage of reversals in the conservative operations
        case $mod in
            SREV) p=100; s="-s"; c="" ;;
            TRANS) p=0; s=""; c="" ;;
            SREV_TRANS) p=50; s="-s"; c="" ;;
            SDCJ) p=100; s="-s"; c="-c" ;;
        esac
        mkdir -p $out_folder/N${n}/${mod}
        for l in $((n / 2)) $n
        do
            for op in $((n / 4)) $((n / 2))
            do
                for indels in 0 $((n / 4))
                do
                    d=$((indels / 2)) # DB-exe applies d deletions and d insertions
                    file=$out_folder/N${n}/${mod}/N${n}_${mod}${op}_I$((2 * d))_L${l}_F0.txt
                    [ -f $file ] && continue
                    ./DB-exe $s -d $d -k $k -n $n -p $p -r $op $c -o $file RandDB -l $l
                done
            done
        done
    done
done

cd - > /dev/null
