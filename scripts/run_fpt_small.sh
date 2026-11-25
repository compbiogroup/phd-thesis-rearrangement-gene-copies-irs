# Experiment with the FPT on the small data set of gen_db_fpt.sh (db_fpt/): all
# the instances of N25 and N50, with F = 0 (equality relation, for which the FPT
# is exact), a time limit of TIME_LIMIT seconds per instance (default 600), and
# Greedy, SOARComb and Approx on the same instances. The data set (INPUT_FOLDER)
# is generated first if it does not exist. The results go to OUT_FOLDER (default
# TESE_OUT/ANS_TESE_FPT_SMALL, see run_fpt.sh), and the table is made with
#   python3 table_part_tese.py ../TESE_OUT/ANS_TESE_FPT_SMALL
#
# With the default parameters (64 files with 100 instances each), the estimate
# from a calibration with 60 s per instance is about 30 hours on 16 threads (most
# of it on the instances of N50 that reach the time limit). With N_INST=50 it is
# about half of that.

scripts=$(cd "$(dirname "$0")" && pwd)
input_folder=${INPUT_FOLDER:-$scripts/../db_fpt}
sizes=${SIZES:-"25 50"}

# the output folder of gen_db_fpt.sh is the input folder of the experiment
OUT_FOLDER=$input_folder SIZES=$sizes bash $scripts/gen_db_fpt.sh || exit 1

INPUT_FOLDER=$input_folder \
OUT_FOLDER=${OUT_FOLDER:-$scripts/../TESE_OUT/ANS_TESE_FPT_SMALL} \
SIZES=$sizes \
N_INST=${N_INST:-100} \
TIME_LIMIT=${TIME_LIMIT:-600} \
OP_VALS= I_VALS= L_VALS= F_VALS= \
    bash $scripts/run_fpt.sh
