#!/bin/bash  
# script to create all minimal distances with support for multiple trajectories

# Print usage instructions
usage() {
    echo "Usage: $0 -traj <path> -pdb <file> -min <val> -max <val> -ndx <file> -sys <name> -traj_mode <mode> -mode <mode> [-ext <ext>] [-stride <n>] [-atom_thr on|off]"
    echo "  -traj: trajectory file OR path to folder with trajectory files (alias: -xtc)"
    echo "         supported formats: xtc, trr, dcd, nc, ncdf, mdcrd, xyz, gro, pdb"
    echo "  -pdb: pdb file"
    echo "  -min: minimum threshold (0-1)"
    echo "  -max: maximum threshold (0-1)"
    echo "  -ndx: indices file"
    echo "  -sys: base name for system, for output files"
    echo "  -traj_mode: traj mode [multi|single] trajectories"
    echo "  -mode: threshold mode [overall|per-trajectory]"
    echo "           overall: contacts averaged across all trajectories"
    echo "           per-trajectory: contacts must meet threshold in each trajectory"
    echo "  -ext: (optional, multi mode) only use files with this extension in the folder,"
    echo "        required if the folder contains trajectories of more than one format"
    echo "  -stride: (optional) use only every n-th frame of each trajectory (default: 1 = all frames)"
    echo "  -atom_thr: (optional) apply the min threshold also to the atom pairs in step 4 [on|off] (default: on)"
    echo "           on:  keep residue pairs with at least one atom pair in contact >= min threshold,"
    echo "                distance = minimum over those atom pairs"
    echo "           off: keep all residue pairs selected in step 2,"
    echo "                distance = minimum over all heavy-atom pairs"
    echo ""
    echo "Examples:"
    echo "  Single trajectory:"
    echo "    $0 -traj traj.xtc -pdb system.pdb -min 0.3 -max 0.9 -ndx indices.ndx -sys system_name -traj_mode single -mode overall"
    echo ""
    echo "  Multiple trajectories (overall mode):"
    echo "    $0 -traj /path/to/traj_folder -pdb system.pdb -min 0.3 -max 0.9 -ndx indices.ndx -sys system_name -traj_mode multi -mode overall"
    echo ""
    echo "  Multiple trajectories (per-trajectory mode):"
    echo "    $0 -traj /path/to/traj_folder -pdb system.pdb -min 0.3 -max 0.9 -ndx indices.ndx -sys system_name -traj_mode multi -mode per-trajectory"
    echo ""
    echo "  Multiple .dcd trajectories in a folder that also contains .xtc files:"
    echo "    $0 -traj /path/to/traj_folder -pdb system.pdb -min 0.3 -max 0.9 -ndx indices.ndx -sys system_name -traj_mode multi -mode overall -ext dcd"
    echo ""
    echo "  Every 10th frame only (faster for long trajectories):"
    echo "    $0 -traj traj.xtc -pdb system.pdb -min 0.3 -max 0.9 -ndx indices.ndx -sys system_name -traj_mode single -mode overall -stride 10"
    echo ""
    echo "  Folder should contain trajectory files like:"
    echo "    /path/to/traj_folder/traj1.xtc"
    echo "    /path/to/traj_folder/traj2.xtc"
    echo "    /path/to/traj_folder/traj3.xtc"
} 

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -traj|-xtc) TRAJ="$2"; shift ;;
        -pdb) PDB="$2"; shift ;;
        -min) MIN_THR="$2"; shift ;;
        -max) MAX_THR="$2"; shift ;;
        -ndx) INDEX="$2"; shift ;;
        -sys) SYSTEM="$2"; shift ;;
        -traj_mode) TRAJ_MODE="$2"; shift ;;
        -mode) MODE="$2"; shift ;;
        -ext) TRAJ_EXT="$2"; shift ;;
        -stride) STRIDE="$2"; shift ;;
        -atom_thr) ATOM_THR="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
    shift
done

# Check if all required arguments are provided
if [ -z "$TRAJ" ] || [ -z "$PDB" ] || [ -z "$MIN_THR" ] || [ -z "$MAX_THR" ] || [ -z "$INDEX" ] || [ -z "$SYSTEM" ] || [ -z "$TRAJ_MODE" ] || [ -z "$MODE" ]; then
    echo "Error: Missing required arguments."
    usage
    exit 1
fi

# Validate trajectory mode
if [ "$TRAJ_MODE" != "single" ] && [ "$TRAJ_MODE" != "multi" ]; then
    echo "Error: -traj_mode must be 'single' or 'multi'"
    usage
    exit 1
fi

# Validate contact mode
if [ "$MODE" != "overall" ] && [ "$MODE" != "per-trajectory" ]; then
    echo "Error: -mode must be 'overall' or 'per-trajectory'"
    usage
    exit 1
fi  

if [ "$TRAJ_MODE" == "single" ] && [ "$MODE" == "per-trajectory" ]; then
    echo "Error: threshold mode 'per-trajectory' is only allowed with trajectory mode 'multi'"
    usage
    exit 1
fi

# Check if multi-trajectory mode
if [ "$TRAJ_MODE" == "multi" ]; then
    TRAJ_ARG="--trajectory-list"
    echo "Running in multi-trajectory mode"
    echo "Trajectory folder: $TRAJ"
    echo "Mode: $MODE"
    if [ "$MODE" == "overall" ]; then
        echo "threshold is set over the total length of simulated data"
    elif [ "$MODE" == "per-trajectory" ]; then
        echo "threshold is set over each single trajectory"
    fi
else
    TRAJ_ARG="--trajectory"
    echo "Running in single-trajectory mode"
    echo "Trajectory file: $TRAJ"
fi  

ATOM_THR="${ATOM_THR:-on}"
if [ "$ATOM_THR" == "on" ]; then
    ATOM_THR_ARG="--atom-threshold"
elif [ "$ATOM_THR" == "off" ]; then
    ATOM_THR_ARG="--no-atom-threshold"
else
    echo "Error: -atom_thr must be 'on' or 'off'"
    usage
    exit 1
fi

STRIDE="${STRIDE:-1}"
if ! [[ "$STRIDE" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: -stride must be a positive integer"
    usage
    exit 1
fi

if [ ! -f "$PDB" ] || [ ! -f "$INDEX" ]; then
    echo "Error: PDB or INDEX file does not exist"
    usage
    exit 1
fi

# Check trajectory based on mode
EXT_ARG=()
if [ "$TRAJ_MODE" == "multi" ]; then
    if [ ! -d "$TRAJ" ]; then
        echo "Error: In multi mode, -traj must be a directory containing trajectory files"
        usage
        exit 1
    fi
    if [ -n "$TRAJ_EXT" ]; then
        EXT_ARG=(--traj-ext "$TRAJ_EXT")
    fi
else
    if [ ! -f "$TRAJ" ]; then
        echo "Error: In single mode, -traj must be a trajectory file"
    usage
    exit 1
    fi
fi 

# created files  
THR_SUFFIX="${MIN_THR}-${MAX_THR}"  
IS_MINDIST="${SYSTEM}.is_mindist"
ATOMDIST="${SYSTEM}.all_thr${THR_SUFFIX}_selected_atom_distances"

# estimate all mindist
if [[ ! -e "${IS_MINDIST}" ]]; then
    echo ""
    echo "step 1 - estimate_contacts.py: compute contacts"
    time python estimate_contacts.py \
        --top $PDB \
        $TRAJ_ARG $TRAJ "${EXT_ARG[@]}" \
        --index $INDEX \
        --output $IS_MINDIST \
        --mode $MODE \
        --stride $STRIDE \
        --verbose
else
    # the stride is recorded in the header (only if > 1); never reuse a
    # result computed with a different stride
    if [ "$STRIDE" -gt 1 ]; then
        grep -q "stride ${STRIDE})" "${IS_MINDIST}"; SAME_STRIDE=$?
    else
        ! grep -q "stride [0-9]*)" "${IS_MINDIST}"; SAME_STRIDE=$?
    fi
    if [ "$SAME_STRIDE" -ne 0 ]; then
        echo "Error: ${IS_MINDIST} exists but was computed with a different stride."
        echo "       Delete it or use another -sys name."
        exit 1
    fi
    echo "skipping step 1 (${IS_MINDIST} already exists)"
fi

# select formed mindist
echo ""
echo "step 2 - extract_indices.py: select residue pairs where mindist is contact for more than min threshold time and less than max threshold time"
time python extract_indices.py \
    --is-contacts $IS_MINDIST \
    --threshold $MIN_THR \
    --max-threshold $MAX_THR 

# extract all atom pairwise distances of selected residues
echo ""
echo "step 3 - contacts.py: extract all atom pairwise dists between residues in .ndx file"
time python contacts.py \
    $TRAJ_ARG $TRAJ "${EXT_ARG[@]}" \
    -s $PDB \
    -n ${IS_MINDIST}.thr${THR_SUFFIX}.ndx \
    -o $ATOMDIST \
    --stride $STRIDE

# extract minimal distances between all atom pairs forming a contact more often
# than the the given threshold
echo ""
echo "step 4 - extract_contacts.py: extracting final minimal distances"
time python extract_contacts.py \
    --contacts $ATOMDIST \
    --index ${ATOMDIST}.atom_indices \
    --threshold ${MIN_THR} \
    $ATOM_THR_ARG \
    --output ${SYSTEM}.mindist  

echo ""
echo "Final output files:"
echo "  - ${SYSTEM}.mindist (final minimal distances)"
echo "  - ${SYSTEM}.mindist.ndx (indices of selected residue pairs)"
#end
