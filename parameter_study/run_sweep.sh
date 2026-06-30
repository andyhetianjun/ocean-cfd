#!/bin/bash
# Parameter sweep for waves2Foam monopile study
# Run from: /shared_folder/andyhe/project/waves2foam/
# Usage: bash run_sweep.sh

BASE="/shared_folder/andyhe/project/waves2foam/monopileSolitary"
WORK="/shared_folder/andyhe/project/waves2foam"

source /usr/lib/openfoam/openfoam2412/etc/bashrc
source /shared_folder/andyhe/tools/waves2Foam/bin/bashrc

# ── case definitions ─────────────────────────────────────────────────────
# Format: "caseName  waveType  height  period  notes"
# period only used for regular waves (stokesFirst/stokesSecond)
# for solitary waves period is ignored

CASES=(
    "solitary_H005  solitaryFirst  0.05  -    solitary_H/d=0.125"
    "solitary_H010  solitaryFirst  0.10  -    solitary_H/d=0.250_BASELINE"
    "solitary_H015  solitaryFirst  0.15  -    solitary_H/d=0.375"
    "solitary_H020  solitaryFirst  0.20  -    solitary_H/d=0.500"
    "regular_T2     stokesFirst    0.10  2.0  regular_T=2s"
    "regular_T3     stokesFirst    0.10  3.0  regular_T=3s"
    "regular_T5     stokesFirst    0.10  5.0  regular_T=5s"
    "regular_T8     stokesFirst    0.10  8.0  regular_T=8s"
    "stokes2_T3     stokesSecond   0.10  3.0  stokes2nd_T=3s"
)

# ── function: create one case ─────────────────────────────────────────────
create_case() {
    local NAME=$1
    local WTYPE=$2
    local HEIGHT=$3
    local PERIOD=$4
    local DEST="$WORK/$NAME"

    echo ""
    echo "=================================================="
    echo "Creating case: $NAME"
    echo "  waveType=$WTYPE  H=$HEIGHT  T=$PERIOD"
    echo "=================================================="

    # copy baseline
    rm -rf "$DEST"
    cp -r "$BASE" "$DEST"

    # clean old results but keep mesh and 0/
    cd "$DEST"
    for d in $(ls -d [0-9]* 2>/dev/null); do
        [ "$d" != "0" ] && rm -rf "$d"
    done
    rm -rf postProcessing processor* log.*

    # ── waveProperties ──────────────────────────────────────────────────
    WP="$DEST/constant/waveProperties"

    if [ "$WTYPE" = "solitaryFirst" ]; then
        # solitary: no period needed
        cat > "$WP" << WPEOF
FoamFile { version 2.0; format ascii; class dictionary; object waveProperties; }
seaLevel            0.4;
seaLevelAsReference true;
relaxationNames     ( );
initializationName  inlet;
inletCoeffs
{
    waveType            solitaryFirst;
    height              $HEIGHT;
    depth               0.4;
    direction           ( 1 0 0 );
    x0                  ( 0 0 0 );
    GABCSettings
    {
        preProcessMethod    polynomialDefault;
        defaultRange        "kh000_030";
        shapeFunction       P03Celerity;
        a0                  1.6207; a1 0; a2 -2.7209; a3 1.5429;
        depth               0.4;
    }
}
outletCoeffs
{
    waveType            potentialCurrent;
    U                   ( 0 0 0 );
    Tsoft               2;
    GABCSettings
    {
        preProcessMethod    polynomialDefault;
        defaultRange        "kh000_030";
        shapeFunction       P03Celerity;
        a0                  1.6207; a1 0; a2 -2.7209; a3 1.5429;
        depth               0.4;
    }
}
WPEOF

    else
        # regular waves: need period, compute omega and wavenumber via python
        OMEGA=$(python3 -c "
import math
T=$PERIOD; d=0.4; g=9.81
omega=2*math.pi/T
# iterate dispersion relation: omega^2 = g*k*tanh(k*d)
k=omega**2/g
for _ in range(50): k=omega**2/(g*math.tanh(k*d))
print(f'{omega:.6f}')
")
        WNUM=$(python3 -c "
import math
T=$PERIOD; d=0.4; g=9.81
omega=2*math.pi/T
k=omega**2/g
for _ in range(50): k=omega**2/(g*math.tanh(k*d))
print(f'{k:.6f}')
")
        echo "  omega=$OMEGA  k=$WNUM"

        cat > "$WP" << WPEOF
FoamFile { version 2.0; format ascii; class dictionary; object waveProperties; }
seaLevel            0.4;
seaLevelAsReference true;
relaxationNames     ( );
initializationName  inlet;
inletCoeffs
{
    waveType            $WTYPE;
    height              $HEIGHT;
    depth               0.4;
    period              $PERIOD;
    direction           ( 1 0 0 );
    phi                 0.0;
    GABCSettings
    {
        preProcessMethod    polynomialDefault;
        defaultRange        "kh000_030";
        shapeFunction       P03Celerity;
        a0                  1.6207; a1 0; a2 -2.7209; a3 1.5429;
        depth               0.4;
    }
}
outletCoeffs
{
    waveType            potentialCurrent;
    U                   ( 0 0 0 );
    Tsoft               2;
    GABCSettings
    {
        preProcessMethod    polynomialDefault;
        defaultRange        "kh000_030";
        shapeFunction       P03Celerity;
        a0                  1.6207; a1 0; a2 -2.7209; a3 1.5429;
        depth               0.4;
    }
}
WPEOF
    fi

    echo "  waveProperties written"

    # ── reset 0/ fields for fresh setWaveField ───────────────────────────
    # restore clean initial fields (uniform)
    for fld in U p_rgh alpha.water; do
        # just reset internalField to uniform — setWaveField overwrites anyway
        python3 -c "
import re
f=open('0/$fld').read()
# reset nonuniform internal fields to uniform
f=re.sub(r'internalField\s+nonuniform[^;]+;','internalField   uniform 0;',f,flags=re.DOTALL)
open('0/$fld','w').write(f)
" 2>/dev/null || true
    done
    # U needs vector uniform
    python3 -c "
import re
f=open('0/U').read()
f=re.sub(r'internalField\s+uniform\s+0;','internalField   uniform (0 0 0);',f)
open('0/U','w').write(f)
" 2>/dev/null || true

    echo "  fields reset"
}

# ── function: run one case ────────────────────────────────────────────────
run_case() {
    local NAME=$1
    local DEST="$WORK/$NAME"
    cd "$DEST"

    echo ""
    echo "Running: $NAME"
    echo "  setWaveField ..."
    setWaveField > log.setWaveField 2>&1
    if grep -q "FATAL" log.setWaveField; then
        echo "  ERROR in setWaveField — check log.setWaveField"
        return 1
    fi
    echo "  setWaveField OK"

    echo "  waveFoam ..."
    waveFoam > log.waveFoam 2>&1
    if grep -q "FATAL" log.waveFoam; then
        echo "  ERROR in waveFoam — check log.waveFoam"
        return 1
    fi
    echo "  waveFoam OK — $(tail -3 log.waveFoam | grep ExecutionTime)"
    return 0
}

# ── main loop ─────────────────────────────────────────────────────────────
cd "$WORK"
echo "Starting parameter sweep: ${#CASES[@]} cases"
echo "Baseline (already run): monopileSolitary"
echo ""

FAILED=()
for entry in "${CASES[@]}"; do
    read NAME WTYPE HEIGHT PERIOD NOTES <<< "$entry"
    # skip baseline — already exists
    if [ "$NAME" = "solitary_H010" ]; then
        echo "Skipping $NAME (same as baseline monopileSolitary)"
        continue
    fi
    create_case "$NAME" "$WTYPE" "$HEIGHT" "$PERIOD"
    run_case "$NAME" || FAILED+=("$NAME")
done

echo ""
echo "=================================================="
echo "SWEEP COMPLETE"
echo "Successful cases:"
for entry in "${CASES[@]}"; do
    read NAME _ <<< "$entry"
    [ "$NAME" != "solitary_H010" ] && ! [[ " ${FAILED[@]} " =~ " $NAME " ]] && echo "  ✓ $NAME"
done
if [ ${#FAILED[@]} -gt 0 ]; then
    echo "Failed cases:"
    for n in "${FAILED[@]}"; do echo "  ✗ $n"; done
fi
echo "=================================================="
