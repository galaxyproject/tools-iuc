#!/usr/bin/bash

# turn off visualization
sed -i 's/"plot_image": true,/"plot_image": false,/g' config.json
sed -i 's/"plot_transcripts": true,/"plot_transcripts": false,/g' config.json
sed -i 's/"plot_cell_seg": true,/"plot_cell_seg": false,/g' config.json
sed -i 's/"plot_nuclear_seg": true,/"plot_nuclear_seg": false,/g' config.json
sed -i 's/"plot_allgenes": true,/"plot_allgenes": false,/g' config.json

# start time
START_TIME=$(date +%s)

echo "Time limit set to ${TIME_LIMIT_SECONDS} seconds"

touch './bellavista.log'
/opt/bellavista/bellavista/bin/bellavista "./config.json" 2>&1 | tee './bellavista.log' &
TOOL_PID=$!

# function to handle complete cleanup
cleanup() {
    # Kill all processes
    pkill -P $TOOL_PID 2>/dev/null || true
    kill -9 $TOOL_PID 2>/dev/null || true

    # Make sure any background processes from this script are terminated
    jobs -p | xargs -r kill -9 2>/dev/null || true

    echo "Bella Vista stopped."
    exit 0
}

# Set trap to ensure cleanup on script exit
trap cleanup EXIT INT TERM

while sleep 1; do
    # Check if the job is finished
    if grep -q "Bella Vista input files created!" './bellavista.log' && \
    [ -d "./BellaVista_output/OMEzarrImages" ] && \
    [ -f "./BellaVista_output/OMEzarrImages/.zgroup" ]; then
        echo "Bella Vista input files created! Stopping the tool..."
        kill -INT $TOOL_PID
        sleep 1
        # If still running, use stronger signal
        if ps -p $TOOL_PID > /dev/null; then
            kill -9 $TOOL_PID
        fi
        break
    elif grep -q "Bella Vista input files created!" './bellavista.log'; then
        echo "Log indicates completion but output directory structure is incomplete. Continuing..."
    fi

    # Check timeout
    CURRENT_TIME=$(date +%s)
    ELAPSED_TIME=$((CURRENT_TIME - START_TIME))

    if [[ "$ELAPSED_TIME" =~ ^[0-9]+$ ]] && [[ "$TIME_LIMIT_SECONDS" =~ ^[0-9]+$ ]]; then
        if [ $ELAPSED_TIME -ge $TIME_LIMIT_SECONDS ]; then
            echo "Time limit of ${TIME_LIMIT_SECONDS} seconds reached. Please contact admins. Stopping the tool..."
            kill -INT $TOOL_PID
            sleep 1
            # If still running, use SIGKILL (force quit)
            if ps -p $TOOL_PID > /dev/null; then
                kill -9 $TOOL_PID
            fi
            break
        fi
    fi
done

# Final cleanup
cleanup

# check this line in test. It should not be printed
echo "Script completed."