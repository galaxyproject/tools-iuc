#!/usr/bin/bash

# turn off visualization
sed -i 's/"plot_image": true,/"plot_image": false,/g' config.json
sed -i 's/"plot_transcripts": true,/"plot_transcripts": false,/g' config.json
sed -i 's/"plot_cell_seg": true,/"plot_cell_seg": false,/g' config.json
sed -i 's/"plot_nuclear_seg": true,/"plot_nuclear_seg": false,/g' config.json
sed -i 's/"plot_allgenes": true,/"plot_allgenes": false,/g' config.json

touch './bellavista.log' &&
/opt/bellavista/bellavista/bin/bellavista "./config.json" 2>&1 | tee './bellavista.log' &&
TOOL_PID=$!

while sleep 1; do
    if grep -q "Bella Vista input files created!" './bellavista.log'; then
        echo "Bella Vista input files created! Stopping the tool..."
        kill -INT $TOOL_PID
        break
    fi
done &&

wait $TOOL_PID &&
echo "Bella Vista stopped."