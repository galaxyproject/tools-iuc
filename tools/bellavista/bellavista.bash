#!/usr/bin/bash

cd "$BELLAVISTA_DIR" &&
touch './bellavista.log' &&
/opt/bellavista/bellavista/bin/bellavista "./config.json" | tee './bellavista.log' &&
TOOL_PID=$!

while sleep 1; do
    if grep -q "OME-Zarr image saved successfully" './bellavista.log'; then
        echo "Bella Vista input files created! Stopping the tool..."
        kill -INT $TOOL_PID
        break
    fi
done &&

wait $TOOL_PID &&
echo "Bella Vista stopped."