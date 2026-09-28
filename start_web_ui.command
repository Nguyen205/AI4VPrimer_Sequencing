#!/bin/bash
# Move to the script's directory
cd "$(dirname "$0")"

echo "=========================================================="
echo "  Starting AI4VPrimer Amplicon-Sanger Suite Web Interface "
echo "=========================================================="
echo ""
echo "Open your browser at: http://127.0.0.1:5001"
echo "Press Ctrl+C in this terminal to stop the server."
echo ""

# Open browser after 1 second
(sleep 1 && open "http://127.0.0.1:5001") &

# Start python app
python3 app.py
