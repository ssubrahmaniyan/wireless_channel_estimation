#!/bin/bash

# Define output files to store results
OUTPUT_FILE="results.txt"
PILOT_RATIOS="pilot_ratios.txt"
MSES="mses.txt"
BERS="bers.txt"

# Clear the files before appending new results
> $OUTPUT_FILE
> $PILOT_RATIOS
> $MSES
> $BERS

# Run the Python script for quanta values from 1 to 10
for quanta in {1..10}; do
    echo "Running simulation for quanta=$quanta..."
    
    # Update the quanta value dynamically and run the Python script
    python3 -c "
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.api import VAR
import json

quanta = $quanta
$(cat datadriverVAR.py)
    " > $OUTPUT_FILE
    
    # Extract results
    PILOT_RATIO=$(grep "Pilot ratio" $OUTPUT_FILE | awk '{print $NF}')
    AVG_MSE=$(grep "Average MSE" $OUTPUT_FILE | awk '{print $NF}')
    AVG_BER=$(grep "Average BERS" $OUTPUT_FILE | awk '{print $NF}')
    
    # Append to respective files
    echo $PILOT_RATIO >> $PILOT_RATIOS
    echo $AVG_MSE >> $MSES
    echo $AVG_BER >> $BERS
done

# Plot the results
python3 <<EOF
import matplotlib.pyplot as plt

# Read data from files
quanta = range(1, 11)
pilot_ratios = [float(line.strip()) for line in open("$PILOT_RATIOS")]
mses = [float(line.strip()) for line in open("$MSES")]
bers = [float(line.strip()) for line in open("$BERS")]

# Create plots
plt.figure(figsize=(12, 8))

# Plot Pilot Ratio
plt.subplot(3, 1, 1)
plt.plot(quanta, pilot_ratios, marker='o', label="Pilot Ratio")
plt.xlabel("Quanta")
plt.ylabel("Pilot Ratio")
plt.title("Pilot Ratio vs Quanta")
plt.grid(True)
plt.legend()

# Plot Average MSE
plt.subplot(3, 1, 2)
plt.plot(quanta, mses, marker='o', label="Average MSE", color='orange')
plt.xlabel("Quanta")
plt.ylabel("MSE")
plt.title("MSE vs Quanta")
plt.grid(True)
plt.legend()

# Plot BER
plt.subplot(3, 1, 3)
plt.plot(quanta, bers, marker='o', label="Average BER", color='green')
plt.xlabel("Quanta")
plt.ylabel("BER")
plt.title("BER vs Quanta")
plt.grid(True)
plt.legend()

# Adjust layout and save plot
plt.tight_layout()
plt.savefig("results_plot.png")
plt.show()
EOF

echo "Simulations completed. Results plotted and saved as 'results_plot.png'."
