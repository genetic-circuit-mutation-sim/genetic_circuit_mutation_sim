$ErrorActionPreference = "Stop"

Write-Host "Starting New Mutation Scenarios with Varied Toggles..."

# 1. Strong Toggle (Robust)
# High Hill coefficient (cooperativity) and low Kd (strong binding)
Write-Host "`n[1/3] Running Strong Toggle Scenario (n=4, Kd=0.1)..."
python main.py --n 1000 --base-n 4.0 --base-Kd 0.1 --output-dir results/strong_toggle --quiet
Write-Host "Completed Strong Toggle."

# 2. Weak Toggle (Fragile)
# Low Hill coefficient (weak cooperativity) and high Kd (weak binding)
Write-Host "`n[2/3] Running Weak Toggle Scenario (n=1.5, Kd=5.0)..."
python main.py --n 1000 --base-n 1.5 --base-Kd 5.0 --output-dir results/weak_toggle --quiet
Write-Host "Completed Weak Toggle."

# 3. Extreme Mutation Load
# 5 Mutations per variant
Write-Host "`n[3/3] Running Extreme Mutation Load (5 muts/variant)..."
python main.py --n 1000 --mutations-per-variant 5 --output-dir results/extreme_load --quiet
Write-Host "Completed Extreme Load."

Write-Host "`nAll new scenarios completed successfully."
