$ErrorActionPreference = "Stop"

Write-Host "Starting Batch Mutation Scenarios..."

# 1. High Substitution Scenario
Write-Host "`n[1/4] Running High Substitution Scenario..."
python main.py --n 1000 --substitution-rate 0.9 --insertion-rate 0.05 --deletion-rate 0.05 --output-dir results/substitution_dominant --quiet
Write-Host "Completed High Substitution."

# 2. High Insertion Scenario
Write-Host "`n[2/4] Running High Insertion Scenario..."
python main.py --n 1000 --substitution-rate 0.1 --insertion-rate 0.8 --deletion-rate 0.1 --output-dir results/insertion_dominant --quiet
Write-Host "Completed High Insertion."

# 3. High Deletion Scenario
Write-Host "`n[3/4] Running High Deletion Scenario..."
python main.py --n 1000 --substitution-rate 0.1 --insertion-rate 0.1 --deletion-rate 0.8 --output-dir results/deletion_dominant --quiet
Write-Host "Completed High Deletion."

# 4. High Mutational Load Scenario
Write-Host "`n[4/4] Running High Mutational Load (3 mutations/variant)..."
python main.py --n 1000 --mutations-per-variant 3 --output-dir results/multi_hit --quiet
Write-Host "Completed High Mutational Load."

Write-Host "`nAll scenarios completed successfully."
