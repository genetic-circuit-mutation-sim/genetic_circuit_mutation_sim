# 🧪 BioCircuit Mutation Simulator

[![Version](https://img.shields.io/badge/version-1.0.0--alpha-blue.svg)](https://github.com/your-repo/biocircuit-sim)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Biophysics](https://img.shields.io/badge/Domain-Synthetic%20Biology-orange.svg)]()

A high-fidelity **biophysical simulation framework** for analyzing the mutational robustness of synthetic genetic circuits. It transforms binary threshold-based modeling into **probabilistic landscapes**, allowing researchers to predict long-term evolutionary stability.

---

## 🚀 Key Features

### 🧬 Advanced Biophysical Engine (Python)
- **Probabilistic Classification**: Uses ensemble simulations with Gaussian noise to provide phenotype probability distributions instead of brittle binary classifications.
- **Stochastic Dynamics**: Integrated Gillespie SSA for analyzing boundary-case bimodality and noise-induced state switching.
- **Evolutionary Trajectory Analysis**: Sequential mutation accumulation ("mutation walks") to identify critical failure steps and decay curves.
- **Buffered Parameter Mapping**: Neutral zones [0.7, 1.3] and log-normal scaling ensure realistic responses to minor point mutations.
- **Sensitivity Suite**: Automated threshold sweeps and mutation-rate sensitivity analysis using PRCC and local derivatives.

### 🖥️ Cybernetic Dashboard (React/Vite)
- **Interactive Multi-mutation Walks**: Visualize robustness decay in real-time.
- **Stability Heatmaps**: Map mutants onto 2D parameter stability landscapes.
- **Threshold Tuning**: Interactive sliders to adjust classification sensitivity and view its impact on overall robustness.
- **Validation Suite**: Built-in sanity checks verified against documented biological facts (e.g., Promoter Strength Ordering).

---

## 🛠️ Installation

### Backend (Python)
```bash
# Install core dependencies
pip install -r requirements.txt
```

### Frontend (Web UI)
```bash
cd frontend
npm install
npm run dev
```

---

## 🏃 Quick Start

### 1. Run Probabilistic Analysis (CLI)
Generate 1000 mutant variants with ensemble-based probability:
```bash
python main.py --n 1000 --probabilistic --ensemble-runs 20 --output-dir ./results
```

### 2. Run Evolutionary Walk
Simulate accumulated mutations until the circuit breaks:
```bash
python main.py --multi-mutation --n 100 --output-dir ./evolution_study
```

### 3. Launch the Dashboard
Open [http://localhost:5173](http://localhost:5173) to explore results interactively.

---

## 📊 Analytical Insights

### Mutation-to-Parameter Mapping
| Region | Biological Marker | Numerical Effect |
|--------|-------------------|-------------------|
| **Promoter** | -10/-35 Box PWM | Log-Normal Transcription Multiplier |
| **RBS** | SD-Complementarity | Translation Initiation Rate Shift |
| **CDS** | Synonymous Neutrality | Functional impact (Nonsense = 0.0) |
| **Operator** | Binding Affinity | $K_d$ and Hill Coefficient perturbations |

### Robustness Definition
Robustness is quantified as the **probability of phenotype persistence** across an ensemble of noise-perturbed parameters, calculated using Wilson Score confidence intervals.

---

## 📂 Project Structure
```
.
├── analysis/           # Probabilistic classifier, evolutionary walks, sensitivity
├── mapping/            # Biophysical mapping (Buffered & Log-Normal)
├── simulation/         # Deterministic ODE & Stochastic SSA engines
├── visualization/      # Static plot generation
├── frontend/           # Modern React/Vite dashboard
├── docs/               # Reproducibility appendix & technical specs
└── main.py             # Feature-rich CLI entry point
```

---

## 📝 Reproducibility
Full details on Hill coefficients, dissociation constants, and mutation severity profiles are documented in the [Reproducibility Appendix](docs/reproducibility.md).

---

## 📜 License
MIT License.

## 🤝 Citation
If you use this simulator for your research, please cite the **BioCircuit Mutation Engine**. 
*"Analyzing the mutational landscape of synthetic toggle switches via ensemble simulation."*
