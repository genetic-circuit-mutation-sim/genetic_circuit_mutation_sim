"""
Simulation Module

Provides simulation engines for the toggle switch genetic circuit.
Supports both deterministic ODE-based simulation (via Tellurium/RoadRunner)
and stochastic Gillespie algorithm simulation.

Public API:
    - simulate_deterministic: ODE-based simulation with steady-state analysis
    - simulate_gillespie: Stochastic simulation using Gillespie SSA
"""

from .deterministic import simulate_deterministic
from .stochastic import simulate_gillespie

__all__ = ['simulate_deterministic', 'simulate_gillespie']
