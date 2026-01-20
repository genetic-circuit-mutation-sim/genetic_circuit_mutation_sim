/**
 * Genetic Circuit Mutation Simulator - Deterministic Engine (JS Port)
 * 
 * Based on the Python implementation:
 * d(mRNA_A)/dt = tx_A * (Kd^n / (Kd^n + B^n)) + leakiness - deg_mRNA * mRNA_A
 * d(mRNA_B)/dt = tx_B * (Kd^n / (Kd^n + A^n)) + leakiness - deg_mRNA * mRNA_B
 * d(A)/dt = tl_A * mRNA_A - deg_protein * A
 * d(B)/dt = tl_B * mRNA_B - deg_protein * B
 */

const DEFAULT_PARAMS = {
    tx_A: 1.0,
    tx_B: 1.0,
    tl_A: 1.0,
    tl_B: 1.0,
    Kd: 1.0,
    n: 2.0,
    deg_mRNA: 0.1,
    deg_protein: 0.01,
    leakiness: 0.01
};

export const SIMULATION_CONSTANTS = {
    DEFAULT_SIM_TIME: 1000,
    DEFAULT_NUM_POINTS: 1000,
    STEADY_STATE_WINDOW: 50,
    STEADY_STATE_TOLERANCE: 0.01
};

/**
 * Run a deterministic simulation of the toggle switch.
 * @param {Object} params - Circuit parameters
 * @param {Object} initialConditions - Initial conditions {A, B, mRNA_A, mRNA_B}
 * @param {number} simTime - Total simulation time
 * @param {number} numPoints - Number of time points
 */
export function simulateDeterministic(
    params = {},
    initialConditions = { A: 0, B: 0, mRNA_A: 0, mRNA_B: 0 },
    simTime = SIMULATION_CONSTANTS.DEFAULT_SIM_TIME,
    numPoints = SIMULATION_CONSTANTS.DEFAULT_NUM_POINTS
) {
    const p = { ...DEFAULT_PARAMS, ...params };
    const dt = simTime / numPoints;

    // State: [mRNA_A, mRNA_B, A, B]
    let state = [
        initialConditions.mRNA_A || 0,
        initialConditions.mRNA_B || 0,
        initialConditions.A || 0,
        initialConditions.B || 0
    ];

    const timePoints = [0];
    const history = {
        mRNA_A: [state[0]],
        mRNA_B: [state[1]],
        A: [state[2]],
        B: [state[3]]
    };

    for (let t = 0; t < simTime; t += dt) {
        state = rk4(state, dt, p);

        timePoints.push(t + dt);
        history.mRNA_A.push(state[0]);
        history.mRNA_B.push(state[1]);
        history.A.push(state[2]);
        history.B.push(state[3]);
    }

    return {
        time: timePoints,
        ...history,
        finalState: {
            mRNA_A: state[0],
            mRNA_B: state[1],
            A: state[2],
            B: state[3]
        }
    };
}

/**
 * Run 4th Order Runge-Kutta step
 */
function rk4(state, dt, p) {
    const k1 = derivatives(state, p);
    const k2 = derivatives(
        state.map((v, i) => v + k1[i] * dt * 0.5),
        p
    );
    const k3 = derivatives(
        state.map((v, i) => v + k2[i] * dt * 0.5),
        p
    );
    const k4 = derivatives(
        state.map((v, i) => v + k3[i] * dt),
        p
    );

    return state.map((v, i) => v + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * dt / 6);
}

/**
 * Calculate derivatives for the toggle switch system
 */
function derivatives([mRNA_A, mRNA_B, A, B], p) {
    // Ensure non-negative
    const safe_A = Math.max(0, A);
    const safe_B = Math.max(0, B);
    const safe_mRNA_A = Math.max(0, mRNA_A);
    const safe_mRNA_B = Math.max(0, mRNA_B);

    const Kd_n = Math.pow(p.Kd, p.n);
    const hill_B = Kd_n / (Kd_n + Math.pow(safe_B, p.n) + 1e-10);
    const hill_A = Kd_n / (Kd_n + Math.pow(safe_A, p.n) + 1e-10);

    const d_mRNA_A = p.tx_A * hill_B + p.leakiness - p.deg_mRNA * safe_mRNA_A;
    const d_mRNA_B = p.tx_B * hill_A + p.leakiness - p.deg_mRNA * safe_mRNA_B;
    const d_A = p.tl_A * safe_mRNA_A - p.deg_protein * safe_A;
    const d_B = p.tl_B * safe_mRNA_B - p.deg_protein * safe_B;

    return [d_mRNA_A, d_mRNA_B, d_A, d_B];
}
