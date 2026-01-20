import React, { useState, useEffect, useMemo } from 'react';
import { simulateDeterministic } from '../lib/simulation';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import { Play, RotateCcw, Sliders } from 'lucide-react';

const PARAM_CONFIG = [
    { name: 'tx_A', label: 'Transcription A', min: 0.1, max: 5.0, step: 0.1 },
    { name: 'tx_B', label: 'Transcription B', min: 0.1, max: 5.0, step: 0.1 },
    { name: 'tl_A', label: 'Translation A', min: 0.1, max: 5.0, step: 0.1 },
    { name: 'tl_B', label: 'Translation B', min: 0.1, max: 5.0, step: 0.1 },
    { name: 'Kd', label: 'Dissociation (Kd)', min: 0.1, max: 10.0, step: 0.1 },
    { name: 'n', label: 'Hill Coeff (n)', min: 1.0, max: 4.0, step: 0.1 },
];

const PRESETS = {
    wildtype: { tx_A: 1.0, tx_B: 1.0, tl_A: 1.0, tl_B: 1.0, Kd: 1.0, n: 2.0 },
    leaky: { tx_A: 1.0, tx_B: 1.0, tl_A: 1.0, tl_B: 1.0, Kd: 5.0, n: 1.0 }, // High Kd, Low n -> Leaky
    broken_A: { tx_A: 0.1, tx_B: 1.0, tl_A: 1.0, tl_B: 1.0, Kd: 1.0, n: 2.0 },
    strong_repression: { tx_A: 2.0, tx_B: 2.0, tl_A: 1.0, tl_B: 1.0, Kd: 0.5, n: 4.0 },
};

export default function ResultsSimulator() {
    const [params, setParams] = useState(PRESETS.wildtype);
    const [initialCondition, setInitialCondition] = useState('highA'); // highA or highB

    // Run simulation whenever params or initial conditions change
    const simResult = useMemo(() => {
        const init = initialCondition === 'highA'
            ? { A: 10, B: 0.1, mRNA_A: 1.0, mRNA_B: 0.0 }
            : { A: 0.1, B: 10, mRNA_A: 0.0, mRNA_B: 1.0 };

        return simulateDeterministic(params, init);
    }, [params, initialCondition]);

    // Format data for Recharts
    const chartData = useMemo(() => {
        return simResult.time.map((t, i) => ({
            time: t.toFixed(1),
            ProteinA: simResult.A[i],
            ProteinB: simResult.B[i]
        })).filter((_, i) => i % 5 === 0); // Downsample for performance
    }, [simResult]);

    const handleParamChange = (name, value) => {
        setParams(prev => ({ ...prev, [name]: parseFloat(value) }));
    };

    const applyPreset = (presetName) => {
        setParams(PRESETS[presetName]);
    };

    return (
        <div className="flex gap-6 h-[calc(100vh-8rem)]">
            {/* Controls Panel */}
            <div className="w-80 bg-white p-6 rounded-xl shadow-sm border border-slate-100 overflow-y-auto">
                <div className="flex items-center gap-2 mb-6">
                    <Sliders className="text-blue-500" size={20} />
                    <h2 className="font-bold text-slate-800">Parameters</h2>
                </div>

                {/* Presets */}
                <div className="mb-8">
                    <h3 className="text-xs font-semibold text-slate-400 uppercase tracking-wider mb-3">Presets</h3>
                    <div className="grid grid-cols-2 gap-2">
                        {Object.keys(PRESETS).map(key => (
                            <button
                                key={key}
                                onClick={() => applyPreset(key)}
                                className="px-3 py-2 text-xs font-medium bg-slate-50 hover:bg-slate-100 text-slate-600 rounded border border-slate-200 capitalize transition-colors"
                            >
                                {key.replace('_', ' ')}
                            </button>
                        ))}
                    </div>
                </div>

                {/* Initial Conditions Toggle */}
                <div className="mb-8 bg-slate-50 p-3 rounded-lg border border-slate-200">
                    <h3 className="text-xs font-semibold text-slate-500 mb-2">Initial State</h3>
                    <div className="flex rounded-md shadow-sm" role="group">
                        <button
                            type="button"
                            onClick={() => setInitialCondition('highA')}
                            className={`px-4 py-2 text-xs font-medium rounded-l-lg border ${initialCondition === 'highA'
                                ? 'bg-blue-600 text-white border-blue-600'
                                : 'bg-white text-slate-700 border-slate-200 hover:bg-slate-50'
                                }`}
                        >
                            High A
                        </button>
                        <button
                            type="button"
                            onClick={() => setInitialCondition('highB')}
                            className={`px-4 py-2 text-xs font-medium rounded-r-lg border-t border-b border-r ${initialCondition === 'highB'
                                ? 'bg-emerald-600 text-white border-emerald-600'
                                : 'bg-white text-slate-700 border-slate-200 hover:bg-slate-50'
                                }`}
                        >
                            High B
                        </button>
                    </div>
                </div>

                {/* Advanced Simulation Controls */}
                <div className="mb-8 p-4 bg-blue-50/50 rounded-xl border border-blue-100">
                    <h3 className="text-xs font-bold text-blue-800 uppercase tracking-wider mb-4">Physics Settings</h3>

                    <div className="flex items-center justify-between mb-4">
                        <span className="text-xs font-medium text-slate-700">Stochastic (SSA)</span>
                        <label className="relative inline-flex items-center cursor-pointer">
                            <input type="checkbox" className="sr-only peer" />
                            <div className="w-11 h-6 bg-slate-200 peer-focus:outline-none rounded-full peer peer-checked:after:translate-x-full peer-checked:after:border-white after:content-[''] after:absolute after:top-[2px] after:left-[2px] after:bg-white after:border-gray-300 after:border after:rounded-full after:h-5 after:w-5 after:transition-all peer-checked:bg-blue-600"></div>
                        </label>
                    </div>

                    <div>
                        <div className="flex justify-between mb-1">
                            <label className="text-xs font-medium text-slate-700">Ensemble Runs</label>
                            <span className="text-xs font-mono text-blue-600 font-bold">20</span>
                        </div>
                        <input
                            type="range"
                            min="1"
                            max="100"
                            step="1"
                            defaultValue="20"
                            className="w-full h-1 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                        />
                        <p className="text-[10px] text-slate-400 mt-1">Number of traces for probabilistic mean</p>
                    </div>
                </div>

                {/* Sliders */}
                <div className="space-y-6">
                    {PARAM_CONFIG.map(config => (
                        <div key={config.name}>
                            <div className="flex justify-between mb-1">
                                <label className="text-xs font-medium text-slate-600">{config.label}</label>
                                <span className="text-xs font-mono text-slate-500">{params[config.name].toFixed(1)}</span>
                            </div>
                            <input
                                type="range"
                                min={config.min}
                                max={config.max}
                                step={config.step}
                                value={params[config.name]}
                                onChange={(e) => handleParamChange(config.name, e.target.value)}
                                className="w-full h-1 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                            />
                        </div>
                    ))}
                </div>
            </div>

            {/* Visualization Panel */}
            <div className="flex-1 flex flex-col gap-6">
                <div className="flex-1 bg-white p-6 rounded-xl shadow-sm border border-slate-100">
                    <div className="flex justify-between items-center mb-6">
                        <h3 className="font-bold text-slate-800">Simulation Output</h3>
                        <div className="flex items-center gap-4">
                            <div className="flex items-center gap-2">
                                <span className="w-3 h-3 rounded-full bg-blue-500"></span>
                                <span className="text-sm text-slate-600 font-medium">Protein A</span>
                                <span className="text-xs text-slate-400 font-mono">
                                    {simResult.finalState.A.toFixed(2)}
                                </span>
                            </div>
                            <div className="flex items-center gap-2">
                                <span className="w-3 h-3 rounded-full bg-emerald-500"></span>
                                <span className="text-sm text-slate-600 font-medium">Protein B</span>
                                <span className="text-xs text-slate-400 font-mono">
                                    {simResult.finalState.B.toFixed(2)}
                                </span>
                            </div>
                        </div>
                    </div>

                    <div className="h-[calc(100%-4rem)] w-full">
                        <ResponsiveContainer width="100%" height="100%">
                            <LineChart data={chartData} margin={{ top: 5, right: 20, bottom: 5, left: 0 }}>
                                <CartesianGrid stroke="#f1f5f9" strokeDasharray="3 3" />
                                <XAxis
                                    dataKey="time"
                                    label={{ value: 'Time', position: 'insideBottomRight', offset: -5 }}
                                    tick={{ fontSize: 12, fill: '#94a3b8' }}
                                />
                                <YAxis
                                    label={{ value: 'Concentration', angle: -90, position: 'insideLeft' }}
                                    tick={{ fontSize: 12, fill: '#94a3b8' }}
                                />
                                <Tooltip
                                    contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }}
                                />
                                <Legend />
                                <Line
                                    type="monotone"
                                    dataKey="ProteinA"
                                    stroke="#3b82f6"
                                    strokeWidth={2}
                                    dot={false}
                                    activeDot={{ r: 6 }}
                                />
                                <Line
                                    type="monotone"
                                    dataKey="ProteinB"
                                    stroke="#10b981"
                                    strokeWidth={2}
                                    dot={false}
                                    activeDot={{ r: 6 }}
                                />
                            </LineChart>
                        </ResponsiveContainer>
                    </div>
                </div>
            </div>
        </div>
    );
}
