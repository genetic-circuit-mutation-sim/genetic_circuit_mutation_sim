import React, { useState, useMemo } from 'react';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, ZAxis } from 'recharts';
import { Layers, Info, Filter } from 'lucide-react';

// Mock parameter space data
const PARAMETERS = [
    { id: 'k_tx_a', name: 'Transcription Rate A (k_tx_a)', range: [0.1, 1.5] },
    { id: 'k_tx_b', name: 'Transcription Rate B (k_tx_b)', range: [0.1, 1.5] },
    { id: 'k_tl_a', name: 'Translation Rate A (k_tl_a)', range: [0.1, 1.5] },
    { id: 'k_tl_b', name: 'Translation Rate B (k_tl_b)', range: [0.1, 1.5] },
    { id: 'kd_a', name: 'Dissociation Const A (kd_a)', range: [1, 1000] },
    { id: 'kd_b', name: 'Dissociation Const B (kd_b)', range: [1, 1000] },
];

function generateSpaceData(xAxis, yAxis) {
    const data = [];
    const xParam = PARAMETERS.find(p => p.id === xAxis);
    const yParam = PARAMETERS.find(p => p.id === yAxis);

    // Generate a grid of points
    for (let x = 0; x < 20; x++) {
        for (let y = 0; y < 20; y++) {
            const xVal = xParam.range[0] + (x / 20) * (xParam.range[1] - xParam.range[0]);
            const yVal = yParam.range[0] + (y / 20) * (yParam.range[1] - yParam.range[0]);

            // Bistability region mock (roughly central-diagonal)
            const dist = Math.abs(xVal - yVal) + Math.random() * 0.1;
            const isBistable = dist < 0.3 && xVal > 0.4 && yVal > 0.4;

            data.push({
                x: xVal,
                y: yVal,
                z: isBistable ? 100 : 20,
                status: isBistable ? 'Bistable' : 'Monostable'
            });
        }
    }
    return data;
}

export default function ParameterSpaceExplorer() {
    const [xAxis, setXAxis] = useState('k_tx_a');
    const [yAxis, setYAxis] = useState('k_tx_b');

    const data = useMemo(() => generateSpaceData(xAxis, yAxis), [xAxis, yAxis]);

    return (
        <div className="space-y-6">
            <div>
                <h2 className="text-2xl font-bold text-slate-900">Parameter Space Explorer</h2>
                <p className="text-slate-500 mt-1">
                    Visualize the stability landscape and mutant distributions
                </p>
            </div>

            {/* Controls */}
            <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                <div className="flex items-center gap-4">
                    <div className="flex-1">
                        <label className="block text-sm font-medium text-slate-700 mb-1">X-Axis Parameter</label>
                        <select
                            value={xAxis}
                            onChange={(e) => setXAxis(e.target.value)}
                            className="w-full bg-slate-50 border border-slate-200 rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
                        >
                            {PARAMETERS.map(p => <option key={p.id} value={p.id}>{p.name}</option>)}
                        </select>
                    </div>
                    <div className="flex-1">
                        <label className="block text-sm font-medium text-slate-700 mb-1">Y-Axis Parameter</label>
                        <select
                            value={yAxis}
                            onChange={(e) => setYAxis(e.target.value)}
                            className="w-full bg-slate-50 border border-slate-200 rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
                        >
                            {PARAMETERS.map(p => <option key={p.id} value={p.id}>{p.name}</option>)}
                        </select>
                    </div>
                    <div className="flex items-end pb-1">
                        <button className="p-2 bg-slate-100 hover:bg-slate-200 rounded-lg transition-colors">
                            <Filter size={18} className="text-slate-600" />
                        </button>
                    </div>
                </div>
            </div>

            {/* Heatmap Area */}
            <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                <div className="flex items-center justify-between mb-6">
                    <h3 className="text-lg font-semibold text-slate-800 flex items-center gap-2">
                        <Layers size={20} className="text-blue-500" />
                        Stability Landscape
                    </h3>
                    <div className="flex items-center gap-4 text-xs font-medium">
                        <div className="flex items-center gap-1">
                            <div className="w-3 h-3 bg-blue-500 rounded-full" />
                            <span>Bistable Region</span>
                        </div>
                        <div className="flex items-center gap-1">
                            <div className="w-3 h-3 bg-slate-200 rounded-full" />
                            <span>Monostable / Failed</span>
                        </div>
                    </div>
                </div>

                <div className="h-[500px] w-full">
                    <ResponsiveContainer width="100%" height="100%">
                        <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
                            <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
                            <XAxis
                                type="number"
                                dataKey="x"
                                name={PARAMETERS.find(p => p.id === xAxis)?.name}
                                domain={PARAMETERS.find(p => p.id === xAxis)?.range}
                                label={{ value: xAxis, position: 'insideBottom', offset: -10 }}
                                tick={{ fontSize: 11 }}
                            />
                            <YAxis
                                type="number"
                                dataKey="y"
                                name={PARAMETERS.find(p => p.id === yAxis)?.name}
                                domain={PARAMETERS.find(p => p.id === yAxis)?.range}
                                label={{ value: yAxis, angle: -90, position: 'insideLeft' }}
                                tick={{ fontSize: 11 }}
                            />
                            <ZAxis type="number" dataKey="z" range={[50, 400]} />
                            <Tooltip
                                cursor={{ strokeDasharray: '3 3' }}
                                contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }}
                            />
                            <Scatter
                                name="Parameter Points"
                                data={data}
                                fill="#3b82f6"
                            >
                                {data.map((entry, index) => (
                                    <Cell
                                        key={`cell-${index}`}
                                        fill={entry.status === 'Bistable' ? '#3b82f6' : '#e2e8f0'}
                                        fillOpacity={entry.status === 'Bistable' ? 0.7 : 0.4}
                                    />
                                ))}
                            </Scatter>
                        </ScatterChart>
                    </ResponsiveContainer>
                </div>
            </div>

            {/* Guidance Box */}
            <div className="bg-slate-50 border border-slate-200 rounded-xl p-6">
                <div className="flex items-start gap-3">
                    <Info className="text-slate-400 mt-1" size={20} />
                    <div>
                        <h4 className="font-semibold text-slate-800">Understanding the Landscape</h4>
                        <p className="text-sm text-slate-600 mt-1">
                            The blue region represents the set of parameter values where the toggle switch circuit exhibits <strong>robust bistability</strong>. Mutations map sequences onto this landscape; if a mutation shifts the parameters into the gray region, the functional phenotype is lost.
                        </p>
                    </div>
                </div>
            </div>
        </div>
    );
}

import { Cell } from 'recharts';
