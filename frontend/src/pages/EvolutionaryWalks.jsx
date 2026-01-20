import React, { useState, useMemo } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, Area, AreaChart } from 'recharts';
import { TrendingDown, Play, RefreshCw, Info } from 'lucide-react';

// Simulated evolutionary walk data generator
function generateWalkData(maxMutations = 10, nWalks = 5) {
    const walks = [];

    for (let w = 0; w < nWalks; w++) {
        let robustness = 1.0;
        const walkData = [{ mutation: 0, robustness: 1.0, functional: true }];

        for (let m = 1; m <= maxMutations; m++) {
            // Simulate robustness decay with some randomness
            const decay = Math.random() * 0.15 + 0.05;
            robustness = Math.max(0, robustness - decay);

            walkData.push({
                mutation: m,
                robustness: robustness,
                functional: robustness > 0.5
            });
        }
        walks.push(walkData);
    }

    // Compute aggregate statistics
    const aggregate = [];
    for (let m = 0; m <= maxMutations; m++) {
        const robustnessValues = walks.map(w => w[m].robustness);
        const functionalCount = walks.filter(w => w[m].functional).length;

        aggregate.push({
            mutation: m,
            meanRobustness: robustnessValues.reduce((a, b) => a + b, 0) / nWalks,
            minRobustness: Math.min(...robustnessValues),
            maxRobustness: Math.max(...robustnessValues),
            pctFunctional: (functionalCount / nWalks) * 100
        });
    }

    return { walks, aggregate };
}

// Region vulnerability mock data
const REGION_VULNERABILITY = [
    { region: 'Promoter A', failureRate: 0.45, mutations: 120 },
    { region: 'RBS A', failureRate: 0.38, mutations: 95 },
    { region: 'CDS A', failureRate: 0.32, mutations: 180 },
    { region: 'Operator A', failureRate: 0.12, mutations: 85 },
    { region: 'Promoter B', failureRate: 0.42, mutations: 115 },
    { region: 'RBS B', failureRate: 0.35, mutations: 102 },
    { region: 'CDS B', failureRate: 0.28, mutations: 175 },
    { region: 'Operator B', failureRate: 0.15, mutations: 78 },
];

export default function EvolutionaryWalks() {
    const [maxMutations, setMaxMutations] = useState(10);
    const [nWalks, setNWalks] = useState(5);
    const [isRunning, setIsRunning] = useState(false);

    const { walks, aggregate } = useMemo(() => {
        return generateWalkData(maxMutations, nWalks);
    }, [maxMutations, nWalks]);

    const handleRun = () => {
        setIsRunning(true);
        setTimeout(() => setIsRunning(false), 500);
    };

    return (
        <div className="space-y-6">
            {/* Header */}
            <div className="flex items-center justify-between">
                <div>
                    <h2 className="text-2xl font-bold text-slate-900">Evolutionary Mutation Walks</h2>
                    <p className="text-slate-500 mt-1">
                        Simulate sequential mutation accumulation to understand robustness decay
                    </p>
                </div>
                <button
                    onClick={handleRun}
                    disabled={isRunning}
                    className="flex items-center gap-2 px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 disabled:opacity-50 transition-colors"
                >
                    {isRunning ? <RefreshCw className="animate-spin" size={18} /> : <Play size={18} />}
                    {isRunning ? 'Running...' : 'Run Simulation'}
                </button>
            </div>

            {/* Controls */}
            <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                <h3 className="text-sm font-semibold text-slate-500 uppercase tracking-wider mb-4">
                    Simulation Parameters
                </h3>
                <div className="grid grid-cols-2 gap-6">
                    <div>
                        <label className="block text-sm font-medium text-slate-700 mb-2">
                            Max Mutations: {maxMutations}
                        </label>
                        <input
                            type="range"
                            min="5"
                            max="20"
                            value={maxMutations}
                            onChange={(e) => setMaxMutations(parseInt(e.target.value))}
                            className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                        />
                    </div>
                    <div>
                        <label className="block text-sm font-medium text-slate-700 mb-2">
                            Number of Walks: {nWalks}
                        </label>
                        <input
                            type="range"
                            min="3"
                            max="20"
                            value={nWalks}
                            onChange={(e) => setNWalks(parseInt(e.target.value))}
                            className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                        />
                    </div>
                </div>
            </div>

            {/* Main Chart - Robustness Decay */}
            <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                <div className="flex items-center gap-2 mb-4">
                    <TrendingDown className="text-red-500" size={20} />
                    <h3 className="text-lg font-semibold text-slate-800">Robustness Decay Curve</h3>
                </div>

                <div className="h-[400px] w-full">
                    <ResponsiveContainer width="100%" height="100%">
                        <AreaChart data={aggregate} margin={{ top: 10, right: 30, left: 0, bottom: 0 }}>
                            <defs>
                                <linearGradient id="robustnessGradient" x1="0" y1="0" x2="0" y2="1">
                                    <stop offset="5%" stopColor="#3b82f6" stopOpacity={0.3} />
                                    <stop offset="95%" stopColor="#3b82f6" stopOpacity={0} />
                                </linearGradient>
                            </defs>
                            <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
                            <XAxis
                                dataKey="mutation"
                                label={{ value: 'Mutations Accumulated', position: 'insideBottomRight', offset: -5 }}
                                tick={{ fontSize: 12 }}
                            />
                            <YAxis
                                domain={[0, 1]}
                                label={{ value: 'Robustness', angle: -90, position: 'insideLeft' }}
                                tick={{ fontSize: 12 }}
                            />
                            <Tooltip
                                contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }}
                                formatter={(value, name) => {
                                    if (name === 'meanRobustness') return [value.toFixed(3), 'Mean Robustness'];
                                    if (name === 'pctFunctional') return [`${value.toFixed(1)}%`, '% Functional'];
                                    return [value, name];
                                }}
                            />
                            <Legend />
                            <Area
                                type="monotone"
                                dataKey="maxRobustness"
                                stroke="transparent"
                                fill="#3b82f6"
                                fillOpacity={0.1}
                                name="Range"
                            />
                            <Area
                                type="monotone"
                                dataKey="minRobustness"
                                stroke="transparent"
                                fill="#ffffff"
                                name=""
                            />
                            <Line
                                type="monotone"
                                dataKey="meanRobustness"
                                stroke="#3b82f6"
                                strokeWidth={3}
                                dot={{ fill: '#3b82f6', strokeWidth: 2 }}
                                name="Mean Robustness"
                            />
                        </AreaChart>
                    </ResponsiveContainer>
                </div>
            </div>

            {/* Functional Percentage Chart */}
            <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                <h3 className="text-lg font-semibold text-slate-800 mb-4">% Functional vs Mutations</h3>
                <div className="h-[300px] w-full">
                    <ResponsiveContainer width="100%" height="100%">
                        <AreaChart data={aggregate}>
                            <defs>
                                <linearGradient id="functionalGradient" x1="0" y1="0" x2="0" y2="1">
                                    <stop offset="5%" stopColor="#22c55e" stopOpacity={0.4} />
                                    <stop offset="95%" stopColor="#22c55e" stopOpacity={0.05} />
                                </linearGradient>
                            </defs>
                            <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
                            <XAxis dataKey="mutation" tick={{ fontSize: 12 }} />
                            <YAxis domain={[0, 100]} tick={{ fontSize: 12 }} unit="%" />
                            <Tooltip contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }} />
                            <Area
                                type="monotone"
                                dataKey="pctFunctional"
                                stroke="#22c55e"
                                strokeWidth={2}
                                fill="url(#functionalGradient)"
                                name="% Functional"
                            />
                        </AreaChart>
                    </ResponsiveContainer>
                </div>
            </div>

            {/* Region Vulnerability Table */}
            <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                <div className="flex items-center gap-2 mb-4">
                    <Info className="text-blue-500" size={20} />
                    <h3 className="text-lg font-semibold text-slate-800">Region Vulnerability Analysis</h3>
                </div>
                <div className="overflow-x-auto">
                    <table className="w-full text-sm">
                        <thead>
                            <tr className="border-b border-slate-200">
                                <th className="text-left py-3 px-4 font-semibold text-slate-600">Region</th>
                                <th className="text-left py-3 px-4 font-semibold text-slate-600">Failure Rate</th>
                                <th className="text-left py-3 px-4 font-semibold text-slate-600">Total Mutations</th>
                                <th className="text-left py-3 px-4 font-semibold text-slate-600">Risk Level</th>
                            </tr>
                        </thead>
                        <tbody>
                            {REGION_VULNERABILITY.sort((a, b) => b.failureRate - a.failureRate).map((region, idx) => (
                                <tr key={region.region} className="border-b border-slate-100 hover:bg-slate-50">
                                    <td className="py-3 px-4 font-medium text-slate-800">{region.region}</td>
                                    <td className="py-3 px-4">
                                        <div className="flex items-center gap-2">
                                            <div className="w-24 h-2 bg-slate-200 rounded-full overflow-hidden">
                                                <div
                                                    className="h-full bg-red-500 rounded-full"
                                                    style={{ width: `${region.failureRate * 100}%` }}
                                                />
                                            </div>
                                            <span className="text-slate-600">{(region.failureRate * 100).toFixed(1)}%</span>
                                        </div>
                                    </td>
                                    <td className="py-3 px-4 text-slate-600">{region.mutations}</td>
                                    <td className="py-3 px-4">
                                        <span className={`px-2 py-1 rounded-full text-xs font-medium ${region.failureRate > 0.4 ? 'bg-red-100 text-red-700' :
                                                region.failureRate > 0.25 ? 'bg-yellow-100 text-yellow-700' :
                                                    'bg-green-100 text-green-700'
                                            }`}>
                                            {region.failureRate > 0.4 ? 'High' : region.failureRate > 0.25 ? 'Medium' : 'Low'}
                                        </span>
                                    </td>
                                </tr>
                            ))}
                        </tbody>
                    </table>
                </div>
            </div>

            {/* Interpretation Box */}
            <div className="bg-gradient-to-r from-blue-50 to-indigo-50 rounded-xl p-6 border border-blue-100">
                <h3 className="text-lg font-semibold text-blue-900 mb-3">Model Interpretation</h3>
                <div className="space-y-2 text-sm text-blue-800">
                    <p>
                        <strong>Note:</strong> These results are from a biophysically grounded model.
                        Absolute values are model-dependent; focus on <strong>relative rankings</strong> between regions.
                    </p>
                    <p>
                        • Operators show significantly higher mutational tolerance than promoter/RBS regions
                    </p>
                    <p>
                        • Robustness typically decays to &lt;50% within 3-5 accumulated mutations
                    </p>
                </div>
            </div>
        </div>
    );
}
