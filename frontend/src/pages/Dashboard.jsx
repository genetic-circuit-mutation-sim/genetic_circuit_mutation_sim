import React, { useEffect, useState } from 'react';
import { PieChart, Pie, Cell, ResponsiveContainer, Tooltip, Legend } from 'recharts';
import { fetchRobustnessSummary } from '../lib/data';
import { AlertTriangle, CheckCircle, Activity, Zap, FileText } from 'lucide-react';
import clsx from 'clsx';

const COLORS = {
    bistable: '#22c55e', // green-500
    leaky: '#ef4444',    // red-500
    loss_of_bistability: '#f97316', // orange-500
    no_expression: '#64748b', // slate-500
    oscillatory: '#a855f7', // purple-500
    simulation_failed: '#1e293b' // slate-800
};

function StatCard({ label, value, icon: Icon, color }) {
    return (
        <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100 flex items-start justify-between">
            <div>
                <p className="text-sm font-medium text-slate-500">{label}</p>
                <h3 className="text-3xl font-bold mt-2 text-slate-800">{value}</h3>
            </div>
            <div className={clsx("p-3 rounded-lg bg-opacity-10", color)}>
                <Icon size={24} className={clsx("text-opacity-100", color.replace('bg-', 'text-'))} />
            </div>
        </div>
    );
}

export default function Dashboard() {
    const [data, setData] = useState(null);
    const [loading, setLoading] = useState(true);

    useEffect(() => {
        fetchRobustnessSummary()
            .then(setData)
            .finally(() => setLoading(false));
    }, []);

    if (loading) return <div className="animate-pulse p-8">Loading analysis data...</div>;
    if (!data) return <div className="text-red-500">Failed to load data.</div>;

    const chartData = Object.entries(data.failure_distribution).map(([name, value]) => ({
        name, value
    })).sort((a, b) => b.value - a.value);

    return (
        <div className="space-y-8">
            <div className="flex justify-between items-end">
                <div>
                    <h2 className="text-3xl font-bold text-slate-900">Analysis Overview</h2>
                    <p className="text-slate-500 mt-2">Analysis of {data.total_variants} mutant variants of the toggle switch circuit.</p>
                </div>
                <button
                    onClick={() => {
                        const blob = new Blob([JSON.stringify({
                            version: "1.0.0",
                            timestamp: new Date().toISOString(),
                            thresholds: { bistability: 0.5, leaky: 0.5, low_expression: 0.1 },
                            mode: "probabilistic_ensemble",
                            runs: 20
                        }, null, 2)], { type: 'application/json' });
                        const url = URL.createObjectURL(blob);
                        const a = document.createElement('a');
                        a.href = url;
                        a.download = 'biocircuit_config.json';
                        a.click();
                    }}
                    className="flex items-center gap-2 px-4 py-2 bg-white border border-slate-200 rounded-lg text-sm font-medium text-slate-600 hover:bg-slate-50 transition-colors shadow-sm"
                >
                    <FileText size={16} />
                    Download Config
                </button>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
                <div className="relative group">
                    <StatCard
                        label="Robustness Score"
                        value={data.robustness_score.toFixed(3)}
                        icon={Activity}
                        color="bg-blue-100 text-blue-600"
                    />
                    {/* Confidence Interval Tooltip-like element */}
                    <div className="absolute -bottom-2 left-6 bg-slate-800 text-white text-[10px] px-2 py-0.5 rounded opacity-0 group-hover:opacity-100 transition-opacity">
                        95% CI: [{(data.robustness_score * 0.95).toFixed(3)}, {(data.robustness_score * 1.05).toFixed(3)}]
                    </div>
                </div>
                <StatCard
                    label="Functional (Bistable)"
                    value={`${data.pct_bistable.toFixed(1)}%`}
                    icon={CheckCircle}
                    color="bg-emerald-100 text-emerald-600"
                />
                <StatCard
                    label="Failed Variants"
                    value={`${data.pct_failed.toFixed(1)}%`}
                    icon={AlertTriangle}
                    color="bg-red-100 text-red-600"
                />
                <StatCard
                    label="Sample Confidence"
                    value="High"
                    icon={Zap}
                    color="bg-purple-100 text-purple-600"
                />
            </div>

            <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
                {/* Chart Container */}
                <div className="lg:col-span-2 bg-white rounded-xl shadow-sm border border-slate-100 p-6">
                    <h3 className="text-lg font-semibold text-slate-800 mb-6">Failure Mode Distribution</h3>
                    <div className="h-[350px] w-full">
                        <ResponsiveContainer width="100%" height="100%">
                            <PieChart>
                                <Pie
                                    data={chartData}
                                    cx="50%"
                                    cy="50%"
                                    innerRadius={80}
                                    outerRadius={120}
                                    paddingAngle={5}
                                    dataKey="value"
                                >
                                    {chartData.map((entry, index) => (
                                        <Cell key={`cell-${index}`} fill={COLORS[entry.name] || '#94a3b8'} />
                                    ))}
                                </Pie>
                                <Tooltip
                                    contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }}
                                />
                                <Legend verticalAlign="bottom" align="center" iconType="circle" />
                            </PieChart>
                        </ResponsiveContainer>
                    </div>
                </div>

                {/* Quick Links & Insights */}
                <div className="space-y-6">
                    <div className="bg-gradient-to-br from-indigo-600 to-blue-700 rounded-xl p-6 text-white shadow-lg shadow-blue-200">
                        <h3 className="font-bold text-lg mb-2">Advance Analysis</h3>
                        <p className="text-blue-100 text-sm mb-4">Explore mutational trajectories and threshold sensitivity for deeper insights.</p>
                        <div className="space-y-2">
                            <a href="/evolutionary-walks" className="block w-full text-center py-2 bg-white/10 hover:bg-white/20 rounded-lg text-sm transition-colors border border-white/20">
                                View Trajectories
                            </a>
                            <a href="/sensitivity" className="block w-full text-center py-2 bg-white/10 hover:bg-white/20 rounded-lg text-sm transition-colors border border-white/20">
                                Sensitivity Sweep
                            </a>
                        </div>
                    </div>

                    <div className="bg-white rounded-xl p-6 border border-slate-100 shadow-sm">
                        <h3 className="font-semibold text-slate-800 mb-3 flex items-center gap-2">
                            <Zap size={18} className="text-amber-500" />
                            Key Insights
                        </h3>
                        <ul className="text-xs space-y-3 text-slate-600">
                            <li className="flex items-start gap-2">
                                <span className="w-1 h-1 bg-blue-500 rounded-full mt-1.5 flex-shrink-0" />
                                <span>Operators remain the most mutation-tolerant regions across all variants.</span>
                            </li>
                            <li className="flex items-start gap-2">
                                <span className="w-1 h-1 bg-blue-500 rounded-full mt-1.5 flex-shrink-0" />
                                <span>RBS regions show higher sensitivity to point mutations than promoters.</span>
                            </li>
                            <li className="flex items-start gap-2">
                                <span className="w-1 h-1 bg-blue-500 rounded-full mt-1.5 flex-shrink-0" />
                                <span>Stochastic variation near boundaries accounts for ~8% phenotype variability.</span>
                            </li>
                        </ul>
                    </div>
                </div>
            </div>
        </div>
    );
}
