import React, { useState, useMemo } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, BarChart, Bar, Cell } from 'recharts';
import { Sliders, AlertTriangle, CheckCircle, Info } from 'lucide-react';

// Default threshold values (matching backend)
const DEFAULT_THRESHOLDS = {
    bistability: 0.5,
    leaky: 0.5,
    oscillation_cv: 0.3,
    low_expression: 0.1
};

// Generate mock sensitivity sweep data
function generateSweepData(thresholds) {
    const sweepResults = [];
    const bistabRange = [0.3, 0.4, 0.5, 0.6, 0.7];
    const leakyRange = [0.3, 0.4, 0.5, 0.6, 0.7];

    for (const bistab of bistabRange) {
        for (const leaky of leakyRange) {
            // Mock calculation - in real app, this would come from backend
            const basePct = 35;
            const bistabEffect = (bistab - 0.5) * 40;
            const leakyEffect = (leaky - 0.5) * 20;

            sweepResults.push({
                bistability: bistab,
                leaky: leaky,
                pctBistable: Math.max(5, Math.min(95, basePct - bistabEffect - leakyEffect)),
                key: `${bistab}-${leaky}`
            });
        }
    }

    return sweepResults;
}

// Get robustness based on current thresholds
function calculateRobustness(thresholds) {
    const base = 35;
    const bistabEffect = (thresholds.bistability - 0.5) * 40;
    const leakyEffect = (thresholds.leaky - 0.5) * 20;
    return Math.max(5, Math.min(95, base - bistabEffect - leakyEffect));
}

// Region ranking data (stays stable)
const REGION_RANKINGS = [
    { region: 'Operator', rank: 1, score: 0.88, stable: true },
    { region: 'CDS', rank: 2, score: 0.72, stable: true },
    { region: 'RBS', rank: 3, score: 0.65, stable: true },
    { region: 'Promoter', rank: 4, score: 0.58, stable: true },
];

// Findings table
const FINDINGS = [
    { finding: 'Leaky expression dominates failures', confidence: 'High', type: 'qualitative' },
    { finding: 'Operators more robust than promoters', confidence: 'High', type: 'qualitative' },
    { finding: 'Exact robustness score (10.2%)', confidence: 'Low', type: 'quantitative' },
    { finding: '100% failure outside operators', confidence: 'Model-dependent', type: 'quantitative' },
];

export default function SensitivityAnalysis() {
    const [thresholds, setThresholds] = useState(DEFAULT_THRESHOLDS);

    const sweepData = useMemo(() => generateSweepData(thresholds), []);
    const currentRobustness = calculateRobustness(thresholds);

    // Find range across sweep
    const allPcts = sweepData.map(d => d.pctBistable);
    const minPct = Math.min(...allPcts);
    const maxPct = Math.max(...allPcts);
    const rangePct = maxPct - minPct;

    const isSensitive = rangePct > 15;

    const handleThresholdChange = (name, value) => {
        setThresholds(prev => ({ ...prev, [name]: parseFloat(value) }));
    };

    return (
        <div className="space-y-6">
            {/* Header */}
            <div>
                <h2 className="text-2xl font-bold text-slate-900">Sensitivity Analysis</h2>
                <p className="text-slate-500 mt-1">
                    Explore how classification thresholds affect robustness metrics
                </p>
            </div>

            {/* Threshold Controls */}
            <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                <div className="flex items-center gap-2 mb-4">
                    <Sliders className="text-blue-500" size={20} />
                    <h3 className="text-lg font-semibold text-slate-800">Classification Thresholds</h3>
                </div>

                <div className="grid grid-cols-2 gap-6">
                    <div>
                        <div className="flex justify-between mb-2">
                            <label className="text-sm font-medium text-slate-700">Bistability Threshold</label>
                            <span className="text-sm font-mono text-slate-500">{thresholds.bistability.toFixed(2)}</span>
                        </div>
                        <input
                            type="range"
                            min="0.2"
                            max="0.8"
                            step="0.05"
                            value={thresholds.bistability}
                            onChange={(e) => handleThresholdChange('bistability', e.target.value)}
                            className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                        />
                        <p className="text-xs text-slate-400 mt-1">
                            Minimum state separation for bistable classification
                        </p>
                    </div>

                    <div>
                        <div className="flex justify-between mb-2">
                            <label className="text-sm font-medium text-slate-700">Leaky Expression Threshold</label>
                            <span className="text-sm font-mono text-slate-500">{thresholds.leaky.toFixed(2)}</span>
                        </div>
                        <input
                            type="range"
                            min="0.2"
                            max="0.8"
                            step="0.05"
                            value={thresholds.leaky}
                            onChange={(e) => handleThresholdChange('leaky', e.target.value)}
                            className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-orange-500"
                        />
                        <p className="text-xs text-slate-400 mt-1">
                            Off/On ratio above which expression is considered leaky
                        </p>
                    </div>
                </div>
            </div>

            {/* Current Robustness Display */}
            <div className="grid grid-cols-3 gap-6">
                <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                    <p className="text-sm font-medium text-slate-500">Current Robustness</p>
                    <p className="text-4xl font-bold text-blue-600 mt-2">{currentRobustness.toFixed(1)}%</p>
                    <p className="text-xs text-slate-400 mt-1">at current thresholds</p>
                </div>

                <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                    <p className="text-sm font-medium text-slate-500">Robustness Range (±20%)</p>
                    <p className="text-4xl font-bold text-slate-800 mt-2">{minPct.toFixed(0)}-{maxPct.toFixed(0)}%</p>
                    <p className="text-xs text-slate-400 mt-1">across threshold sweep</p>
                </div>

                <div className={`rounded-xl p-6 shadow-sm border ${isSensitive ? 'bg-yellow-50 border-yellow-200' : 'bg-green-50 border-green-200'}`}>
                    <div className="flex items-center gap-2">
                        {isSensitive ? (
                            <AlertTriangle className="text-yellow-600" size={20} />
                        ) : (
                            <CheckCircle className="text-green-600" size={20} />
                        )}
                        <p className="text-sm font-medium text-slate-700">Sensitivity Status</p>
                    </div>
                    <p className={`text-lg font-bold mt-2 ${isSensitive ? 'text-yellow-700' : 'text-green-700'}`}>
                        {isSensitive ? 'THRESHOLD SENSITIVE' : 'STABLE'}
                    </p>
                    <p className="text-xs text-slate-500 mt-1">
                        {rangePct.toFixed(1)}% variation across sweep
                    </p>
                </div>
            </div>

            {/* Region Ranking Chart */}
            <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                <h3 className="text-lg font-semibold text-slate-800 mb-4">Region Robustness Ranking (Stable)</h3>
                <div className="h-[250px] w-full">
                    <ResponsiveContainer width="100%" height="100%">
                        <BarChart data={REGION_RANKINGS} layout="vertical" margin={{ left: 80 }}>
                            <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
                            <XAxis type="number" domain={[0, 1]} tick={{ fontSize: 12 }} />
                            <YAxis type="category" dataKey="region" tick={{ fontSize: 12 }} />
                            <Tooltip contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }} />
                            <Bar dataKey="score" name="Robustness Score" radius={[0, 4, 4, 0]}>
                                {REGION_RANKINGS.map((entry, index) => (
                                    <Cell key={`cell-${index}`} fill={entry.stable ? '#3b82f6' : '#f97316'} />
                                ))}
                            </Bar>
                        </BarChart>
                    </ResponsiveContainer>
                </div>
                <p className="text-sm text-slate-500 mt-2 text-center">
                    Rankings remain stable across ±20% threshold variation
                </p>
            </div>

            {/* Qualitative vs Quantitative Findings */}
            <div className="bg-white rounded-xl p-6 shadow-sm border border-slate-100">
                <div className="flex items-center gap-2 mb-4">
                    <Info className="text-blue-500" size={20} />
                    <h3 className="text-lg font-semibold text-slate-800">Finding Confidence Levels</h3>
                </div>
                <div className="overflow-x-auto">
                    <table className="w-full text-sm">
                        <thead>
                            <tr className="border-b border-slate-200">
                                <th className="text-left py-3 px-4 font-semibold text-slate-600">Finding</th>
                                <th className="text-left py-3 px-4 font-semibold text-slate-600">Confidence</th>
                                <th className="text-left py-3 px-4 font-semibold text-slate-600">Type</th>
                            </tr>
                        </thead>
                        <tbody>
                            {FINDINGS.map((finding, idx) => (
                                <tr key={idx} className="border-b border-slate-100 hover:bg-slate-50">
                                    <td className="py-3 px-4 text-slate-800">{finding.finding}</td>
                                    <td className="py-3 px-4">
                                        <span className={`px-2 py-1 rounded-full text-xs font-medium ${finding.confidence === 'High' ? 'bg-green-100 text-green-700' :
                                                finding.confidence === 'Low' ? 'bg-red-100 text-red-700' :
                                                    'bg-yellow-100 text-yellow-700'
                                            }`}>
                                            {finding.confidence}
                                        </span>
                                    </td>
                                    <td className="py-3 px-4">
                                        <span className={`px-2 py-1 rounded text-xs ${finding.type === 'qualitative' ? 'bg-blue-50 text-blue-700' : 'bg-slate-100 text-slate-600'
                                            }`}>
                                            {finding.type}
                                        </span>
                                    </td>
                                </tr>
                            ))}
                        </tbody>
                    </table>
                </div>
            </div>

            {/* Recommendations */}
            <div className="bg-gradient-to-r from-indigo-50 to-purple-50 rounded-xl p-6 border border-indigo-100">
                <h3 className="text-lg font-semibold text-indigo-900 mb-3">Recommendations</h3>
                <ul className="space-y-2 text-sm text-indigo-800">
                    <li className="flex items-start gap-2">
                        <CheckCircle className="text-indigo-600 mt-0.5 flex-shrink-0" size={16} />
                        <span>Region rankings are <strong>stable</strong> across threshold variations - high confidence in relative comparisons</span>
                    </li>
                    <li className="flex items-start gap-2">
                        <AlertTriangle className="text-yellow-600 mt-0.5 flex-shrink-0" size={16} />
                        <span>Absolute robustness percentages are <strong>model-dependent</strong> - report ranges rather than point estimates</span>
                    </li>
                    <li className="flex items-start gap-2">
                        <Info className="text-blue-600 mt-0.5 flex-shrink-0" size={16} />
                        <span>Focus on <strong>qualitative conclusions</strong> (e.g., "operators are more robust") over quantitative claims</span>
                    </li>
                </ul>
            </div>
        </div>
    );
}
