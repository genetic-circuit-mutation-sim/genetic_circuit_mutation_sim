import React, { useEffect, useState, useMemo } from 'react';
import { fetchResultsCsv } from '../lib/data';
import { ScatterChart, Scatter, XAxis, YAxis, ZAxis, Tooltip, ResponsiveContainer, Cell, ReferenceLine } from 'recharts';

export default function Visualizations() {
    const [data, setData] = useState([]);
    const [loading, setLoading] = useState(true);

    useEffect(() => {
        fetchResultsCsv()
            .then(setData)
            .finally(() => setLoading(false));
    }, []);

    // Prepare Heatmap Data (Matrix of Region vs Failure Mode)
    const heatmapData = useMemo(() => {
        if (!data.length) return [];

        const regions = [...new Set(data.map(d => d.region))].filter(Boolean);
        const modes = ['bistable', 'leaky', 'loss_of_bistability', 'no_expression', 'oscillatory'];

        const matrix = [];
        regions.forEach(region => {
            const regionRows = data.filter(d => d.region === region);
            modes.forEach(mode => {
                const count = regionRows.filter(d => d.failure_label === mode).length;
                const rate = count / regionRows.length;
                if (rate > 0) {
                    matrix.push({ x: mode.replace(/_/g, '\n'), y: region, z: rate });
                }
            });
        });
        return matrix;
    }, [data]);

    // Position failure rate
    const positionData = useMemo(() => {
        if (!data.length) return [];
        // Flatten all mutation positions
        const allPositions = [];
        data.forEach(row => {
            if (row.mutation_positions && typeof row.mutation_positions === 'string') {
                // If it's a string representation from CSV like "[1, 2]"
                try {
                    const pos = JSON.parse(row.mutation_positions);
                    if (Array.isArray(pos)) allPositions.push(...pos);
                } catch (e) { }
            }
        });

        // Count failures per position (heuristic: assume if variants fail, positions involved "failed")
        // This is a simplified view compared to the Python script which maps index to rate precisely
        // Since we don't have the exact gene length easily, we'll just plot density of mutations in failed variants

        const failedRows = data.filter(d => d.failure_label !== 'bistable');
        const failedPositions = [];
        failedRows.forEach(row => {
            if (row.mutation_positions && typeof row.mutation_positions === 'string') {
                try {
                    const pos = JSON.parse(row.mutation_positions);
                    if (Array.isArray(pos)) failedPositions.push(...pos);
                } catch (e) { }
            }
        });

        // Binning
        const bins = {};
        failedPositions.forEach(p => {
            const bin = Math.floor(p / 10) * 10;
            bins[bin] = (bins[bin] || 0) + 1;
        });

        return Object.entries(bins).map(([bin, count]) => ({
            bin: parseInt(bin),
            count
        })).sort((a, b) => a.bin - b.bin);

    }, [data]);


    if (loading) return <div className="p-8 text-slate-500">Loading visualization data...</div>;

    return (
        <div className="space-y-8">
            <div>
                <h2 className="text-2xl font-bold text-slate-900">Advanced Visualizations</h2>
                <p className="text-slate-500 mt-1">Deep dive into failure modes and regional sensitivity.</p>
            </div>

            <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
                {/* Heatmap-like Scatter Plot */}
                <div className="bg-white p-6 rounded-xl shadow-sm border border-slate-100">
                    <h3 className="font-bold text-slate-800 mb-6">Region vs Failure Mode</h3>
                    <div className="h-[400px]">
                        <ResponsiveContainer width="100%" height="100%">
                            <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 40 }}>
                                <XAxis type="category" dataKey="x" name="Mode" tick={{ fontSize: 10 }} />
                                <YAxis type="category" dataKey="y" name="Region" width={80} tick={{ fontSize: 10 }} />
                                <ZAxis type="number" dataKey="z" range={[0, 1000]} name="Rate" />
                                <Tooltip cursor={{ strokeDasharray: '3 3' }} />
                                <Scatter data={heatmapData} shape="circle">
                                    {heatmapData.map((entry, index) => (
                                        <Cell key={`cell-${index}`} fill={entry.z > 0.5 ? '#ef4444' : '#f97316'} fillOpacity={entry.z} />
                                    ))}
                                </Scatter>
                            </ScatterChart>
                        </ResponsiveContainer>
                    </div>
                </div>

                {/* Simulated Mutation Density */}
                <div className="bg-white p-6 rounded-xl shadow-sm border border-slate-100">
                    <h3 className="font-bold text-slate-800 mb-6">Mutation Density in Failed Variants</h3>
                    <div className="h-[400px]">
                        <ResponsiveContainer width="100%" height="100%">
                            <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 10 }}>
                                <XAxis dataKey="bin" name="Position" unit="bp" />
                                <YAxis dataKey="count" name="Frequency" />
                                <Tooltip cursor={{ strokeDasharray: '3 3' }} />
                                <Scatter name="Failures" data={positionData} fill="#8884d8" line shape="circle" />
                            </ScatterChart>
                        </ResponsiveContainer>
                    </div>
                </div>
            </div>
        </div>
    );
}
