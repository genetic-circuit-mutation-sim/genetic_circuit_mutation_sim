import React, { useEffect, useState, useMemo } from 'react';
import { fetchResultsCsv } from '../lib/data';
import { Search, ChevronLeft, ChevronRight, AlertCircle, CheckCircle } from 'lucide-react';
import clsx from 'clsx';

export default function ResultsExplorer() {
    const [data, setData] = useState([]);
    const [loading, setLoading] = useState(true);
    const [search, setSearch] = useState('');
    const [currentPage, setCurrentPage] = useState(1);
    const rowsPerPage = 20;

    useEffect(() => {
        fetchResultsCsv()
            .then(setData)
            .finally(() => setLoading(false));
    }, []);

    const filteredData = useMemo(() => {
        return data.filter(row =>
            row.region?.toLowerCase().includes(search.toLowerCase()) ||
            row.failure_label?.toLowerCase().includes(search.toLowerCase()) ||
            row.mutation_type?.toLowerCase().includes(search.toLowerCase())
        );
    }, [data, search]);

    const totalPages = Math.ceil(filteredData.length / rowsPerPage);
    const paginatedData = filteredData.slice((currentPage - 1) * rowsPerPage, currentPage * rowsPerPage);

    if (loading) return <div className="p-8 text-slate-500">Loading results data...</div>;

    return (
        <div className="space-y-6">
            <div className="flex justify-between items-center">
                <div>
                    <h2 className="text-2xl font-bold text-slate-900">Results Explorer</h2>
                    <p className="text-slate-500 text-sm mt-1">
                        Viewing {filteredData.length} simulation records
                    </p>
                </div>
                <div className="relative">
                    <Search className="absolute left-3 top-1/2 -translate-y-1/2 text-slate-400" size={18} />
                    <input
                        type="text"
                        placeholder="Search region, type, or failure mode..."
                        value={search}
                        onChange={(e) => setSearch(e.target.value)}
                        className="pl-10 pr-4 py-2 border border-slate-200 rounded-lg text-sm w-80 focus:outline-none focus:ring-2 focus:ring-blue-500/20 focus:border-blue-500"
                    />
                </div>
            </div>

            <div className="bg-white rounded-xl shadow-sm border border-slate-100 overflow-hidden">
                <div className="overflow-x-auto">
                    <table className="w-full text-sm text-left">
                        <thead className="bg-slate-50 border-b border-slate-100 text-slate-500 font-medium">
                            <tr>
                                <th className="px-6 py-4">ID</th>
                                <th className="px-6 py-4">Region</th>
                                <th className="px-6 py-4">Type</th>
                                <th className="px-6 py-4">Mutations</th>
                                <th className="px-6 py-4">Params (Relative)</th>
                                <th className="px-6 py-4">Outcome</th>
                            </tr>
                        </thead>
                        <tbody className="divide-y divide-slate-50">
                            {paginatedData.map((row) => (
                                <tr key={row.variant_id} className="hover:bg-slate-50/50">
                                    <td className="px-6 py-4 font-mono text-slate-400">#{row.variant_id}</td>
                                    <td className="px-6 py-4 font-medium text-slate-700">{row.region || '-'}</td>
                                    <td className="px-6 py-4">
                                        <span className={clsx(
                                            "inline-flex items-center px-2 py-1 rounded-full text-xs font-medium capitalize",
                                            row.mutation_type === 'substitution' ? "bg-blue-50 text-blue-700" :
                                                row.mutation_type === 'deletion' ? "bg-red-50 text-red-700" :
                                                    row.mutation_type === 'insertion' ? "bg-amber-50 text-amber-700" :
                                                        "bg-slate-100 text-slate-600"
                                        )}>
                                            {row.mutation_type || 'N/A'}
                                        </span>
                                    </td>
                                    <td className="px-6 py-4 text-slate-600 font-mono text-xs max-w-xs truncate" title={row.mutations}>
                                        {row.mutations}
                                    </td>
                                    <td className="px-6 py-4">
                                        <div className="flex flex-col gap-1 text-xs text-slate-500">
                                            <span>tx_A: {row.param_tx_A?.toFixed(2)}</span>
                                            <span>Kd: {row.param_Kd?.toFixed(2)}</span>
                                        </div>
                                    </td>
                                    <td className="px-6 py-4">
                                        <div className="flex items-center gap-2">
                                            {row.failure_label === 'bistable' ? (
                                                <CheckCircle size={16} className="text-emerald-500" />
                                            ) : (
                                                <AlertCircle size={16} className="text-red-500" />
                                            )}
                                            <span className={clsx(
                                                "font-medium capitalize",
                                                row.failure_label === 'bistable' ? "text-emerald-700" : "text-slate-600"
                                            )}>
                                                {row.failure_label?.replace(/_/g, ' ')}
                                            </span>
                                        </div>
                                    </td>
                                </tr>
                            ))}
                        </tbody>
                    </table>
                </div>

                <div className="px-6 py-4 border-t border-slate-100 flex items-center justify-between">
                    <div className="text-sm text-slate-500">
                        Page {currentPage} of {totalPages}
                    </div>
                    <div className="flex gap-2">
                        <button
                            disabled={currentPage === 1}
                            onClick={() => setCurrentPage(c => c - 1)}
                            className="p-1 rounded hover:bg-slate-100 disabled:opacity-50 disabled:cursor-not-allowed"
                        >
                            <ChevronLeft size={20} className="text-slate-500" />
                        </button>
                        <button
                            disabled={currentPage === totalPages}
                            onClick={() => setCurrentPage(c => c + 1)}
                            className="p-1 rounded hover:bg-slate-100 disabled:opacity-50 disabled:cursor-not-allowed"
                        >
                            <ChevronRight size={20} className="text-slate-500" />
                        </button>
                    </div>
                </div>
            </div>
        </div>
    );
}
