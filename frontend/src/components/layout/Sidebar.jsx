import React from 'react';
import { LayoutDashboard, Play, FileText, Activity, TrendingDown, Sliders, Layers, ShieldCheck } from 'lucide-react';
import { Link, useLocation } from 'react-router-dom';
import clsx from 'clsx';

const NAV_ITEMS = [
    { label: 'Dashboard', path: '/', icon: LayoutDashboard },
    { label: 'Simulator', path: '/simulator', icon: Play },
    { label: 'Evolutionary Walks', path: '/evolutionary-walks', icon: TrendingDown },
    { label: 'Sensitivity Analysis', path: '/sensitivity', icon: Sliders },
    { label: 'Parameter Space', path: '/parameter-space', icon: Layers },
    { label: 'Validation', path: '/validation', icon: ShieldCheck },
    { label: 'Visualizations', path: '/visualizations', icon: Activity },
    { label: 'Results Explorer', path: '/results', icon: FileText },
];

export default function Sidebar() {
    const location = useLocation();

    return (
        <aside className="w-64 bg-slate-900 text-white min-h-screen fixed left-0 top-0 flex flex-col border-r border-slate-800">
            <div className="p-6 border-b border-slate-800">
                <h1 className="text-xl font-bold bg-gradient-to-r from-blue-400 to-emerald-400 bg-clip-text text-transparent">
                    BioCircuit Sim
                </h1>
                <p className="text-xs text-slate-400 mt-1">Mutation Robustness Analysis</p>
            </div>

            <nav className="flex-1 py-6 px-3 space-y-1">
                {NAV_ITEMS.map((item) => {
                    const isActive = location.pathname === item.path;
                    const Icon = item.icon;
                    return (
                        <Link
                            key={item.path}
                            to={item.path}
                            className={clsx(
                                "flex items-center gap-3 px-4 py-3 rounded-lg text-sm font-medium transition-colors",
                                isActive
                                    ? "bg-slate-800 text-white border border-slate-700"
                                    : "text-slate-400 hover:bg-slate-800/50 hover:text-white"
                            )}
                        >
                            <Icon size={18} className={clsx(isActive ? "text-emerald-400" : "text-slate-500 group-hover:text-slate-400")} />
                            {item.label}
                        </Link>
                    );
                })}
            </nav>

            <div className="p-4 border-t border-slate-800">
                <div className="text-xs text-slate-500 text-center">
                    v1.0.0 Alpha
                </div>
            </div>
        </aside>
    );
}
