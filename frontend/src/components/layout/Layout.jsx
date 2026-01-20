import React from 'react';
import Sidebar from './Sidebar';

export default function Layout({ children }) {
    return (
        <div className="flex min-h-screen bg-slate-50">
            <Sidebar />
            <main className="flex-1 ml-64 p-8 flex flex-col">
                <div className="max-w-7xl mx-auto flex-1 w-full">
                    {children}
                </div>

                <footer className="mt-12 pt-8 border-t border-slate-200 text-slate-400 text-xs pb-4">
                    <div className="max-w-7xl mx-auto flex justify-between items-center px-4">
                        <div className="max-w-2xl">
                            <strong>Model Disclaimer:</strong> This tool provides model-based insights for <em>in silico</em> genetic circuit design.
                            Results represent biophysical predictions and should be validated with laboratory experiments.
                            Qualitative rankings are more robust than absolute numerical percentages.
                        </div>
                        <div className="flex gap-4">
                            <a href="/docs/reproducibility.md" className="hover:text-blue-500 underline transition-colors">Reproducibility Appendix</a>
                            <span>BioCircuit v1.0 Alpha</span>
                        </div>
                    </div>
                </footer>
            </main>
        </div>
    );
}
