import Papa from 'papaparse';

export const fetchRobustnessSummary = async () => {
    const response = await fetch('/data/robustness_summary.json');
    if (!response.ok) {
        throw new Error('Failed to fetch robustness summary');
    }
    return await response.json();
};

export const fetchResultsCsv = async () => {
    return new Promise((resolve, reject) => {
        Papa.parse('/data/results.csv', {
            download: true,
            header: true,
            dynamicTyping: true,
            skipEmptyLines: true,
            complete: (results) => {
                resolve(results.data);
            },
            error: (error) => {
                reject(error);
            }
        });
    });
};
export const fetchEvolutionaryWalks = async () => {
    // Mock for now, will fetch from /data/evolutionary_walks.json
    try {
        const response = await fetch('/data/evolutionary_walks.json');
        if (response.ok) return await response.json();
    } catch (e) { console.warn('Using mock walks data'); }
    return { walks: [], aggregate: [] };
};

export const fetchSensitivityAnalysis = async () => {
    try {
        const response = await fetch('/data/sensitivity_analysis.json');
        if (response.ok) return await response.json();
    } catch (e) { console.warn('Using mock sensitivity data'); }
    return { sweep_results: [] };
};

export const fetchValidationResults = async () => {
    try {
        const response = await fetch('/data/validation_report.json');
        if (response.ok) return await response.json();
    } catch (e) { console.warn('Using mock validation data'); }
    return { results: [], overall_status: 'NOT_RUN' };
};
