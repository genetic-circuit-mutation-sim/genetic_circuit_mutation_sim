import React from 'react';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import Layout from './components/layout/Layout';
import Dashboard from './pages/Dashboard';
import ResultsSimulator from './pages/ResultsSimulator';
import ResultsExplorer from './pages/ResultsExplorer';
import Visualizations from './pages/Visualizations';
import EvolutionaryWalks from './pages/EvolutionaryWalks';
import SensitivityAnalysis from './pages/SensitivityAnalysis';
import ParameterSpaceExplorer from './pages/ParameterSpaceExplorer';
import ValidationBenchmarks from './pages/ValidationBenchmarks';

function App() {
  return (
    <Router>
      <Layout>
        <Routes>
          <Route path="/" element={<Dashboard />} />
          <Route path="/simulator" element={<ResultsSimulator />} />
          <Route path="/visualizations" element={<Visualizations />} />
          <Route path="/results" element={<ResultsExplorer />} />
          <Route path="/evolutionary-walks" element={<EvolutionaryWalks />} />
          <Route path="/sensitivity" element={<SensitivityAnalysis />} />
          <Route path="/parameter-space" element={<ParameterSpaceExplorer />} />
          <Route path="/validation" element={<ValidationBenchmarks />} />
        </Routes>
      </Layout>
    </Router>
  );
}

export default App;
