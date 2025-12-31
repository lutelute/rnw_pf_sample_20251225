import React, { useState, useMemo, useCallback } from 'react';

// ============================================
// MATPOWER-style Power Flow Calculator
// Newton-Raphson Method Implementation
// ============================================

// Complex number operations
const Complex = {
  add: (a, b) => ({ re: a.re + b.re, im: a.im + b.im }),
  sub: (a, b) => ({ re: a.re - b.re, im: a.im - b.im }),
  mul: (a, b) => ({
    re: a.re * b.re - a.im * b.im,
    im: a.re * b.im + a.im * b.re
  }),
  div: (a, b) => {
    const denom = b.re * b.re + b.im * b.im;
    return {
      re: (a.re * b.re + a.im * b.im) / denom,
      im: (a.im * b.re - a.re * b.im) / denom
    };
  },
  conj: (a) => ({ re: a.re, im: -a.im }),
  abs: (a) => Math.sqrt(a.re * a.re + a.im * a.im),
  fromPolar: (mag, angle) => ({
    re: mag * Math.cos(angle),
    im: mag * Math.sin(angle)
  }),
  zero: () => ({ re: 0, im: 0 })
};

// Matrix operations
const Matrix = {
  zeros: (rows, cols) => Array(rows).fill(null).map(() => Array(cols).fill(0)),
  
  solve: (A, b) => {
    // Gaussian elimination with partial pivoting
    const n = A.length;
    const augmented = A.map((row, i) => [...row, b[i]]);
    
    for (let col = 0; col < n; col++) {
      // Find pivot
      let maxRow = col;
      for (let row = col + 1; row < n; row++) {
        if (Math.abs(augmented[row][col]) > Math.abs(augmented[maxRow][col])) {
          maxRow = row;
        }
      }
      [augmented[col], augmented[maxRow]] = [augmented[maxRow], augmented[col]];
      
      if (Math.abs(augmented[col][col]) < 1e-12) continue;
      
      // Eliminate
      for (let row = col + 1; row < n; row++) {
        const factor = augmented[row][col] / augmented[col][col];
        for (let j = col; j <= n; j++) {
          augmented[row][j] -= factor * augmented[col][j];
        }
      }
    }
    
    // Back substitution
    const x = Array(n).fill(0);
    for (let row = n - 1; row >= 0; row--) {
      let sum = augmented[row][n];
      for (let col = row + 1; col < n; col++) {
        sum -= augmented[row][col] * x[col];
      }
      x[row] = Math.abs(augmented[row][row]) > 1e-12 ? sum / augmented[row][row] : 0;
    }
    return x;
  }
};

// MATPOWER-style bus types
const BUS_TYPE = {
  PQ: 1,    // Load bus
  PV: 2,    // Generator bus
  SLACK: 3  // Slack/Swing bus
};

// Default IEEE 9-bus system (simplified)
const defaultBusData = [
  { bus: 1, type: BUS_TYPE.SLACK, Pd: 0, Qd: 0, Vm: 1.04, Va: 0, baseKV: 16.5, name: 'Gen1' },
  { bus: 2, type: BUS_TYPE.PV, Pd: 0, Qd: 0, Vm: 1.025, Va: 0, baseKV: 18, name: 'Gen2' },
  { bus: 3, type: BUS_TYPE.PV, Pd: 0, Qd: 0, Vm: 1.025, Va: 0, baseKV: 13.8, name: 'Gen3' },
  { bus: 4, type: BUS_TYPE.PQ, Pd: 0, Qd: 0, Vm: 1.0, Va: 0, baseKV: 230, name: 'Bus4' },
  { bus: 5, type: BUS_TYPE.PQ, Pd: 125, Qd: 50, Vm: 1.0, Va: 0, baseKV: 230, name: 'Load A' },
  { bus: 6, type: BUS_TYPE.PQ, Pd: 90, Qd: 30, Vm: 1.0, Va: 0, baseKV: 230, name: 'Load B' },
  { bus: 7, type: BUS_TYPE.PQ, Pd: 0, Qd: 0, Vm: 1.0, Va: 0, baseKV: 230, name: 'Bus7' },
  { bus: 8, type: BUS_TYPE.PQ, Pd: 100, Qd: 35, Vm: 1.0, Va: 0, baseKV: 230, name: 'Load C' },
  { bus: 9, type: BUS_TYPE.PQ, Pd: 0, Qd: 0, Vm: 1.0, Va: 0, baseKV: 230, name: 'Bus9' },
];

const defaultBranchData = [
  { from: 1, to: 4, r: 0, x: 0.0576, b: 0, tap: 1 },
  { from: 2, to: 7, r: 0, x: 0.0625, b: 0, tap: 1 },
  { from: 3, to: 9, r: 0, x: 0.0586, b: 0, tap: 1 },
  { from: 4, to: 5, r: 0.01, x: 0.085, b: 0.176, tap: 1 },
  { from: 4, to: 6, r: 0.017, x: 0.092, b: 0.158, tap: 1 },
  { from: 5, to: 7, r: 0.032, x: 0.161, b: 0.306, tap: 1 },
  { from: 6, to: 9, r: 0.039, x: 0.170, b: 0.358, tap: 1 },
  { from: 7, to: 8, r: 0.0085, x: 0.072, b: 0.149, tap: 1 },
  { from: 8, to: 9, r: 0.0119, x: 0.1008, b: 0.209, tap: 1 },
];

const defaultGenData = [
  { bus: 1, Pg: 0, Qg: 0, Qmax: 300, Qmin: -300, Vg: 1.04, mBase: 100, status: 1 },
  { bus: 2, Pg: 163, Qg: 0, Qmax: 300, Qmin: -300, Vg: 1.025, mBase: 100, status: 1 },
  { bus: 3, Pg: 85, Qg: 0, Qmax: 300, Qmin: -300, Vg: 1.025, mBase: 100, status: 1 },
];

// Build Ybus matrix
function buildYbus(busData, branchData) {
  const n = busData.length;
  const Y = Array(n).fill(null).map(() => Array(n).fill(null).map(() => Complex.zero()));
  
  branchData.forEach(branch => {
    const i = branch.from - 1;
    const k = branch.to - 1;
    const tap = branch.tap || 1;
    
    // Series admittance
    const z = { re: branch.r, im: branch.x };
    const y = Complex.div({ re: 1, im: 0 }, z);
    
    // Shunt admittance (line charging)
    const bShunt = branch.b / 2;
    
    // Π model with tap ratio
    // Y_ii += y/tap^2 + jb/2
    Y[i][i] = Complex.add(Y[i][i], { 
      re: y.re / (tap * tap), 
      im: y.im / (tap * tap) + bShunt 
    });
    
    // Y_kk += y + jb/2
    Y[k][k] = Complex.add(Y[k][k], { 
      re: y.re, 
      im: y.im + bShunt 
    });
    
    // Y_ik = Y_ki = -y/tap
    const yMutual = { re: -y.re / tap, im: -y.im / tap };
    Y[i][k] = Complex.add(Y[i][k], yMutual);
    Y[k][i] = Complex.add(Y[k][i], yMutual);
  });
  
  return Y;
}

// Calculate power injections
function calcPowerInjection(V, Y) {
  const n = V.length;
  const S = [];
  
  for (let i = 0; i < n; i++) {
    let sum = Complex.zero();
    for (let k = 0; k < n; k++) {
      const Vik = Complex.mul(V[i], Complex.conj(V[k]));
      sum = Complex.add(sum, Complex.mul(Y[i][k], Complex.conj(V[k])));
    }
    const Si = Complex.mul(V[i], Complex.conj(sum));
    S.push(Si);
  }
  
  return S;
}

// Newton-Raphson Power Flow
function runPowerFlow(busData, branchData, genData, options = {}) {
  const { maxIter = 50, tolerance = 1e-6 } = options;
  const n = busData.length;
  const baseMVA = 100;
  
  // Build Ybus
  const Y = buildYbus(busData, branchData);
  
  // Initialize voltage
  const V = busData.map(bus => Complex.fromPolar(bus.Vm, bus.Va * Math.PI / 180));
  const Vm = busData.map(bus => bus.Vm);
  const Va = busData.map(bus => bus.Va * Math.PI / 180);
  
  // Specified power
  const Pspec = busData.map(bus => {
    const gen = genData.find(g => g.bus === bus.bus);
    const Pg = gen ? gen.Pg / baseMVA : 0;
    const Pd = bus.Pd / baseMVA;
    return Pg - Pd;
  });
  
  const Qspec = busData.map(bus => {
    const gen = genData.find(g => g.bus === bus.bus);
    const Qg = gen ? gen.Qg / baseMVA : 0;
    const Qd = bus.Qd / baseMVA;
    return Qg - Qd;
  });
  
  // Find PQ and PV buses
  const pqBuses = busData.map((bus, i) => bus.type === BUS_TYPE.PQ ? i : -1).filter(i => i >= 0);
  const pvBuses = busData.map((bus, i) => bus.type === BUS_TYPE.PV ? i : -1).filter(i => i >= 0);
  const nonSlackBuses = [...pqBuses, ...pvBuses].sort((a, b) => a - b);
  
  const history = [];
  let converged = false;
  let iter = 0;
  
  // Newton-Raphson iteration
  for (iter = 0; iter < maxIter; iter++) {
    // Update voltage phasors
    for (let i = 0; i < n; i++) {
      V[i] = Complex.fromPolar(Vm[i], Va[i]);
    }
    
    // Calculate power injections
    const Pcalc = [];
    const Qcalc = [];
    
    for (let i = 0; i < n; i++) {
      let Pi = 0, Qi = 0;
      for (let k = 0; k < n; k++) {
        const Gik = Y[i][k].re;
        const Bik = Y[i][k].im;
        const angleIK = Va[i] - Va[k];
        Pi += Vm[i] * Vm[k] * (Gik * Math.cos(angleIK) + Bik * Math.sin(angleIK));
        Qi += Vm[i] * Vm[k] * (Gik * Math.sin(angleIK) - Bik * Math.cos(angleIK));
      }
      Pcalc.push(Pi);
      Qcalc.push(Qi);
    }
    
    // Power mismatch
    const dP = nonSlackBuses.map(i => Pspec[i] - Pcalc[i]);
    const dQ = pqBuses.map(i => Qspec[i] - Qcalc[i]);
    const mismatch = [...dP, ...dQ];
    
    const maxMismatch = Math.max(...mismatch.map(Math.abs));
    history.push({ iter, maxMismatch, Vm: [...Vm], Va: Va.map(a => a * 180 / Math.PI) });
    
    if (maxMismatch < tolerance) {
      converged = true;
      break;
    }
    
    // Build Jacobian
    const nP = nonSlackBuses.length;
    const nQ = pqBuses.length;
    const J = Matrix.zeros(nP + nQ, nP + nQ);
    
    // J1 = dP/dθ
    for (let ii = 0; ii < nP; ii++) {
      const i = nonSlackBuses[ii];
      for (let kk = 0; kk < nP; kk++) {
        const k = nonSlackBuses[kk];
        if (i === k) {
          let sum = 0;
          for (let j = 0; j < n; j++) {
            if (j !== i) {
              const Gij = Y[i][j].re;
              const Bij = Y[i][j].im;
              sum += Vm[j] * (Gij * Math.sin(Va[i] - Va[j]) - Bij * Math.cos(Va[i] - Va[j]));
            }
          }
          J[ii][kk] = -Vm[i] * sum;
        } else {
          const Gik = Y[i][k].re;
          const Bik = Y[i][k].im;
          J[ii][kk] = Vm[i] * Vm[k] * (Gik * Math.sin(Va[i] - Va[k]) - Bik * Math.cos(Va[i] - Va[k]));
        }
      }
    }
    
    // J2 = dP/dV
    for (let ii = 0; ii < nP; ii++) {
      const i = nonSlackBuses[ii];
      for (let kk = 0; kk < nQ; kk++) {
        const k = pqBuses[kk];
        if (i === k) {
          let sum = 0;
          for (let j = 0; j < n; j++) {
            const Gij = Y[i][j].re;
            const Bij = Y[i][j].im;
            sum += Vm[j] * (Gij * Math.cos(Va[i] - Va[j]) + Bij * Math.sin(Va[i] - Va[j]));
          }
          J[ii][nP + kk] = sum + Y[i][i].re * Vm[i];
        } else {
          const Gik = Y[i][k].re;
          const Bik = Y[i][k].im;
          J[ii][nP + kk] = Vm[i] * (Gik * Math.cos(Va[i] - Va[k]) + Bik * Math.sin(Va[i] - Va[k]));
        }
      }
    }
    
    // J3 = dQ/dθ
    for (let ii = 0; ii < nQ; ii++) {
      const i = pqBuses[ii];
      for (let kk = 0; kk < nP; kk++) {
        const k = nonSlackBuses[kk];
        if (i === k) {
          let sum = 0;
          for (let j = 0; j < n; j++) {
            if (j !== i) {
              const Gij = Y[i][j].re;
              const Bij = Y[i][j].im;
              sum += Vm[j] * (Gij * Math.cos(Va[i] - Va[j]) + Bij * Math.sin(Va[i] - Va[j]));
            }
          }
          J[nP + ii][kk] = Vm[i] * sum;
        } else {
          const Gik = Y[i][k].re;
          const Bik = Y[i][k].im;
          J[nP + ii][kk] = -Vm[i] * Vm[k] * (Gik * Math.cos(Va[i] - Va[k]) + Bik * Math.sin(Va[i] - Va[k]));
        }
      }
    }
    
    // J4 = dQ/dV
    for (let ii = 0; ii < nQ; ii++) {
      const i = pqBuses[ii];
      for (let kk = 0; kk < nQ; kk++) {
        const k = pqBuses[kk];
        if (i === k) {
          let sum = 0;
          for (let j = 0; j < n; j++) {
            const Gij = Y[i][j].re;
            const Bij = Y[i][j].im;
            sum += Vm[j] * (Gij * Math.sin(Va[i] - Va[j]) - Bij * Math.cos(Va[i] - Va[j]));
          }
          J[nP + ii][nP + kk] = sum - Y[i][i].im * Vm[i];
        } else {
          const Gik = Y[i][k].re;
          const Bik = Y[i][k].im;
          J[nP + ii][nP + kk] = Vm[i] * (Gik * Math.sin(Va[i] - Va[k]) - Bik * Math.cos(Va[i] - Va[k]));
        }
      }
    }
    
    // Solve J * dx = mismatch
    const dx = Matrix.solve(J, mismatch);
    
    // Update variables
    for (let ii = 0; ii < nP; ii++) {
      Va[nonSlackBuses[ii]] += dx[ii];
    }
    for (let ii = 0; ii < nQ; ii++) {
      Vm[pqBuses[ii]] += dx[nP + ii];
    }
  }
  
  // Final power calculation
  for (let i = 0; i < n; i++) {
    V[i] = Complex.fromPolar(Vm[i], Va[i]);
  }
  
  const Pcalc = [];
  const Qcalc = [];
  for (let i = 0; i < n; i++) {
    let Pi = 0, Qi = 0;
    for (let k = 0; k < n; k++) {
      const Gik = Y[i][k].re;
      const Bik = Y[i][k].im;
      const angleIK = Va[i] - Va[k];
      Pi += Vm[i] * Vm[k] * (Gik * Math.cos(angleIK) + Bik * Math.sin(angleIK));
      Qi += Vm[i] * Vm[k] * (Gik * Math.sin(angleIK) - Bik * Math.cos(angleIK));
    }
    Pcalc.push(Pi * baseMVA);
    Qcalc.push(Qi * baseMVA);
  }
  
  // Line flows
  const lineFlows = branchData.map(branch => {
    const i = branch.from - 1;
    const k = branch.to - 1;
    const tap = branch.tap || 1;
    
    const z = { re: branch.r, im: branch.x };
    const y = Complex.div({ re: 1, im: 0 }, z);
    const bShunt = branch.b / 2;
    
    // Current from i to k
    const Vi = V[i];
    const Vk = V[k];
    
    const Iik_series = Complex.mul(y, Complex.sub(
      Complex.div(Vi, { re: tap, im: 0 }),
      Vk
    ));
    const Iik_shunt = Complex.mul({ re: 0, im: bShunt }, Complex.div(Vi, { re: tap * tap, im: 0 }));
    const Iik = Complex.add(Iik_series, Iik_shunt);
    
    const Sik = Complex.mul(Complex.div(Vi, { re: tap, im: 0 }), Complex.conj(Iik));
    
    // Current from k to i
    const Iki_series = Complex.mul(y, Complex.sub(Vk, Complex.div(Vi, { re: tap, im: 0 })));
    const Iki_shunt = Complex.mul({ re: 0, im: bShunt }, Vk);
    const Iki = Complex.add(Iki_series, Iki_shunt);
    
    const Ski = Complex.mul(Vk, Complex.conj(Iki));
    
    const Ploss = (Sik.re + Ski.re) * baseMVA;
    const Qloss = (Sik.im + Ski.im) * baseMVA;
    
    return {
      from: branch.from,
      to: branch.to,
      Pij: Sik.re * baseMVA,
      Qij: Sik.im * baseMVA,
      Pji: Ski.re * baseMVA,
      Qji: Ski.im * baseMVA,
      Ploss,
      Qloss
    };
  });
  
  return {
    converged,
    iterations: iter + 1,
    Vm,
    Va: Va.map(a => a * 180 / Math.PI),
    Pcalc,
    Qcalc,
    lineFlows,
    history,
    Y
  };
}

// Main Component
export default function PowerFlowCalculator() {
  const [busData, setBusData] = useState(defaultBusData);
  const [branchData, setBranchData] = useState(defaultBranchData);
  const [genData, setGenData] = useState(defaultGenData);
  const [results, setResults] = useState(null);
  const [activeTab, setActiveTab] = useState('bus');
  const [showYbus, setShowYbus] = useState(false);
  const [options, setOptions] = useState({ maxIter: 50, tolerance: 1e-6 });

  const runCalculation = useCallback(() => {
    const result = runPowerFlow(busData, branchData, genData, options);
    setResults(result);
  }, [busData, branchData, genData, options]);

  const updateBus = (index, field, value) => {
    const newData = [...busData];
    newData[index] = { ...newData[index], [field]: parseFloat(value) || value };
    setBusData(newData);
  };

  const updateBranch = (index, field, value) => {
    const newData = [...branchData];
    newData[index] = { ...newData[index], [field]: parseFloat(value) };
    setBranchData(newData);
  };

  const updateGen = (index, field, value) => {
    const newData = [...genData];
    newData[index] = { ...newData[index], [field]: parseFloat(value) };
    setGenData(newData);
  };

  const addBus = () => {
    const newBus = {
      bus: busData.length + 1,
      type: BUS_TYPE.PQ,
      Pd: 0,
      Qd: 0,
      Vm: 1.0,
      Va: 0,
      baseKV: 230,
      name: `Bus${busData.length + 1}`
    };
    setBusData([...busData, newBus]);
  };

  const addBranch = () => {
    const newBranch = { from: 1, to: 2, r: 0.01, x: 0.1, b: 0, tap: 1 };
    setBranchData([...branchData, newBranch]);
  };

  const addGen = () => {
    const newGen = { bus: 1, Pg: 0, Qg: 0, Qmax: 300, Qmin: -300, Vg: 1.0, mBase: 100, status: 1 };
    setGenData([...genData, newGen]);
  };

  const totalLosses = useMemo(() => {
    if (!results) return { P: 0, Q: 0 };
    return results.lineFlows.reduce((acc, line) => ({
      P: acc.P + line.Ploss,
      Q: acc.Q + line.Qloss
    }), { P: 0, Q: 0 });
  }, [results]);

  return (
    <div style={styles.container}>
      {/* Header */}
      <header style={styles.header}>
        <div style={styles.headerContent}>
          <div style={styles.logo}>
            <svg width="32" height="32" viewBox="0 0 32 32" fill="none">
              <circle cx="16" cy="16" r="14" stroke="#00ff88" strokeWidth="2"/>
              <circle cx="16" cy="8" r="3" fill="#00ff88"/>
              <circle cx="8" cy="22" r="3" fill="#ff6b35"/>
              <circle cx="24" cy="22" r="3" fill="#4da6ff"/>
              <line x1="16" y1="11" x2="10" y2="19" stroke="#00ff88" strokeWidth="1.5"/>
              <line x1="16" y1="11" x2="22" y2="19" stroke="#4da6ff" strokeWidth="1.5"/>
              <line x1="11" y1="22" x2="21" y2="22" stroke="#ff6b35" strokeWidth="1.5"/>
            </svg>
          </div>
          <div>
            <h1 style={styles.title}>潮流計算ツール</h1>
            <p style={styles.subtitle}>Newton-Raphson Power Flow Calculator</p>
          </div>
        </div>
        <div style={styles.headerRight}>
          <span style={styles.badge}>MATPOWER Style</span>
        </div>
      </header>

      <div style={styles.mainContent}>
        {/* Left Panel - Data Input */}
        <div style={styles.leftPanel}>
          {/* Tabs */}
          <div style={styles.tabs}>
            {['bus', 'branch', 'gen', 'options'].map(tab => (
              <button
                key={tab}
                style={{
                  ...styles.tab,
                  ...(activeTab === tab ? styles.tabActive : {})
                }}
                onClick={() => setActiveTab(tab)}
              >
                {tab === 'bus' ? '母線 (Bus)' :
                 tab === 'branch' ? '送電線 (Branch)' :
                 tab === 'gen' ? '発電機 (Gen)' : '設定'}
              </button>
            ))}
          </div>

          {/* Data Tables */}
          <div style={styles.tableContainer}>
            {activeTab === 'bus' && (
              <>
                <div style={styles.tableHeader}>
                  <h3 style={styles.tableTitle}>Bus Data</h3>
                  <button style={styles.addBtn} onClick={addBus}>+ 追加</button>
                </div>
                <div style={styles.tableWrapper}>
                  <table style={styles.table}>
                    <thead>
                      <tr>
                        <th style={styles.th}>#</th>
                        <th style={styles.th}>Name</th>
                        <th style={styles.th}>Type</th>
                        <th style={styles.th}>Pd (MW)</th>
                        <th style={styles.th}>Qd (MVAr)</th>
                        <th style={styles.th}>Vm (pu)</th>
                        <th style={styles.th}>Va (°)</th>
                      </tr>
                    </thead>
                    <tbody>
                      {busData.map((bus, i) => (
                        <tr key={i} style={styles.tr}>
                          <td style={styles.td}>{bus.bus}</td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              value={bus.name}
                              onChange={e => updateBus(i, 'name', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <select
                              style={styles.select}
                              value={bus.type}
                              onChange={e => updateBus(i, 'type', parseInt(e.target.value))}
                            >
                              <option value={1}>PQ</option>
                              <option value={2}>PV</option>
                              <option value={3}>Slack</option>
                            </select>
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              value={bus.Pd}
                              onChange={e => updateBus(i, 'Pd', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              value={bus.Qd}
                              onChange={e => updateBus(i, 'Qd', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              step="0.01"
                              value={bus.Vm}
                              onChange={e => updateBus(i, 'Vm', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              value={bus.Va}
                              onChange={e => updateBus(i, 'Va', e.target.value)}
                            />
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </>
            )}

            {activeTab === 'branch' && (
              <>
                <div style={styles.tableHeader}>
                  <h3 style={styles.tableTitle}>Branch Data</h3>
                  <button style={styles.addBtn} onClick={addBranch}>+ 追加</button>
                </div>
                <div style={styles.tableWrapper}>
                  <table style={styles.table}>
                    <thead>
                      <tr>
                        <th style={styles.th}>From</th>
                        <th style={styles.th}>To</th>
                        <th style={styles.th}>R (pu)</th>
                        <th style={styles.th}>X (pu)</th>
                        <th style={styles.th}>B (pu)</th>
                        <th style={styles.th}>Tap</th>
                      </tr>
                    </thead>
                    <tbody>
                      {branchData.map((br, i) => (
                        <tr key={i} style={styles.tr}>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              value={br.from}
                              onChange={e => updateBranch(i, 'from', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              value={br.to}
                              onChange={e => updateBranch(i, 'to', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              step="0.001"
                              value={br.r}
                              onChange={e => updateBranch(i, 'r', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              step="0.001"
                              value={br.x}
                              onChange={e => updateBranch(i, 'x', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              step="0.001"
                              value={br.b}
                              onChange={e => updateBranch(i, 'b', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              step="0.01"
                              value={br.tap}
                              onChange={e => updateBranch(i, 'tap', e.target.value)}
                            />
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </>
            )}

            {activeTab === 'gen' && (
              <>
                <div style={styles.tableHeader}>
                  <h3 style={styles.tableTitle}>Generator Data</h3>
                  <button style={styles.addBtn} onClick={addGen}>+ 追加</button>
                </div>
                <div style={styles.tableWrapper}>
                  <table style={styles.table}>
                    <thead>
                      <tr>
                        <th style={styles.th}>Bus</th>
                        <th style={styles.th}>Pg (MW)</th>
                        <th style={styles.th}>Qg (MVAr)</th>
                        <th style={styles.th}>Vg (pu)</th>
                        <th style={styles.th}>Qmax</th>
                        <th style={styles.th}>Qmin</th>
                      </tr>
                    </thead>
                    <tbody>
                      {genData.map((gen, i) => (
                        <tr key={i} style={styles.tr}>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              value={gen.bus}
                              onChange={e => updateGen(i, 'bus', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              value={gen.Pg}
                              onChange={e => updateGen(i, 'Pg', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              value={gen.Qg}
                              onChange={e => updateGen(i, 'Qg', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              step="0.01"
                              value={gen.Vg}
                              onChange={e => updateGen(i, 'Vg', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              value={gen.Qmax}
                              onChange={e => updateGen(i, 'Qmax', e.target.value)}
                            />
                          </td>
                          <td style={styles.td}>
                            <input
                              style={styles.input}
                              type="number"
                              value={gen.Qmin}
                              onChange={e => updateGen(i, 'Qmin', e.target.value)}
                            />
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </>
            )}

            {activeTab === 'options' && (
              <div style={styles.optionsPanel}>
                <h3 style={styles.tableTitle}>計算設定</h3>
                <div style={styles.optionRow}>
                  <label style={styles.optionLabel}>最大反復回数</label>
                  <input
                    style={styles.optionInput}
                    type="number"
                    value={options.maxIter}
                    onChange={e => setOptions({...options, maxIter: parseInt(e.target.value)})}
                  />
                </div>
                <div style={styles.optionRow}>
                  <label style={styles.optionLabel}>収束判定 (tolerance)</label>
                  <input
                    style={styles.optionInput}
                    type="number"
                    step="0.0000001"
                    value={options.tolerance}
                    onChange={e => setOptions({...options, tolerance: parseFloat(e.target.value)})}
                  />
                </div>
                <div style={styles.optionRow}>
                  <label style={styles.optionLabel}>Base MVA</label>
                  <input style={styles.optionInput} type="number" value={100} disabled />
                </div>
              </div>
            )}
          </div>

          {/* Run Button */}
          <button style={styles.runBtn} onClick={runCalculation}>
            <svg width="20" height="20" viewBox="0 0 24 24" fill="none" style={{ marginRight: '8px' }}>
              <path d="M8 5v14l11-7L8 5z" fill="currentColor"/>
            </svg>
            潮流計算実行
          </button>
        </div>

        {/* Right Panel - Results */}
        <div style={styles.rightPanel}>
          {results ? (
            <>
              {/* Status Banner */}
              <div style={{
                ...styles.statusBanner,
                background: results.converged 
                  ? 'linear-gradient(135deg, #0d3320 0%, #1a4d3a 100%)'
                  : 'linear-gradient(135deg, #3d1515 0%, #5c2020 100%)'
              }}>
                <div style={styles.statusIcon}>
                  {results.converged ? '✓' : '✗'}
                </div>
                <div>
                  <div style={styles.statusTitle}>
                    {results.converged ? '収束完了' : '収束失敗'}
                  </div>
                  <div style={styles.statusDetail}>
                    {results.iterations} iterations
                  </div>
                </div>
                <div style={styles.lossInfo}>
                  <div style={styles.lossLabel}>Total Losses</div>
                  <div style={styles.lossValue}>
                    P: {totalLosses.P.toFixed(2)} MW
                  </div>
                  <div style={styles.lossValue}>
                    Q: {totalLosses.Q.toFixed(2)} MVAr
                  </div>
                </div>
              </div>

              {/* Results Tabs */}
              <div style={styles.resultsTabs}>
                <button
                  style={{...styles.resultsTab, ...(showYbus ? {} : styles.resultsTabActive)}}
                  onClick={() => setShowYbus(false)}
                >
                  計算結果
                </button>
                <button
                  style={{...styles.resultsTab, ...(showYbus ? styles.resultsTabActive : {})}}
                  onClick={() => setShowYbus(true)}
                >
                  Ybus行列
                </button>
              </div>

              {!showYbus ? (
                <>
                  {/* Bus Results */}
                  <div style={styles.resultSection}>
                    <h3 style={styles.resultTitle}>
                      <span style={styles.resultIcon}>⚡</span>
                      母線結果 (Bus Results)
                    </h3>
                    <div style={styles.tableWrapper}>
                      <table style={styles.table}>
                        <thead>
                          <tr>
                            <th style={styles.th}>Bus</th>
                            <th style={styles.th}>Name</th>
                            <th style={styles.th}>Vm (pu)</th>
                            <th style={styles.th}>Va (°)</th>
                            <th style={styles.th}>P (MW)</th>
                            <th style={styles.th}>Q (MVAr)</th>
                          </tr>
                        </thead>
                        <tbody>
                          {busData.map((bus, i) => (
                            <tr key={i} style={styles.tr}>
                              <td style={styles.td}>{bus.bus}</td>
                              <td style={{...styles.td, color: '#00ff88'}}>{bus.name}</td>
                              <td style={styles.tdNum}>{results.Vm[i].toFixed(4)}</td>
                              <td style={styles.tdNum}>{results.Va[i].toFixed(2)}</td>
                              <td style={styles.tdNum}>{results.Pcalc[i].toFixed(2)}</td>
                              <td style={styles.tdNum}>{results.Qcalc[i].toFixed(2)}</td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </div>

                  {/* Line Flow Results */}
                  <div style={styles.resultSection}>
                    <h3 style={styles.resultTitle}>
                      <span style={styles.resultIcon}>⟷</span>
                      線路潮流 (Line Flows)
                    </h3>
                    <div style={styles.tableWrapper}>
                      <table style={styles.table}>
                        <thead>
                          <tr>
                            <th style={styles.th}>From</th>
                            <th style={styles.th}>To</th>
                            <th style={styles.th}>P_ij (MW)</th>
                            <th style={styles.th}>Q_ij (MVAr)</th>
                            <th style={styles.th}>P_ji (MW)</th>
                            <th style={styles.th}>Q_ji (MVAr)</th>
                            <th style={styles.th}>Loss (MW)</th>
                          </tr>
                        </thead>
                        <tbody>
                          {results.lineFlows.map((line, i) => (
                            <tr key={i} style={styles.tr}>
                              <td style={styles.td}>{line.from}</td>
                              <td style={styles.td}>{line.to}</td>
                              <td style={{...styles.tdNum, color: line.Pij >= 0 ? '#4da6ff' : '#ff6b35'}}>
                                {line.Pij.toFixed(2)}
                              </td>
                              <td style={styles.tdNum}>{line.Qij.toFixed(2)}</td>
                              <td style={{...styles.tdNum, color: line.Pji >= 0 ? '#4da6ff' : '#ff6b35'}}>
                                {line.Pji.toFixed(2)}
                              </td>
                              <td style={styles.tdNum}>{line.Qji.toFixed(2)}</td>
                              <td style={{...styles.tdNum, color: '#ff6b35'}}>
                                {line.Ploss.toFixed(3)}
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </div>

                  {/* Convergence History */}
                  <div style={styles.resultSection}>
                    <h3 style={styles.resultTitle}>
                      <span style={styles.resultIcon}>📈</span>
                      収束履歴
                    </h3>
                    <div style={styles.convergenceChart}>
                      {results.history.map((h, i) => (
                        <div key={i} style={styles.convergenceBar}>
                          <div
                            style={{
                              ...styles.convergenceBarFill,
                              height: `${Math.min(100, Math.max(5, -Math.log10(h.maxMismatch) * 10))}%`
                            }}
                          />
                          <span style={styles.convergenceLabel}>{i + 1}</span>
                        </div>
                      ))}
                    </div>
                    <div style={styles.convergenceInfo}>
                      Final mismatch: {results.history[results.history.length - 1]?.maxMismatch.toExponential(4)}
                    </div>
                  </div>
                </>
              ) : (
                /* Ybus Matrix */
                <div style={styles.resultSection}>
                  <h3 style={styles.resultTitle}>
                    <span style={styles.resultIcon}>🔢</span>
                    アドミタンス行列 (Ybus)
                  </h3>
                  <div style={styles.matrixWrapper}>
                    <table style={styles.matrixTable}>
                      <thead>
                        <tr>
                          <th style={styles.matrixTh}></th>
                          {busData.map((_, j) => (
                            <th key={j} style={styles.matrixTh}>{j + 1}</th>
                          ))}
                        </tr>
                      </thead>
                      <tbody>
                        {results.Y.map((row, i) => (
                          <tr key={i}>
                            <td style={styles.matrixTh}>{i + 1}</td>
                            {row.map((cell, j) => (
                              <td key={j} style={{
                                ...styles.matrixTd,
                                background: i === j ? 'rgba(0,255,136,0.1)' : 'transparent'
                              }}>
                                <div style={styles.complexNum}>
                                  {cell.re.toFixed(3)}
                                </div>
                                <div style={{
                                  ...styles.complexNum,
                                  color: cell.im >= 0 ? '#4da6ff' : '#ff6b35'
                                }}>
                                  {cell.im >= 0 ? '+' : ''}{cell.im.toFixed(3)}j
                                </div>
                              </td>
                            ))}
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
              )}
            </>
          ) : (
            <div style={styles.placeholder}>
              <svg width="80" height="80" viewBox="0 0 80 80" fill="none" style={{ opacity: 0.3 }}>
                <circle cx="40" cy="40" r="35" stroke="#00ff88" strokeWidth="2" strokeDasharray="8 4"/>
                <circle cx="40" cy="20" r="6" fill="#00ff88"/>
                <circle cx="20" cy="55" r="6" fill="#ff6b35"/>
                <circle cx="60" cy="55" r="6" fill="#4da6ff"/>
                <line x1="40" y1="26" x2="24" y2="49" stroke="#888" strokeWidth="2"/>
                <line x1="40" y1="26" x2="56" y2="49" stroke="#888" strokeWidth="2"/>
                <line x1="26" y1="55" x2="54" y2="55" stroke="#888" strokeWidth="2"/>
              </svg>
              <p style={styles.placeholderText}>
                左側でデータを入力し<br/>
                「潮流計算実行」をクリック
              </p>
            </div>
          )}
        </div>
      </div>

      {/* Footer */}
      <footer style={styles.footer}>
        <div>Newton-Raphson Power Flow | IEEE 9-Bus Test System</div>
        <div style={{ opacity: 0.6 }}>Base MVA: 100 | Per-unit system</div>
      </footer>
    </div>
  );
}

// Styles
const styles = {
  container: {
    minHeight: '100vh',
    background: 'linear-gradient(180deg, #0a0f14 0%, #111820 50%, #0d1218 100%)',
    color: '#e0e0e0',
    fontFamily: "'JetBrains Mono', 'SF Mono', 'Fira Code', monospace",
    display: 'flex',
    flexDirection: 'column',
  },
  header: {
    padding: '16px 24px',
    background: 'linear-gradient(90deg, rgba(0,255,136,0.08) 0%, transparent 50%, rgba(77,166,255,0.08) 100%)',
    borderBottom: '1px solid rgba(0,255,136,0.2)',
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
  },
  headerContent: {
    display: 'flex',
    alignItems: 'center',
    gap: '16px',
  },
  logo: {
    width: '48px',
    height: '48px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
  },
  title: {
    margin: 0,
    fontSize: '24px',
    fontWeight: 700,
    background: 'linear-gradient(135deg, #00ff88 0%, #4da6ff 100%)',
    WebkitBackgroundClip: 'text',
    WebkitTextFillColor: 'transparent',
  },
  subtitle: {
    margin: 0,
    fontSize: '12px',
    color: '#888',
    letterSpacing: '1px',
  },
  headerRight: {
    display: 'flex',
    gap: '12px',
  },
  badge: {
    padding: '6px 12px',
    background: 'rgba(0,255,136,0.15)',
    border: '1px solid rgba(0,255,136,0.3)',
    borderRadius: '4px',
    fontSize: '11px',
    color: '#00ff88',
    letterSpacing: '1px',
  },
  mainContent: {
    flex: 1,
    display: 'flex',
    gap: '1px',
    background: 'rgba(0,255,136,0.1)',
    padding: '1px',
  },
  leftPanel: {
    width: '50%',
    background: '#0d1218',
    display: 'flex',
    flexDirection: 'column',
    padding: '16px',
  },
  rightPanel: {
    width: '50%',
    background: '#0d1218',
    display: 'flex',
    flexDirection: 'column',
    padding: '16px',
    overflowY: 'auto',
  },
  tabs: {
    display: 'flex',
    gap: '4px',
    marginBottom: '16px',
  },
  tab: {
    flex: 1,
    padding: '10px 8px',
    background: 'rgba(255,255,255,0.03)',
    border: '1px solid rgba(255,255,255,0.08)',
    borderRadius: '6px',
    color: '#888',
    fontSize: '12px',
    cursor: 'pointer',
    transition: 'all 0.2s',
  },
  tabActive: {
    background: 'rgba(0,255,136,0.15)',
    borderColor: 'rgba(0,255,136,0.4)',
    color: '#00ff88',
  },
  tableContainer: {
    flex: 1,
    overflowY: 'auto',
  },
  tableHeader: {
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    marginBottom: '12px',
  },
  tableTitle: {
    margin: 0,
    fontSize: '14px',
    color: '#00ff88',
    fontWeight: 600,
  },
  addBtn: {
    padding: '6px 12px',
    background: 'rgba(77,166,255,0.15)',
    border: '1px solid rgba(77,166,255,0.3)',
    borderRadius: '4px',
    color: '#4da6ff',
    fontSize: '12px',
    cursor: 'pointer',
  },
  tableWrapper: {
    overflowX: 'auto',
  },
  table: {
    width: '100%',
    borderCollapse: 'collapse',
    fontSize: '12px',
  },
  th: {
    padding: '10px 8px',
    background: 'rgba(0,0,0,0.3)',
    borderBottom: '1px solid rgba(255,255,255,0.1)',
    textAlign: 'left',
    color: '#888',
    fontWeight: 500,
    whiteSpace: 'nowrap',
  },
  tr: {
    borderBottom: '1px solid rgba(255,255,255,0.05)',
  },
  td: {
    padding: '6px 8px',
  },
  tdNum: {
    padding: '6px 8px',
    textAlign: 'right',
    fontFamily: "'JetBrains Mono', monospace",
  },
  input: {
    width: '100%',
    padding: '6px 8px',
    background: 'rgba(255,255,255,0.05)',
    border: '1px solid rgba(255,255,255,0.1)',
    borderRadius: '4px',
    color: '#e0e0e0',
    fontSize: '12px',
    fontFamily: "'JetBrains Mono', monospace",
  },
  select: {
    width: '100%',
    padding: '6px 8px',
    background: 'rgba(255,255,255,0.05)',
    border: '1px solid rgba(255,255,255,0.1)',
    borderRadius: '4px',
    color: '#e0e0e0',
    fontSize: '12px',
  },
  optionsPanel: {
    padding: '16px',
  },
  optionRow: {
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    marginBottom: '16px',
    padding: '12px',
    background: 'rgba(255,255,255,0.03)',
    borderRadius: '6px',
  },
  optionLabel: {
    fontSize: '13px',
    color: '#aaa',
  },
  optionInput: {
    width: '120px',
    padding: '8px 12px',
    background: 'rgba(255,255,255,0.05)',
    border: '1px solid rgba(255,255,255,0.1)',
    borderRadius: '4px',
    color: '#e0e0e0',
    fontSize: '13px',
    fontFamily: "'JetBrains Mono', monospace",
    textAlign: 'right',
  },
  runBtn: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    marginTop: '16px',
    padding: '14px 24px',
    background: 'linear-gradient(135deg, #00cc6a 0%, #00ff88 100%)',
    border: 'none',
    borderRadius: '8px',
    color: '#0a0f14',
    fontSize: '15px',
    fontWeight: 700,
    cursor: 'pointer',
    transition: 'transform 0.2s, box-shadow 0.2s',
    boxShadow: '0 4px 20px rgba(0,255,136,0.3)',
  },
  statusBanner: {
    display: 'flex',
    alignItems: 'center',
    gap: '16px',
    padding: '16px 20px',
    borderRadius: '8px',
    marginBottom: '16px',
  },
  statusIcon: {
    width: '40px',
    height: '40px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    background: 'rgba(0,0,0,0.3)',
    borderRadius: '50%',
    fontSize: '20px',
  },
  statusTitle: {
    fontSize: '16px',
    fontWeight: 700,
    color: '#fff',
  },
  statusDetail: {
    fontSize: '12px',
    color: 'rgba(255,255,255,0.6)',
  },
  lossInfo: {
    marginLeft: 'auto',
    textAlign: 'right',
  },
  lossLabel: {
    fontSize: '11px',
    color: 'rgba(255,255,255,0.5)',
    marginBottom: '4px',
  },
  lossValue: {
    fontSize: '12px',
    color: '#ff6b35',
    fontFamily: "'JetBrains Mono', monospace",
  },
  resultsTabs: {
    display: 'flex',
    gap: '8px',
    marginBottom: '16px',
  },
  resultsTab: {
    padding: '8px 16px',
    background: 'transparent',
    border: '1px solid rgba(255,255,255,0.1)',
    borderRadius: '4px',
    color: '#888',
    fontSize: '12px',
    cursor: 'pointer',
  },
  resultsTabActive: {
    background: 'rgba(77,166,255,0.15)',
    borderColor: 'rgba(77,166,255,0.4)',
    color: '#4da6ff',
  },
  resultSection: {
    marginBottom: '20px',
  },
  resultTitle: {
    display: 'flex',
    alignItems: 'center',
    gap: '8px',
    margin: '0 0 12px 0',
    fontSize: '14px',
    color: '#00ff88',
    fontWeight: 600,
  },
  resultIcon: {
    fontSize: '16px',
  },
  convergenceChart: {
    display: 'flex',
    alignItems: 'flex-end',
    gap: '4px',
    height: '80px',
    padding: '12px',
    background: 'rgba(0,0,0,0.2)',
    borderRadius: '6px',
  },
  convergenceBar: {
    flex: 1,
    height: '100%',
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
    justifyContent: 'flex-end',
  },
  convergenceBarFill: {
    width: '100%',
    background: 'linear-gradient(180deg, #00ff88 0%, #00cc6a 100%)',
    borderRadius: '2px 2px 0 0',
    transition: 'height 0.3s',
  },
  convergenceLabel: {
    fontSize: '9px',
    color: '#666',
    marginTop: '4px',
  },
  convergenceInfo: {
    marginTop: '8px',
    fontSize: '11px',
    color: '#888',
    textAlign: 'center',
  },
  matrixWrapper: {
    overflowX: 'auto',
    background: 'rgba(0,0,0,0.2)',
    borderRadius: '6px',
    padding: '12px',
  },
  matrixTable: {
    borderCollapse: 'collapse',
    fontSize: '10px',
  },
  matrixTh: {
    padding: '8px',
    background: 'rgba(0,255,136,0.1)',
    color: '#00ff88',
    fontWeight: 600,
    textAlign: 'center',
    minWidth: '60px',
  },
  matrixTd: {
    padding: '6px',
    textAlign: 'center',
    border: '1px solid rgba(255,255,255,0.05)',
  },
  complexNum: {
    fontSize: '10px',
    lineHeight: '1.4',
  },
  placeholder: {
    flex: 1,
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
    justifyContent: 'center',
    gap: '24px',
  },
  placeholderText: {
    textAlign: 'center',
    color: '#555',
    fontSize: '14px',
    lineHeight: '1.8',
  },
  footer: {
    padding: '12px 24px',
    borderTop: '1px solid rgba(255,255,255,0.05)',
    display: 'flex',
    justifyContent: 'space-between',
    fontSize: '11px',
    color: '#555',
  },
};
