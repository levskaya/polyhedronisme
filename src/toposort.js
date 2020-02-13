// Polyhédronisme
// ================================================================================================
//
// A toy for constructing and manipulating polyhedra.
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License
//

/* eslint-disable camelcase */
/* eslint-disable standard/array-bracket-even-spacing */
/* eslint-disable no-multi-spaces */

import { __range__ } from './geo';

// Single-shot class for performing Kahn's algorithm to topologically sort nodes in a DAG.
export class TopoSorter {
  // constructor of initially empty directed graph
  constructor (num_nodes) {
    // graph
    this.fwdgraph = __range__(0, num_nodes, false).map(i => []);
    // inverted graph
    this.revgraph = __range__(0, num_nodes, false).map(i => []);
    this.num_nodes = num_nodes;
    this.num_edges = 0;
  }

  addEdge (n, m) {
    this.fwdgraph[n].push(m);
    this.revgraph[m].push(n);
    this.num_edges++;
  }

  removeEdge (n, m) {
    const idxm = this.fwdgraph[n].indexOf(m);
    const idxn = this.revgraph[m].indexOf(n);
    if (idxm > -1 && idxn > -1) {
      this.fwdgraph[n].splice(idxm, 1);
      this.revgraph[m].splice(idxn, 1);
      this.num_edges--;
    } else {
      throw new Error('Malformed graph.');
    }
  }

  // Get topologically sorted nodes, destructively!
  sort () {
    let sorted_list = [];
    let unconstrained_set = [];
    // Initialize unconstrained_set to nodes without inbound edges.
    for (let i = 0; i < this.num_nodes; i++) {
      if (this.revgraph[i].length === 0) {
        unconstrained_set.push(i);
      }
    }
    // Iterative topological sort.
    while (unconstrained_set.length > 0) {
      let n = unconstrained_set.shift();
      sorted_list.push(n);
      for (let m of this.fwdgraph[n].slice()) { // must work with a copy of data
        this.removeEdge(n, m);
        if (this.revgraph[m].length === 0) {
          unconstrained_set.push(m);
        }
      }
    }
    // Check for cycles
    for (let i = 0; i < this.num_nodes; i++) {
      if (this.fwdgraph[i].length !== 0 || this.revgraph[i].length !== 0) {
        throw new Error('Topological sort failed, a cycle was detected' +
                        ' due to self-intersecting geometry.');
      }
    }
    return sorted_list;
  }
};
