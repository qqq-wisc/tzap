//! Anchor-frontier computation for incremental runs over an evolving circuit.
//!
//! A pass instance that has already analyzed a previous version of its input
//! only needs to anchor windows near gates that changed since then. A window
//! whose gates are all unchanged was analyzed against the same synthesis
//! MURM in an earlier run and selected nothing — had it selected a rewrite,
//! its gates would have been replaced and would differ now — so skipping it
//! leaves the output identical to a full sweep.

use crate::circuit::Circuit;

use super::unique_qubits;

/// Gate indices per qubit; a multi-qubit gate appears in each operand's stream.
fn qubit_streams(circuit: &Circuit) -> Vec<Vec<usize>> {
    let mut streams = vec![Vec::new(); circuit.num_qubits];
    for (index, gate) in circuit.gates.iter().enumerate() {
        for qubit in unique_qubits(gate) {
            streams[qubit as usize].push(index);
        }
    }
    streams
}

/// Compute the anchor frontier for `circuit` given the instance's previous
/// input: a bitmap marking every gate allowed to anchor a window. `None`
/// anchors everywhere (no previous input, or qubit-count mismatch).
pub(super) fn anchor_frontier(
    circuit: &Circuit,
    prev: Option<&Circuit>,
    window_gates: usize,
) -> Option<Vec<bool>> {
    let prev = prev?;
    if prev.num_qubits != circuit.num_qubits {
        return None;
    }

    let streams = qubit_streams(circuit);
    let prev_streams = qubit_streams(prev);

    // Seed: per qubit, trim the streams' common prefix and suffix and mark
    // the remainder dirty, widened by one surviving gate on each side so a
    // pure insertion or deletion (empty remainder) still seeds the survivors
    // it made adjacent.
    let mut dirty = vec![false; circuit.gates.len()];
    let mut ring = Vec::new();
    for (stream, prev_stream) in streams.iter().zip(&prev_streams) {
        let same = |c: usize, p: usize| circuit.gates[stream[c]] == prev.gates[prev_stream[p]];
        let limit = stream.len().min(prev_stream.len());
        let mut prefix = 0;
        while prefix < limit && same(prefix, prefix) {
            prefix += 1;
        }
        if prefix == stream.len() && prefix == prev_stream.len() {
            continue;
        }
        let mut suffix = 0;
        while suffix < limit - prefix
            && same(stream.len() - 1 - suffix, prev_stream.len() - 1 - suffix)
        {
            suffix += 1;
        }
        let lo = prefix.saturating_sub(1);
        let hi = (stream.len() - suffix + 1).min(stream.len());
        for &gate_index in &stream[lo..hi] {
            if !dirty[gate_index] {
                dirty[gate_index] = true;
                ring.push(gate_index);
            }
        }
    }

    // Dilate: a window holds at most `window_gates` connected gates, and
    // window gates consecutive on a shared qubit are also consecutive in that
    // qubit's stream, so every member of a window containing a dirty gate —
    // its anchor in particular — is within `window_gates - 1` stream steps of
    // one. Marking that ball makes anchoring only inside it exhaustive.
    let mut frontier = dirty;
    for _ in 1..window_gates {
        if ring.is_empty() {
            break;
        }
        let mut next = Vec::new();
        for &gate_index in &ring {
            for qubit in unique_qubits(&circuit.gates[gate_index]) {
                let stream = &streams[qubit as usize];
                let position = stream
                    .binary_search(&gate_index)
                    .expect("gate is in each of its qubits' streams");
                let mut mark = |neighbor: usize| {
                    if !frontier[neighbor] {
                        frontier[neighbor] = true;
                        next.push(neighbor);
                    }
                };
                if position > 0 {
                    mark(stream[position - 1]);
                }
                if position + 1 < stream.len() {
                    mark(stream[position + 1]);
                }
            }
        }
        ring = next;
    }
    Some(frontier)
}
