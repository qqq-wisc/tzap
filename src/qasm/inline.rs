use super::scan::nninteger;
use super::{
    Evaluation, Operands,
    definitions::{self, Callable, Registry},
};
use crate::circuit::Gate;
use smallvec::SmallVec;

#[derive(Clone, Copy)]
pub(super) enum Operand {
    Indexed(u32),
    Register { offset: usize, width: usize },
}
impl Operand {
    fn at(self, lane: usize) -> u32 {
        match self {
            Self::Indexed(q) => q,
            Self::Register { offset, .. } => (offset + lane) as u32,
        }
    }
}

/// Check for aliasing in any broadcast lane without iterating over the lanes.
/// Register/register equality depends on their starting offsets; an indexed
/// qubit aliases a register exactly when it lies in that register's range.
fn distinct_operands(args: &[Operand], name: &str, line: usize) -> Result<(), String> {
    let mut indexed = SmallVec::<[u32; 3]>::new();
    let mut registers = SmallVec::<[(usize, usize); 3]>::new();
    for arg in args {
        match *arg {
            Operand::Indexed(q) => indexed.push(q),
            Operand::Register { offset, width } => registers.push((offset, width)),
        }
    }
    indexed.sort_unstable();
    registers.sort_unstable();
    let duplicates = indexed.windows(2).any(|pair| pair[0] == pair[1])
        || registers.windows(2).any(|pair| pair[0].0 == pair[1].0);
    let intersection = indexed.iter().any(|&q| {
        let q = q as usize;
        let next = registers.partition_point(|&(offset, _)| offset <= q);
        next > 0 && q - registers[next - 1].0 < registers[next - 1].1
    });
    if duplicates || intersection {
        return Err(format!(
            "line {line}: {name} requires distinct qubit operands"
        ));
    }
    Ok(())
}

pub(super) fn operand(
    s: &str,
    registers: &[(String, usize, usize)],
    line: usize,
    classical: bool,
) -> Result<Operand, String> {
    // Most large inputs have one register and canonical indexed operands.
    // Resolve those without searching twice for brackets or scanning a name.
    if let [(name, offset, width)] = registers
        && let Some(rest) = s.strip_prefix(name.as_str())
        && let Some(index) = rest.strip_prefix('[').and_then(|r| r.strip_suffix(']'))
        && let Some(index) = nninteger(index)
        && index < *width
    {
        return Ok(Operand::Indexed((offset + index) as u32));
    }
    let kind = if classical { "cbit" } else { "qubit" };
    let register_kind = if classical {
        "classical register"
    } else {
        "register"
    };
    let s = s.trim();
    let (name, index) = if let Some(open) = s.find('[') {
        let close = s[open + 1..]
            .find(']')
            .map(|i| open + 1 + i)
            .ok_or_else(|| format!("line {line}: invalid {kind} operand '{s}'"))?;
        if !s[close + 1..].trim().is_empty() {
            return Err(format!("line {line}: invalid {kind} operand '{s}'"));
        }
        let text = s[open + 1..close].trim();
        let index =
            nninteger(text).ok_or_else(|| format!("line {line}: bad {kind} index '{text}'"))?;
        (s[..open].trim(), Some(index))
    } else {
        (s, None)
    };
    let (_, offset, width) = registers
        .iter()
        .find(|(n, _, _)| n == name)
        .ok_or_else(|| format!("line {line}: unknown {register_kind} '{name}'"))?;
    if let Some(index) = index {
        if index >= *width {
            return Err(format!(
                "line {line}: index {index} out of range for {register_kind} '{name}' (size {width})"
            ));
        }
        Ok(Operand::Indexed((offset + index) as u32))
    } else {
        Ok(Operand::Register {
            offset: *offset,
            width: *width,
        })
    }
}

#[derive(Default)]
pub(super) struct Budget {
    gates: usize,
    work: usize,
}

/// Avoid expression/descriptor machinery for ordinary indexed native gates.
/// A name the program defined itself takes the general path.
pub(super) fn native_indexed(
    name: &str,
    rest: &str,
    line: usize,
    registry: &definitions::Registry,
    registers: &[(String, usize, usize)],
    gates: &mut Vec<Gate>,
) -> Result<bool, String> {
    if registry.contains(name) {
        return Ok(false);
    }
    let Some(native) = definitions::Native::lookup(name) else {
        return Ok(false);
    };
    let signature = native.signature();
    if signature.0 != 0 || rest.starts_with('(') {
        return Ok(false);
    }
    let mut qubits = Operands::new();
    for s in definitions::comma_list(rest) {
        match operand(s, registers, line, false)? {
            Operand::Indexed(q) => qubits.push(q),
            Operand::Register { .. } => return Ok(false),
        }
    }
    definitions::check_signature(name, 0, qubits.len(), signature, line)?;
    definitions::distinct(&qubits, name, line)?;
    native.emit(&[], &qubits, gates);
    Ok(true)
}

pub(super) fn call(
    source: &str,
    line: usize,
    registry: &Registry,
    registers: &[(String, usize, usize)],
    evaluation: Evaluation<'_>,
    budget: &mut Budget,
    gates: &mut Vec<Gate>,
) -> Result<(), String> {
    let Evaluation { limits, angles } = evaluation;
    let (name, params, qubits) = definitions::call_parts(source, line)?;
    let callee = registry.lookup(name, line)?;
    let expressions = definitions::expressions(params, None, limits, line, angles)?;
    let parameters = definitions::evaluate(&expressions, &[], line, angles)?;
    let args: SmallVec<[Operand; 3]> = definitions::comma_list(qubits)
        .map(|s| operand(s, registers, line, false))
        .collect::<Result<_, _>>()?;
    definitions::check_signature(
        name,
        parameters.len(),
        args.len(),
        registry.signature(callee),
        line,
    )?;
    let mut width = None;
    for arg in &args {
        if let Operand::Register { width: n, .. } = arg {
            if width.is_some_and(|w| w != *n) {
                return Err(format!(
                    "line {line}: {name} register operand size mismatch"
                ));
            }
            width = Some(*n);
        }
    }
    let lanes = width.unwrap_or(1);
    // Both custom calls and native register broadcasts expand compact source.
    // Explicit indexed native calls remain outside the expansion budget.
    let charged = if matches!(callee, Callable::Custom(_)) || width.is_some() {
        let cost = registry.cost(callee);
        if cost.gates == usize::MAX || cost.work == usize::MAX {
            return Err(format!("line {line}: '{name}' expansion cost overflow"));
        }
        if cost.depth > limits.call_depth {
            return Err(format!("line {line}: '{name}' call depth limit exceeded"));
        }
        let count = cost
            .gates
            .checked_mul(lanes)
            .and_then(|n| n.checked_add(budget.gates))
            .filter(|n| *n <= limits.expanded_gates)
            .ok_or_else(|| format!("line {line}: '{name}' expanded gate limit exceeded"))?;
        let work = cost
            .work
            .checked_mul(lanes)
            .and_then(|n| n.checked_add(budget.work))
            .filter(|n| *n <= limits.expansion_work)
            .ok_or_else(|| format!("line {line}: '{name}' expansion work limit exceeded"))?;
        Some((count, work))
    } else {
        None
    };
    // No instances means no alias constraint, retaining zero-width behavior.
    if lanes != 0 {
        distinct_operands(&args, name, line)?;
    }
    if let Some((count, work)) = charged {
        gates
            .try_reserve(count - budget.gates)
            .map_err(|e| format!("line {line}: cannot allocate gate expansion: {e}"))?;
        budget.gates = count;
        budget.work = work;
    }
    for lane in 0..lanes {
        let qubits: Operands = args.iter().map(|a| a.at(lane)).collect();
        match callee {
            Callable::Native(n) => n.emit(&parameters, &qubits, gates),
            Callable::Custom(id) => definitions::expand(
                registry,
                id,
                parameters.clone(),
                qubits,
                line,
                gates,
                angles,
            )?,
        }
    }
    Ok(())
}
