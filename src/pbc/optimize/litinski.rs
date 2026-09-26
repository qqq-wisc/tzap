//! Litinski's greedy layering ("A Game of Surface Codes", Sec. 1), as a
//! baseline, and T-depth measurement.

use super::*;

impl Optimizer<'_> {
    /// Litinski's algorithm, per segment: partition the rotations into layers
    /// of mutually commuting rotations (a new layer whenever the next rotation
    /// anticommutes with the current one); then repeatedly move each rotation
    /// of layer i + 1 into layer i if it commutes with all of layer i, until
    /// nothing moves. Equal rotations that meet in a layer combine; a combined
    /// Clifford is commuted to the end (into the frame), which changes later
    /// axes, so layering restarts until no more rotations combine. Returns
    /// whether the frame changed.
    pub(super) fn litinski(
        &mut self,
        items: &mut Vec<Item>,
        frame: &mut [(u32, i8)],
    ) -> Result<bool, PbcError> {
        // Combined Cliffords always go to the end, as in the paper.
        self.options.clifford_to_frame = true;
        let mut frame_changed = false;
        // Initial Cliffords, if any, go to the end first.
        frame_changed |= self.stream(items, frame, false)? > 0;
        loop {
            let mut combined = 0;
            let mut output = Vec::with_capacity(items.len());
            let mut segment = Vec::new();
            for &item in items.iter() {
                match item {
                    Item::Rot(rot) => segment.push(rot),
                    Item::Barrier { .. } => {
                        combined += self.layer_segment(&mut segment);
                        output.extend(segment.drain(..).map(Item::Rot));
                        output.push(item);
                    }
                }
            }
            combined += self.layer_segment(&mut segment);
            output.extend(segment.drain(..).map(Item::Rot));
            *items = output;
            self.stats.merges += combined;
            if combined == 0 {
                break;
            }
            frame_changed |= self.stream(items, frame, false)? > 0;
        }
        Ok(frame_changed)
    }

    /// Layer one segment and write it back in layer order. Returns the number
    /// of rotations combined.
    fn layer_segment(&mut self, rots: &mut Vec<Rot>) -> usize {
        // Naive partition; equal rotations in one layer combine.
        let mut combined = 0;
        let mut layers: Vec<Vec<Rot>> = Vec::new();
        for &rot in rots.iter() {
            match layers.last_mut() {
                Some(layer) if self.commutes_with_all(layer, rot) => {
                    match layer.iter_mut().find(|r| r.axis == rot.axis) {
                        Some(equal) => {
                            equal.k = normalize(i32::from(equal.k) + i32::from(rot.k));
                            combined += 1;
                        }
                        None => layer.push(rot),
                    }
                }
                _ => layers.push(vec![rot]),
            }
        }
        loop {
            let mut moved = false;
            for i in 0..layers.len().saturating_sub(1) {
                let next = std::mem::take(&mut layers[i + 1]);
                let mut stay = Vec::with_capacity(next.len());
                for rot in next {
                    if !self.commutes_with_all(&layers[i], rot) {
                        stay.push(rot);
                        continue;
                    }
                    moved = true;
                    match layers[i].iter_mut().find(|r| r.axis == rot.axis) {
                        Some(equal) => {
                            equal.k = normalize(i32::from(equal.k) + i32::from(rot.k));
                            combined += 1;
                        }
                        None => layers[i].push(rot),
                    }
                }
                layers[i + 1] = stay;
            }
            for layer in &mut layers {
                layer.retain(|r| r.k != 0);
            }
            layers.retain(|layer| !layer.is_empty());
            if !moved {
                break;
            }
        }
        *rots = layers.concat();
        combined
    }

    fn commutes_with_all(&self, layer: &[Rot], rot: Rot) -> bool {
        layer
            .iter()
            .all(|r| r.support & rot.support == 0 || !self.axes.anticommute(r.axis, rot.axis))
    }

    /// T depth: odd-angle rotations placed as early as possible into layers
    /// of mutually commuting rotations (one past the latest layer holding an
    /// anticommuting rotation), summed over segments. Even-angle rotations
    /// are Clifford and ignored.
    pub(super) fn t_depth(&self, items: &[Item]) -> usize {
        let mut total = 0;
        // Per layer: its rotations and the union of their signatures.
        let mut layers: Vec<(u64, Vec<Rot>)> = Vec::new();
        // Per signature bit: the highest layer touching it, plus one.
        let mut top = [0usize; 64];
        for item in items {
            match *item {
                Item::Rot(rot) if rot.k % 2 != 0 => {
                    // No layer above the highest one sharing a bit can hold
                    // an anticommuting rotation.
                    let bound = bits(rot.support).map(|b| top[b]).max().unwrap_or(0);
                    let mut place = 0;
                    for index in (0..bound).rev() {
                        let (signature, layer) = &layers[index];
                        if signature & rot.support != 0 && !self.commutes_with_all(layer, rot) {
                            place = index + 1;
                            break;
                        }
                    }
                    if place == layers.len() {
                        layers.push((0, Vec::new()));
                    }
                    layers[place].0 |= rot.support;
                    layers[place].1.push(rot);
                    for b in bits(rot.support) {
                        top[b] = top[b].max(place + 1);
                    }
                }
                Item::Rot(_) => {}
                Item::Barrier { .. } => {
                    total += layers.len();
                    layers.clear();
                    top = [0; 64];
                }
            }
        }
        total + layers.len()
    }
}
