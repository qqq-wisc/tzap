//! `FuseDecomposeRz`: a composite pass that fuses Clifford+RZ phases FIRST,
//! then dispatches the surviving RZ gates to gridsynth.
//!
//! This pass adds no new logic — it only re-orders calls to existing passes
//! already in the crate. The motivation: tzap's default `--decompose-rz`
//! pipeline runs `DecomposeRz` before `PhaseFoldRand`, which expands every
//! RZ to ~3·log(1/ε) T gates before any fusion can occur. Once expanded,
//! the H gates internal to each gridsynth sequence break parity, so the
//! arbitrary-angle merge in `PhaseFoldRand::record_float` can no longer
//! recover the would-have-been fusion. Running the fold first sees the
//! raw RZ gates and merges them at the symbolic level; the dispatched
//! gridsynth call set is then strictly smaller.
//!
//! Pipeline:
//!   1. CancelPairs           (existing — eliminates self-inverse pairs)
//!   2. PhaseFoldRand         (existing — symbolic RZ fusion via record_float)
//!   3. DecomposeRz           (existing — gridsynth per surviving RZ)
//!   4. CancelPairs           (existing — clean up the T-sequence boundary)
//!   5. PhaseFoldRand         (existing — fold integer-π/4 multiples)

use indicatif::ProgressBar;

use crate::cancel::CancelPairs;
use crate::circuit::Circuit;
use crate::decompose_rz::DecomposeRz;
use crate::pass::Pass;
use crate::phase_fold_rand::PhaseFoldRand;

pub struct FuseDecomposeRz {
    pub epsilon: f64,
}

impl Default for FuseDecomposeRz {
    fn default() -> Self {
        Self { epsilon: 1e-10 }
    }
}

impl Pass for FuseDecomposeRz {
    fn name(&self) -> &str {
        "Fuse-first Rz → Clifford+T (PhaseFold → DecomposeRz → PhaseFold)"
    }

    fn run_with_progress(&self, circuit: &Circuit, _pb: &ProgressBar) -> Circuit {
        // Stage 1: cancel + fuse RZ symbolically while gates are still RZ.
        let c = CancelPairs.run(circuit);
        let c = PhaseFoldRand.run(&c);

        // Stage 2: now decompose only the surviving (fused) RZ angles.
        let c = DecomposeRz { epsilon: self.epsilon }.run(&c);

        // Stage 3: clean up the gridsynth output (integer-π/4 fold + cancel).
        let c = CancelPairs.run(&c);
        let c = PhaseFoldRand.run(&c);

        c
    }
}
