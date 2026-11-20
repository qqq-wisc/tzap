import Lake
open Lake DSL

package «pauli-fold» where
  version := v!"0.1.0"

-- The concrete semantics: tzap's Lean development, in the `lean` directory of the tzap
-- repository, pinned to a commit.
require «tzap-lean» from git
  "https://github.com/qqq-wisc/tzap.git" @ "5b353e78c17c3513733299e8308d6d45e54c7496" / "lean"

@[default_target]
lean_lib PauliFold where
  roots := #[`PauliFold]
