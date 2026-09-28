import Lake
open Lake DSL

package «pauli-fold» where
  version := v!"0.1.0"

require «tzap-lean» from ".." / "lean"

@[default_target]
lean_lib PauliFold where
  roots := #[`PauliFold]
