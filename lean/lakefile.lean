import Lake
open Lake DSL

package «tzap-lean» where
  version := v!"0.1.0"

require "leanprover-community" / "mathlib" @ git "v4.33.1"

@[default_target]
lean_lib TzapLean where
  roots := #[`TzapLean]

@[default_target] lean_exe «tzap-lean» where
  root := `Main
