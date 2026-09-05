import TzapLean.Cli

/-- The `tzap-lean` executable. -/
def main (argv : List String) : IO UInt32 := TzapLean.Cli.main argv
