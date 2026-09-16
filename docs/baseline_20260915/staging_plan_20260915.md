# Final explicit staging plan

Reviewed 2026-09-16. Supersedes earlier 293/294-file recommendations. No staging or commit executed. The sole INCLUDE set is `include_paths.txt`, one literal repository-relative file per line.

```powershell
Set-Location 'C:\Users\Administrator\Desktop\codex-workspaces\Console3\solver_v2'
git --literal-pathspecs add --pathspec-from-file=docs/baseline_20260915/include_paths.txt
```

For later user execution only. No wildcard or whole runs/final_results staging. `include_manifest.json` records bytes and SHA-256; its own hash and size are excluded to avoid self-reference. This is working-tree size, not Git pack size. Current smoke logs replace redundant historic core logs in the minimal set; A1/A2, geometry and mesh frozen evidence remain. final_results README and figure manifest remain historical provenance.
