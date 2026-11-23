#!/usr/bin/env bash
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
mathlib_dir="${1:?Usage: ./verify.sh /absolute/path/to/mathlib4}"
expected_rev=fc4f59e577feb1a4956671bc04ac18aaaf8cd842
actual_rev="$(git -C "$mathlib_dir" rev-parse HEAD)"
if [[ "$actual_rev" != "$expected_rev" ]]; then
  echo "Mathlib revision mismatch: $actual_rev" >&2
  exit 1
fi
cd "$mathlib_dir"
lake env lean --version
lake env lean -DwarningAsError=true "$package_dir/Brownian_Compensator_Verified_Components.lean"
echo 'PASS: all 39 included theorem statements compiled, with warnings treated as errors.'
echo 'Coverage: see VERIFICATION.md; the complete compensator martingale is an open formalization obligation.'
