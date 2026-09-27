# ScanFold2 source and patch

This directory pins [moss-lab/ScanFold2.0](https://github.com/moss-lab/ScanFold2.0)
at commit `c5cc73291fa5a1c06bf5b3c8367a396a3238141c` and applies
`multifasta.patch`. The patch incorporates filename collision fixes from the
[local fork](https://github.com/dgruano/ScanFold2.0) (commits `5e828af` and
`7e474e4`) and subsequent local fixes. Upstream already accepts multiple FASTA
records; the patch gives each record its own auxiliary filenames, validates
record IDs before processing, and records source hashes in run logs. It does
not change the folding or statistical model.

`source-lock.json` records the expected source and model hashes and the patch
checksum. The upstream MIT license is included as `LICENSE`.

`historical-local.patch` and `historical-source.json` preserve an earlier local
source state. They are archival material; new installations use
`multifasta.patch` and `source-lock.json`.
