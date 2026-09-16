# LDPC MATLAB 5G Toolbox

A from-scratch exploration of 5G NR LDPC coding (3GPP TS 38.212, Base Graph 1), built to produce **real, verifiable numbers** for the نبض سیستم (Nabz System) article and video series on LDPC and Polar codes — rather than relying on toy examples or unverified claims.

The project compares MATLAB's official 5G Toolbox LDPC pipeline against a from-scratch min-sum decoder, on the real BG1 parity-check structure, and reports honest findings (including where they diverge).

## Requirements

- MATLAB with the **5G Toolbox** (`nrLDPCEncode`, `nrLDPCDecode`, `nrRateMatchLDPC`, `nrRateRecoverLDPC`, `nrCRCEncode`/`nrCRCDecode`, `nrCodeBlockSegmentLDPC`/`nrCodeBlockDesegmentLDPC`, `nrDLSCHInfo`)
- Python 3 with `scipy` and `numpy` (only for the optional `.npz` export bridge in `export_video_telemetry.m`)

## Project structure

```
data/    Input tables and generated intermediate results (.mat/.npz are git-ignored)
figures/ Generated plots
src/     All MATLAB scripts and functions
```

## Pipeline — run in this order

| Step | Script | Purpose |
|---|---|---|
| 1 | `part1_toolbox_blackbox.m` | Baseline: encode/transmit/decode/measure BLER using **only** the official 5G Toolbox functions. |
| 2 | `part2_extract_bg1.m` | Expand the real 3GPP BG1 shift-value table (`data/bg1_raw_table.txt`) into the full sparse parity-check matrix `H`. |
| 3 | `part3_manual_decoder.m` | Self-test: validates `decode_min_sum.m` against a hand-worked example before trusting it on real data. |
| 4 | `part4_curve_overlay.m` | Overlays the official decoder's BLER curve against the from-scratch decoder on the same real `H`. |
| 5 | `part5_error_floor_montecarlo.m` | Large-scale Monte Carlo search for a real error floor at high Eb/N0. |
| — | `export_video_telemetry.m` | Generates deterministic, single-block telemetry (constellation, LLR history, syndrome history) used for the Manim "Underwaterfall" animation. |

`decode_min_sum.m` is the shared from-scratch min-sum LDPC decoder used by Parts 3–5.

## Important caveats

- **`data/bg1_raw_table.txt` is third-party input, not something we hand-derived.** MATLAB's 5G Toolbox doesn't expose the raw BG1 shift-value table directly, so this file was sourced externally (see comments in `part2_extract_bg1.m` for the two independent sources used) and should be spot-checked against both before being trusted for anything beyond this project.
- **Part 5 is computationally heavy.** The manual decoder's check/variable-node update loops are not fully vectorized; a full Monte Carlo run at realistic block counts can take a long time. Time a small run first (few blocks, one Eb/N0 point) before committing to the full sweep.
- Findings from Part 4 and Part 5 (curve overlap, presence/absence of an error floor) are reported as-observed, whatever they turn out to be — including a prior toy 4-cycle example (referenced in `part3`) that was found to *not* show a persistent error floor, unlike what an earlier draft of the video assumed.

## AI-assisted code disclosure

This codebase was written with the assistance of Gemini Flash 3.8. It has been reviewed, but has **not been independently verified line-by-line against the 3GPP specification** by a domain expert. If you plan to reuse this code for anything beyond educational/illustrative purposes, verify the parity-check construction, puncturing/rate-matching logic, and `bg1_raw_table.txt` contents independently. Issues and corrections are welcome.


## License

This project is licensed under the [MIT License](LICENSE).