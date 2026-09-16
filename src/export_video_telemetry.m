%% EXPORT_VIDEO_TELEMETRY — Ground-Truth Telemetry Exporter for 5G Manim Video
%
% Generates deterministic, single-block telemetry data matching the 4 acts
% of the "Underwaterfall" video series (../5Gvide):
%   Act I  : Channel Noise — Raw complex constellation scatter at Eb/No ~ 1.2 dB
%   Act II : Belief Propagation — Iteration-by-iteration soft LLR tensor evolution
%   Act III: Convergence — Rapid collapse of unsatisfied parity checks (s = H*c' mod 2)
%   Act IV : Recovery — Reconstructed clean constellation, zero bit errors, CRC pass
%
% 3GPP Specification Baseline:
%   - 3GPP TS 38.212 Section 5.3.2 (LDPC coding, Base Graph 1, Zc = 384)
%   - 3GPP TS 38.211 Section 5.1.3 (QPSK modulation mapping)
%   - Code block: K = 8448 bits (TBS = 8424 bits + CRC-24A, 0 filler bits, C = 1)
%   - Code rate R = 1/2, rv = 0
%   - Target operational SNR: Eb/No = 1.2 dB
%
% Output:
%   - ../data/telemetry_underwaterfall.mat
%   - ../data/telemetry_underwaterfall.npz (via Python bridge)

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, '..', 'data');
if ~exist(dataDir, 'dir'); mkdir(dataDir); end
addpath(scriptDir);
rng(42); % Fixed deterministic seed for reproducible video ground truth

%% 1. Transmission Parameters (3GPP TS 38.212 / TS 38.211)
TBS           = 8424;         % Transport block size (bits)
codeRate      = 1/2;          % Target code rate
rv            = 0;            % Redundancy version 0 (initial transmission)
modulation    = 'QPSK';
bitsPerSymbol = 2;            % QPSK
nlayers       = 1;
EbNo          = 1.2;          % dB (waterfall transition region)

cbsInfo = nrDLSCHInfo(TBS, codeRate);
assert(cbsInfo.BGN == 1, 'Expected Base Graph 1.');
assert(cbsInfo.Zc == 384, 'Expected maximum lifting size Zc = 384.');
assert(cbsInfo.K == 8448, 'Expected code block size K = 8448 bits.');
assert(cbsInfo.C == 1, 'Expected single code block (C = 1).');

fprintf('============================================================\n');
fprintf('  5G NR LDPC VIDEO TELEMETRY GENERATION (3GPP TS 38.212)    \n');
fprintf('============================================================\n');
fprintf('Configuration: BG%d, Zc=%d, TBS=%d bits, Code Rate=%.2f\n', ...
        cbsInfo.BGN, cbsInfo.Zc, TBS, codeRate);
fprintf('Total systematic bits K = %d (CRC: %s, %d bits, filler: %d bits)\n', ...
        cbsInfo.K, cbsInfo.CRC, cbsInfo.L, cbsInfo.F);

%% 2. Code Block Generation, CRC Attachment & LDPC Encoding
inBits = randi([0 1], TBS, 1, 'int8');
tbIn   = nrCRCEncode(inBits, cbsInfo.CRC);
cbsIn  = nrCodeBlockSegmentLDPC(tbIn, cbsInfo.BGN);
enc    = nrLDPCEncode(cbsIn, cbsInfo.BGN);

outlen = ceil(TBS / codeRate);
chIn   = nrRateMatchLDPC(enc, outlen, rv, modulation, nlayers);

%% 3. QPSK Modulation & AWGN Channel (Act I: Channel Noise)
% Mapping per 3GPP TS 38.211 Section 5.1.3:
%   d(i) = (1/sqrt(2)) * ((1 - 2*b(2i)) + 1j*(1 - 2*b(2i+1)))
tx_symbols = ((1 - 2*double(chIn(1:2:end))) + 1j*(1 - 2*double(chIn(2:2:end)))) / sqrt(2);

snrdB    = convertSNR(EbNo, 'ebno', 'BitsPerSymbol', bitsPerSymbol, 'CodingRate', TBS/outlen);
noiseVar = 10^(-snrdB/10);
noise    = sqrt(noiseVar/2) * (randn(size(tx_symbols)) + 1j*randn(size(tx_symbols)));
rx_symbols = tx_symbols + noise;

fprintf('Channel: AWGN, Eb/No = %.2f dB (SNR = %.2f dB, noiseVar = %.4f)\n', ...
        EbNo, snrdB, noiseVar);
fprintf('Transmitted/Received complex QPSK symbols: %d symbols\n', length(tx_symbols));

%% 4. Demodulation & Rate Recovery
rxDeint = zeros(size(chIn));
rxDeint(1:2:end) = real(rx_symbols) * sqrt(2);
rxDeint(2:2:end) = imag(rx_symbols) * sqrt(2);
rxLLR = 2 * rxDeint / noiseVar;

raterec = nrRateRecoverLDPC(rxLLR, TBS, codeRate, rv, modulation, nlayers);

%% 5. Belief Propagation Iteration Stepping (Act II & III: Convergence)
% Determine full convergence iteration with early termination:
[decBitsFinal, finalIter, finalParityChecks] = nrLDPCDecode(raterec, cbsInfo.BGN, 25, ...
    'Algorithm', 'Belief propagation', ...
    'OutputFormat', 'info', ...
    'DecisionType', 'hard', ...
    'Termination', 'early');

[blkFinal, crcErr] = nrCodeBlockDesegmentLDPC(decBitsFinal, cbsInfo.BGN, TBS + cbsInfo.L);
[~, tbCrcErr] = nrCRCDecode(blkFinal, cbsInfo.CRC);

fprintf('\nOfficial BP Early Termination:\n');
fprintf('  Converged at iteration: %d\n', finalIter);
fprintf('  Final unsatisfied parity checks: %d\n', sum(finalParityChecks ~= 0));
fprintf('  Transport block CRC error: %d (0 = PASS)\n', tbCrcErr);
assert(tbCrcErr == 0, 'Transport block CRC failed!');

% Step through iterations 1 to finalIter to record full telemetry history:
max_iter_record = finalIter;
N_full = cbsInfo.N + 2*cbsInfo.Zc; % 26112 (full lifted graph nodes)

llr_history      = zeros(max_iter_record, N_full);
syndrome_history = zeros(max_iter_record, 1);

fprintf('\nRecording iteration-by-iteration telemetry:\n');
for it = 1:max_iter_record
    [softCodeword, ~, parityChecks] = nrLDPCDecode(raterec, cbsInfo.BGN, it, ...
        'Algorithm', 'Belief propagation', ...
        'OutputFormat', 'whole', ...
        'DecisionType', 'soft', ...
        'Termination', 'max');
    
    llr_history(it, :)      = softCodeword';
    syndrome_history(it)    = sum(parityChecks ~= 0);
    
    % Systematic hard-decision check at this iteration
    hardCodeword = double(softCodeword < 0);
    sysBitsIter = hardCodeword(1:cbsInfo.K);
    bitErrorsIter = sum(sysBitsIter ~= double(tbIn));
    
    fprintf('  Iter %2d: unsat parity checks = %5d | bit errors = %4d | mean |LLR| = %.2f\n', ...
            it, syndrome_history(it), bitErrorsIter, mean(abs(softCodeword)));
end

%% 6. Recovery & Verification (Act IV: Zero Bit Errors)
decoded_bits = double(decBitsFinal);
tx_systematic_bits = double(tbIn(1:TBS));
recovered_systematic_bits = decoded_bits(1:TBS);
bit_errors = sum(tx_systematic_bits ~= recovered_systematic_bits);

fprintf('\nAct IV Verification:\n');
fprintf('  Systematic payload bits: %d\n', TBS);
fprintf('  Total bit errors: %d\n', bit_errors);
assert(bit_errors == 0, 'Bit error count is nonzero!');
assert(syndrome_history(end) == 0, 'Final syndrome did not converge to zero!');

%% 7. Parity-Check Matrix H Construction (QC-LDPC BG1 Set 2)
bgs = load('nr5g/internal/ldpc/baseGraph');
P_shifts = nr5g.internal.ldpc.calcShiftValues(bgs.BG1S2, cbsInfo.Zc);
H_sparse = ldpcQuasiCyclicMatrix(cbsInfo.Zc, P_shifts);
fprintf('Parity-check matrix H: %d x %d (sparse, %d nonzeros)\n', ...
        size(H_sparse,1), size(H_sparse,2), nnz(H_sparse));

%% 8. Export Telemetry to MAT-file
matFile = fullfile(dataDir, 'telemetry_underwaterfall.mat');
npzFile = fullfile(dataDir, 'telemetry_underwaterfall.npz');
save(matFile, ...
     'tx_symbols', 'rx_symbols', ...
     'llr_history', 'syndrome_history', ...
     'tx_systematic_bits', 'recovered_systematic_bits', ...
     'decoded_bits', 'bit_errors', 'finalIter', ...
     'TBS', 'codeRate', 'EbNo', 'snrdB', 'noiseVar', ...
     'modulation', 'cbsInfo', 'H_sparse');
fprintf('\nSaved MATLAB telemetry to: %s\n', matFile);

%% 9. Bridge: Export NumPy (.npz) for Manim Pipeline
cmd = sprintf('python3 -c "import scipy.io as sio, numpy as np; m = sio.loadmat(''%s''); np.savez_compressed(''%s'', tx_symbols=m[''tx_symbols''], rx_symbols=m[''rx_symbols''], llr_history=m[''llr_history''], syndrome_history=m[''syndrome_history''], tx_bits=m[''tx_systematic_bits''], decoded_bits=m[''decoded_bits''], final_iter=int(m[''finalIter''][0,0]), ebno=float(m[''EbNo''][0,0]), snr_db=float(m[''snrdB''][0,0]), noise_var=float(m[''noiseVar''][0,0])); print(''NPZ export complete.'')" ', matFile, npzFile);
system(cmd);

fprintf('Video ground-truth telemetry export SUCCESSFUL.\n');
