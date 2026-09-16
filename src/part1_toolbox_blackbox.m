%% PART 1 — The Black Box: Official 5G Toolbox LDPC Encode/Decode + BLER
%
% Goal: use ONLY the documented, official 5G Toolbox functions to encode,
% transmit over an AWGN channel, decode, and measure Block Error Rate (BLER)
% vs Eb/No for a real BG1-selected configuration. This is the "trust the
% black box" baseline that Part 4 will later overlay against our own
% from-scratch decoder.
%
% Requires: 5G Toolbox (nrLDPCEncode, nrLDPCDecode, nrCodeBlockSegmentLDPC,
% nrCodeBlockDesegmentLDPC, nrCRCEncode, nrCRCDecode, nrDLSCHInfo,
% nrRateMatchLDPC, nrRateRecoverLDPC, nrDLSCHInfo)
%
% NOTE FOR RECORDING: run this once end-to-end BEFORE writing the
% screen-recording script, and save the printed BGN/Zc values + the final
% figure — the video script's narration should quote the real numbers this
% produces, not the numbers below (those are only a reasonable starting
% guess to force BG1 selection).

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, '..', 'data');
figDir    = fullfile(scriptDir, '..', 'figures');
if ~exist(dataDir, 'dir'); mkdir(dataDir); end
if ~exist(figDir, 'dir');  mkdir(figDir);  end
addpath(scriptDir);
rng(42); % reproducible for the video — same run every time you re-record

%% --- Choose parameters that select BG1 (not BG2) ---
% BG1 is selected for larger transport blocks / higher code rates.
% TBS=5300, R=1/2 selects BG1 with lifting size Zc=256 (Set 0), matching bg1_raw_table.txt.
TBS      = 5300;        % transport block size (bits)
codeRate = 1/2;         % target code rate
rv       = 0;           % redundancy version (first transmission)
modulation = 'QPSK';
nlayers  = 1;

cbsInfo = nrDLSCHInfo(TBS, codeRate);
fprintf('--- DL-SCH coding parameters ---\n');
disp(cbsInfo);

assert(cbsInfo.BGN == 1, ...
    ['Selected BGN is ' num2str(cbsInfo.BGN) ', not 1 — increase TBS or codeRate and re-run.']);

fprintf('Confirmed: this configuration uses Base Graph 1 (BG1).\n');
fprintf('Lifting size Zc = %d\n', cbsInfo.Zc);
fprintf('Codeword length per block N = %d\n\n', cbsInfo.N);

%% --- Sweep Eb/No, measure BLER ---
EbNoRange   = -2:0.5:4;   % dB — adjust range after first run based on where BLER actually falls
maxBlocksPerPoint = 2000; % increase for smoother curve, decrease for faster first pass
minErrorsPerPoint = 50;   % stop early once enough errors collected (standard MC practice)

bitsPerSymbol = 2; % QPSK

BLER = zeros(size(EbNoRange));

for idx = 1:length(EbNoRange)
    EbNo = EbNoRange(idx);
    numBlockErrors = 0;
    numBlocksSimulated = 0;

    while numBlockErrors < minErrorsPerPoint && numBlocksSimulated < maxBlocksPerPoint

        % --- generate + encode ---
        in = randi([0 1], TBS, 1, 'int8');
        tbIn = nrCRCEncode(in, cbsInfo.CRC);
        cbsIn = nrCodeBlockSegmentLDPC(tbIn, cbsInfo.BGN);
        enc = nrLDPCEncode(cbsIn, cbsInfo.BGN);

        outlen = ceil(TBS / codeRate);
        chIn = nrRateMatchLDPC(enc, outlen, rv, modulation, nlayers);

        % --- BPSK-equivalent modulation (QPSK per-bit is equiv. to BPSK
        % per real/imag component for this LLR-based BLER test) ---
        txSym = 1 - 2*double(chIn); % 0 -> +1, 1 -> -1

        snrdB = convertSNR(EbNo, 'ebno', 'BitsPerSymbol', bitsPerSymbol, ...
                            'CodingRate', TBS/outlen);
        noiseVar = 10^(-snrdB/10);
        noise = sqrt(noiseVar/2) * randn(size(txSym));
        rxSym = txSym + noise;

        % --- LLR demodulation (BPSK-equivalent) ---
        rxLLR = 2*rxSym/noiseVar;

        % --- rate recovery + decode ---
        raterec = nrRateRecoverLDPC(rxLLR, TBS, codeRate, rv, modulation, nlayers);
        [decBits, actNumIter, finalParity] = nrLDPCDecode(raterec, cbsInfo.BGN, 25, ...
                                                             'Algorithm', 'Belief propagation');
        [blk, crcErr] = nrCodeBlockDesegmentLDPC(decBits, cbsInfo.BGN, TBS + cbsInfo.L);
        [~, tbCrcErr] = nrCRCDecode(blk, cbsInfo.CRC);

        blockInError = (tbCrcErr ~= 0);
        numBlockErrors = numBlockErrors + blockInError;
        numBlocksSimulated = numBlocksSimulated + 1;
    end

    BLER(idx) = numBlockErrors / numBlocksSimulated;
    fprintf('Eb/No = %5.1f dB | blocks = %5d | errors = %4d | BLER = %.4e\n', ...
            EbNo, numBlocksSimulated, numBlockErrors, BLER(idx));
end

%% --- Plot ---
figure('Color', [1 1 1]);
semilogy(EbNoRange, BLER, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Block Error Rate (BLER)');
title(sprintf('5G NR LDPC BG1, official Toolbox decoder (TBS=%d, R=%.2f)', TBS, codeRate));

% save for use in Part 4's overlay
save(fullfile(dataDir, 'part1_toolbox_bler.mat'), 'EbNoRange', 'BLER', 'TBS', 'codeRate', 'cbsInfo');
saveas(gcf, fullfile(figDir, 'part1_toolbox_bler.png'));

fprintf('\nSaved: %s and %s\n', ...
        fullfile(dataDir, 'part1_toolbox_bler.mat'), fullfile(figDir, 'part1_toolbox_bler.png'));
fprintf('IMPORTANT: note the actual EbNo range where BLER drops from ~1 to ~1e-4 —\n');
fprintf('this range is what the video narration should reference, not a guess.\n');
