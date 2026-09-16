%% PART 4 — The Moment of Truth: Overlay Our Decoder vs. the Official One
%
% Goal: encode and transmit with the OFFICIAL, trusted 5G Toolbox pipeline
% (same as Part 1 — encoding correctness isn't what we're testing), but
% decode every received block TWICE: once with nrLDPCDecode, once with
% our own decode_min_sum from Part 3. If the curves overlap, we've
% confirmed our from-scratch implementation actually behaves like the
% real thing on the real BG1 structure — not just on the small toy example.
%
% Run Part 1, Part 2, and Part 3 (in that order) before this script.

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, '..', 'data');
figDir    = fullfile(scriptDir, '..', 'figures');
if ~exist(dataDir, 'dir'); mkdir(dataDir); end
if ~exist(figDir, 'dir');  mkdir(figDir);  end
addpath(scriptDir);
rng(42); % same seed as Part 1, for a fair side-by-side comparison

load(fullfile(dataDir, 'part1_toolbox_bler.mat'), 'TBS', 'codeRate', 'cbsInfo');
load(fullfile(dataDir, 'part2_H_matrix.mat'), 'H', 'Zc');

assert(Zc == cbsInfo.Zc, ...
    sprintf(['Part 2''s Zc (%d) does not match Part 1''s actual lifting size (%d). ' ...
             'Re-run Part 2 with Zc = %d for a fair comparison.'], Zc, cbsInfo.Zc, cbsInfo.Zc));

rv = 0;
modulation = 'QPSK';
nlayers = 1;
bitsPerSymbol = 2;

EbNoRange = -2:0.5:4; % keep identical to Part 1 for a fair overlay
maxBlocksPerPoint = 500;   % lower than Part 1 — our manual decoder is much
                            % slower than the compiled Toolbox one; raise
                            % this once you've confirmed the loop timing
minErrorsPerPoint = 30;
maxManualIter = 25;        % match nrLDPCDecode's iteration cap from Part 1

BLER_official = zeros(size(EbNoRange));
BLER_manual   = zeros(size(EbNoRange));

for idx = 1:length(EbNoRange)
    EbNo = EbNoRange(idx);
    errOfficial = 0; errManual = 0; numBlocks = 0;

    while numBlocks < maxBlocksPerPoint && ...
          (errOfficial < minErrorsPerPoint || errManual < minErrorsPerPoint)

        in = randi([0 1], TBS, 1, 'int8');
        tbIn = nrCRCEncode(in, cbsInfo.CRC);
        cbsIn = nrCodeBlockSegmentLDPC(tbIn, cbsInfo.BGN);
        enc = nrLDPCEncode(cbsIn, cbsInfo.BGN);

        outlen = ceil(TBS / codeRate);
        chIn = nrRateMatchLDPC(enc, outlen, rv, modulation, nlayers);

        txSym = 1 - 2*double(chIn);
        snrdB = convertSNR(EbNo, 'ebno', 'BitsPerSymbol', bitsPerSymbol, ...
                            'CodingRate', TBS/outlen);
        noiseVar = 10^(-snrdB/10);
        noise = sqrt(noiseVar/2) * randn(size(txSym));
        rxSym = txSym + noise;
        rxLLR = 2*rxSym/noiseVar;

        raterec = nrRateRecoverLDPC(rxLLR, TBS, codeRate, rv, modulation, nlayers);

        % --- official decode ---
        [decBitsOff, ~, ~] = nrLDPCDecode(raterec, cbsInfo.BGN, 25, 'Algorithm', 'Belief propagation');
        [blkOff, ~] = nrCodeBlockDesegmentLDPC(decBitsOff, cbsInfo.BGN, TBS + cbsInfo.L);
        [~, crcOff] = nrCRCDecode(blkOff, cbsInfo.CRC);

        % --- our decode (on the FIRST code block only, for speed —
        % raterec may contain multiple concatenated code blocks; adjust
        % indexing here if cbsInfo.C > 1 for your TBS/codeRate choice) ---
        assert(cbsInfo.C == 1, ...
            'This script assumes a single code block (cbsInfo.C == 1). Adjust indexing if not.');

        % Reinsert 2*Zc zero-LLRs for the punctured first two systematic columns (TS 38.212 §5.3.2)
        % Sign convention: positive LLR = bit 0 (matches both nrLDPCDecode and decode_min_sum)
        llrForManual = [zeros(1, 2*Zc), raterec(1:cbsInfo.N)'];
        [decBitsManual, ~, ~] = decode_min_sum(llrForManual, H, maxManualIter);

        % Extract only the systematic bits (first K = 22*Zc bits) for desegmentation & CRC check
        sysBitsManual = decBitsManual(1:cbsInfo.K)';
        [blkManual, ~] = nrCodeBlockDesegmentLDPC(int8(sysBitsManual), cbsInfo.BGN, TBS + cbsInfo.L);
        [~, crcManual] = nrCRCDecode(blkManual, cbsInfo.CRC);

        errOfficial = errOfficial + (crcOff ~= 0);
        errManual   = errManual   + (crcManual ~= 0);
        numBlocks = numBlocks + 1;
    end

    BLER_official(idx) = errOfficial / numBlocks;
    BLER_manual(idx)   = errManual / numBlocks;
    fprintf('Eb/No=%5.1f dB | blocks=%4d | official BLER=%.3e | manual BLER=%.3e\n', ...
            EbNo, numBlocks, BLER_official(idx), BLER_manual(idx));
end

%% --- Overlay plot ---
figure('Color', [1 1 1]);
semilogy(EbNoRange, BLER_official, 'o-', 'LineWidth', 2, 'DisplayName', 'Official nrLDPCDecode');
hold on;
semilogy(EbNoRange, BLER_manual, 's--', 'LineWidth', 2, 'DisplayName', 'Our from-scratch min-sum');
grid on;
legend('Location', 'southwest');
xlabel('E_b/N_0 (dB)');
ylabel('Block Error Rate (BLER)');
title('Official Toolbox decoder vs. our from-scratch implementation, same real BG1 H');

saveas(gcf, fullfile(figDir, 'part4_curve_overlay.png'));
save(fullfile(dataDir, 'part4_overlay_results.mat'), 'EbNoRange', 'BLER_official', 'BLER_manual');

fprintf('\nSaved: %s and %s\n', ...
        fullfile(figDir, 'part4_curve_overlay.png'), fullfile(dataDir, 'part4_overlay_results.mat'));
fprintf('IMPORTANT: report how closely the two curves actually track — that''s the real\n');
fprintf('finding for the video narration, whatever it turns out to be.\n');
