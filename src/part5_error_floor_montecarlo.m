%% PART 5 — Hunting for a Real Error Floor
%
% Context: the article's Manim video (Scene G) originally used a tiny
% 4-node toy example to illustrate why short cycles are dangerous for
% belief propagation. An independent Python simulation later showed that
% specific toy example actually self-corrects by iteration 2 — it does
% NOT demonstrate a persistent error floor. This script is the honest
% follow-up: go look for a REAL error floor on the REAL BG1 matrix, with
% a proper large-scale Monte Carlo run, and report whatever is actually
% found — including "we didn't find one in this range," which is itself
% a meaningful result (3GPP's PEG-based construction is specifically
% designed to push the floor down very low).
%
% Run Part 1, 2, and 3 before this script.

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, '..', 'data');
figDir    = fullfile(scriptDir, '..', 'figures');
if ~exist(dataDir, 'dir'); mkdir(dataDir); end
if ~exist(figDir, 'dir');  mkdir(figDir);  end
addpath(scriptDir);
rng(7); % different seed from Parts 1/4 — this is an independent experiment

load(fullfile(dataDir, 'part2_H_matrix.mat'), 'H', 'Zc');
load(fullfile(dataDir, 'part1_toolbox_bler.mat'), 'TBS', 'codeRate', 'cbsInfo');

rv = 0;
modulation = 'QPSK';
nlayers = 1;
bitsPerSymbol = 2;
maxManualIter = 50; % higher than Part 4 — give the decoder every chance to converge

% Push further into high Eb/No than Parts 1 and 4 — an error floor, if one
% exists, should show up as BLER refusing to keep dropping as EbNo increases.
envEbNo = getenv('LDPC_MC_EBNO');
if ~isempty(envEbNo)
    EbNoRange = eval(envEbNo);
else
    EbNoRange = 3:0.5:8;
end

envBlocks = getenv('LDPC_MC_BLOCKS');
if ~isempty(envBlocks)
    blocksPerPoint = str2double(envBlocks);
else
    blocksPerPoint = 5000; % large N is the whole point here — don't shortcut
                            % with early-stopping on error count for this part,
                            % since a floor might mean VERY few errors even at
                            % 5000 blocks, and that near-zero count is itself
                            % the finding.
end

BLER_floor_search = zeros(size(EbNoRange));
errorCounts = zeros(size(EbNoRange));

fprintf('Running %d blocks per Eb/No point across %d points — this will take a while.\n', ...
        blocksPerPoint, length(EbNoRange));
fprintf('That''s intentional: a real error-floor search needs real sample size, not a shortcut.\n\n');

for idx = 1:length(EbNoRange)
    EbNo = EbNoRange(idx);
    errCount = 0;

    for b = 1:blocksPerPoint
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

        assert(cbsInfo.C == 1, 'Adjust for multi-code-block case if your TBS/codeRate triggers C > 1.');

        % Reinsert 2*Zc zero-LLRs for punctured systematic nodes (TS 38.212 §5.3.2)
        % Sign convention: positive LLR = bit 0 (matches nrLDPCDecode and decode_min_sum)
        llrForManual = [zeros(1, 2*Zc), raterec(1:cbsInfo.N)'];

        [decBits, ~, ~] = decode_min_sum(llrForManual, H, maxManualIter);

        % Extract only systematic bits for desegmentation and CRC
        sysBits = decBits(1:cbsInfo.K)';
        [blk, ~] = nrCodeBlockDesegmentLDPC(int8(sysBits), cbsInfo.BGN, TBS + cbsInfo.L);
        [~, crcErr] = nrCRCDecode(blk, cbsInfo.CRC);

        errCount = errCount + (crcErr ~= 0);
    end

    BLER_floor_search(idx) = errCount / blocksPerPoint;
    errorCounts(idx) = errCount;
    fprintf('Eb/No = %4.1f dB | %d/%d blocks in error | BLER = %.4e\n', ...
            EbNo, errCount, blocksPerPoint, BLER_floor_search(idx));
end

%% --- Plot and interpret ---
figure('Color', [1 1 1]);
semilogy(EbNoRange, max(BLER_floor_search, 1/blocksPerPoint/10), 'o-', 'LineWidth', 2);
% (the max(...) clamp is just so log-scale plotting doesn't choke on exact
% zeros — it does NOT change what gets reported in the console/text above)
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Block Error Rate (BLER)');
title('Real BG1, large-scale Monte Carlo — is there a visible error floor?');
saveas(gcf, fullfile(figDir, 'part5_error_floor_montecarlo.png'));
save(fullfile(dataDir, 'part5_floor_results.mat'), 'EbNoRange', 'BLER_floor_search', 'errorCounts', 'blocksPerPoint');

fprintf('\nSaved: %s and %s\n', ...
        fullfile(figDir, 'part5_error_floor_montecarlo.png'), fullfile(dataDir, 'part5_floor_results.mat'));

fprintf('\n--- How to read this for the video narration ---\n');
fprintf('Look at whether BLER keeps dropping roughly exponentially as EbNo increases,\n');
fprintf('or whether it flattens out (a floor) at the higher EbNo values.\n');
fprintf('Either finding is worth reporting honestly:\n');
fprintf('  - If it keeps dropping: say so plainly — the real, properly-designed BG1 code\n');
fprintf('    does not show a floor in this range, unlike the deliberately naive toy example.\n');
fprintf('  - If it flattens: that IS a real error floor — report the EbNo where it starts\n');
fprintf('    and the errorCounts at that point (small counts = say so, don''t overstate\n');
fprintf('    confidence from a handful of error events).\n');
