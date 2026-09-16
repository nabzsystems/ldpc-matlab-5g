%% PART 2 — From Base Graph to Real Parity-Check Matrix
%
% Goal: expand a real 5G NR BG1 shift-value table into the full sparse
% parity-check matrix H, using the exact circular-shift mechanism from
% the article and the Manim "Lifting" scene — this time on the real
% structure, not a 2x3 toy example.
%
% IMPORTANT — READ BEFORE RUNNING:
% MATLAB's 5G Toolbox does NOT expose a documented function that returns
% the raw base-graph table directly (confirmed via MATLAB Answers). So
% this script does NOT fabricate the 46x68 shift-value table itself —
% typing ~3100 numbers from memory is exactly the kind of unverifiable
% claim this whole project has been trying to avoid.
%
% Instead, get the real table from ONE of these two independent, public
% sources (cross-check a handful of cells between them before trusting
% either):
%   1. https://github.com/manuts/NR-LDPC-BG
%      Plain-text base-graph tables, one file per (base graph, lifting-size-set).
%      Confirm the exact delimiter/format when you open the file — adjust
%      the loader below if it isn't whitespace-separated integers.
%   2. https://github.com/vodafone-chair/5g-nr-ldpc
%      A working MATLAB implementation (ldpcGenCodeMat.m) that builds H
%      directly. Useful as a second, independently-coded reference to
%      confirm a few specific cells match.
%
% Save the base-graph table from source 1 as a plain CSV/whitespace file
% at: ../data/bg1_raw_table.txt  (rows = 46, columns = 68, entries are
% either -1 for "no connection" or an integer shift amount 0..Zc-1)

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, '..', 'data');
figDir    = fullfile(scriptDir, '..', 'figures');
if ~exist(dataDir, 'dir'); mkdir(dataDir); end
if ~exist(figDir, 'dir');  mkdir(figDir);  end
addpath(scriptDir);

%% --- Load the real BG1 table (you provide this file, see note above) ---
bgFile = fullfile(dataDir, 'bg1_raw_table.txt');
assert(isfile(bgFile), ...
    ['Missing ' bgFile '. The 3GPP TS 38.212 BG1 shift table must be located at this path.']);

BG = readmatrix(bgFile);
[numRowsBG, numColsBG] = size(BG);
fprintf('Loaded base graph: %d rows x %d columns\n', numRowsBG, numColsBG);
assert(isequal(size(BG), [46 68]), ...
    'Expected a 46x68 table for BG1 — got a different size. Check the downloaded file.');

%% --- Sanity checks against known, article-verified BG1 facts ---
% These don't prove every cell is correct, but they catch gross
% transcription/parsing errors before we build a 65000+ entry matrix from it.
numConnections = sum(BG(:) ~= -1);
fprintf('Non-empty cells: %d out of %d (%.1f%% density)\n', ...
        numConnections, numel(BG), 100*numConnections/numel(BG));

Kb = 22; % systematic columns for BG1, verified in the article against Consensus research
fprintf('Expected systematic columns (Kb) = %d — verify columns 1:%d look denser than the rest.\n', Kb, Kb);

%% --- Choose a lifting size and expand ---
% Auto-sync with Part 1 if run previously, or default to Zc = 256 (Set 0)
part1File = fullfile(dataDir, 'part1_toolbox_bler.mat');
if isfile(part1File)
    p1 = load(part1File, 'cbsInfo');
    Zc = p1.cbsInfo.Zc;
    fprintf('Auto-synced lifting size Zc = %d from Part 1 (cbsInfo.Zc).\n', Zc);
else
    Zc = 256; % matches 3GPP BG1 Set 0 (bg1_raw_table.txt)
    fprintf('Using default lifting size Zc = %d (Set 0).\n', Zc);
end

H = expand_base_graph(BG, Zc);
fprintf('\nExpanded H: %d x %d, %d nonzeros (%.3f%% density)\n', ...
        size(H,1), size(H,2), nnz(H), 100*nnz(H)/numel(H));

%% --- Visualize sparsity pattern ---
figure('Color', [1 1 1]);
spy(H);
title(sprintf('Real 5G NR BG1 parity-check matrix, Z_c = %d (%dx%d)', Zc, size(H,1), size(H,2)));
saveas(gcf, fullfile(figDir, 'part2_bg1_sparsity.png'));

save(fullfile(dataDir, 'part2_H_matrix.mat'), 'H', 'BG', 'Zc');
fprintf('\nSaved: %s and %s\n', ...
        fullfile(dataDir, 'part2_H_matrix.mat'), fullfile(figDir, 'part2_bg1_sparsity.png'));

%% --- Local function: the verified circular-shift expansion ---
function H = expand_base_graph(BG, Z)
    % Same mechanism as the article's Section 3 and the Manim "Lifting" scene:
    % each -1 cell becomes a ZxZ zero block; each shift-value cell s becomes
    % a ZxZ identity matrix with rows cyclically shifted by s.
    [numBGrows, numBGcols] = size(BG);
    H = sparse(numBGrows*Z, numBGcols*Z);
    I = eye(Z);
    for r = 1:numBGrows
        for c = 1:numBGcols
            shift = BG(r,c);
            if shift >= 0
                block = circshift(I, [0, mod(shift, Z)]);
                H((r-1)*Z+1:r*Z, (c-1)*Z+1:c*Z) = block;
            end
            % shift == -1: leave as the pre-allocated zero block
        end
    end
end
