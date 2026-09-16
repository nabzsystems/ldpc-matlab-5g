%% PART 3 — Self-Test: Does Our From-Scratch Decoder Match the Article?
%
% The actual decoder lives in utils/decode_min_sum.m (MATLAB requires
% reusable functions to live in their own file, not at the bottom of a
% script, if part4/part5 need to call them too).
%
% This script just runs that function against the article's own
% hand-worked Scene B example and checks the output matches exactly —
% run this BEFORE part4/part5, and don't trust either of those if this
% self-test fails.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

%% --- Self-test on the article's own hand-verifiable example ---
H_test = [1 1 0 1 0 0;
          0 1 1 0 1 0;
          1 0 1 0 0 1];

trueBits = [1 1 0 0 1 1];
corruptedBits = [1 0 0 0 1 1]; % bit 2 flipped, matches article's Scene B

% Our convention: positive LLR = leans bit 0, negative = leans bit 1.
llr_test = -(2*corruptedBits - 1) * 3.0; % confident magnitude, correct sign per convention

fprintf('--- Self-test: min-sum decoder vs. the article''s hand-worked example ---\n');
[decoded, iters, converged] = decode_min_sum(llr_test, H_test, 20);
fprintf('Decoded bits:  %s\n', mat2str(decoded));
fprintf('True bits:     %s\n', mat2str(trueBits));
fprintf('Converged in %d iterations: %d\n', iters, converged);

assert(isequal(decoded, trueBits), ...
    'Self-test FAILED — decoder does not match the article''s hand-worked derivation. Do not proceed to Part 4/5.');
fprintf('Self-test PASSED.\n\n');

%% --- Second self-test: the Scene G trapping-set example ---
% Cross-check against this project's own earlier Python simulation, which
% showed this specific 4-cycle example does NOT produce a persistent
% error floor — both v1 and v2 converge to the correct bit (0) by
% iteration 2 and stay there. This MATLAB decoder should show the same.
H_cycle = [1 1 1 0;   % c1: v1,v2,v3
           1 1 0 1];  % c2: v1,v2,v4
llr_cycle = [0.5 -0.3 2.0 1.8];

fprintf('--- Second self-test: Scene G trapping-set example ---\n');
fprintf('(Should converge to [0 0 0 0] by ~iteration 2, per this project''s Python check.)\n');
for testIter = [1 2 3 5 10]
    [decoded_c, ~, ~] = decode_min_sum(llr_cycle, H_cycle, testIter);
    fprintf('  after %2d iteration(s): %s\n', testIter, mat2str(decoded_c));
end
fprintf('If this does NOT match the earlier Python trace, there is a bug — fix before Part 5.\n');
