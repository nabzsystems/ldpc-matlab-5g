function [decodedBits, iterUsed, converged, telemetry] = decode_min_sum(llrChannel, H, maxIter)
% DECODE_MIN_SUM  Flooding-schedule min-sum LDPC decoder, written from scratch.
%
%   [decodedBits, iterUsed, converged] = decode_min_sum(llrChannel, H, maxIter)
%   [decodedBits, iterUsed, converged, telemetry] = decode_min_sum(llrChannel, H, maxIter)
%
%   llrChannel : 1xN vector of channel LLRs (positive = leans bit 0, our
%                convention throughout this project — matches the article
%                and the Manim "Error Floor" scene)
%   H          : MxN sparse or full parity-check matrix (0/1 entries)
%   maxIter    : maximum number of iterations before giving up
%
%   decodedBits : 1xN hard-decision bits
%   iterUsed    : number of iterations actually run
%   converged   : true if a valid codeword (all parity checks satisfied)
%                 was found before maxIter was reached
%   telemetry   : (optional) struct containing:
%                   .llr_history     : [iterUsed, N] matrix of total LLRs
%                   .syndrome_history: [iterUsed, 1] vector of unsatisfied parity checks
%                   .converged       : logical convergence flag
%                   .iterUsed        : iterations taken
%
%   Algorithm matches the hand-worked derivation in the article's Scene G
%   and this project's independent Python verification — see
%   part3_manual_decoder.m for a self-test against that hand derivation
%   before trusting this on real BG1-sized matrices.

    [M, N] = size(H);
    llrChannel = llrChannel(:)'; % row vector

    % Persistent cache for precomputed bipartite edge indexing.
    % Avoids repeating O(E) accumarray/sorting across thousands of Monte Carlo blocks.
    persistent cached_M cached_N cached_nnz cached_checkEdges cached_varEdges cached_varIdx cached_numEdges;

    current_nnz = nnz(H);
    if isempty(cached_M) || cached_M ~= M || cached_N ~= N || cached_nnz ~= current_nnz
        [checkIdx, varIdx] = find(H);
        numEdges = length(checkIdx);

        % Vectorized O(E) edge grouping via accumarray
        edgeIdx = (1:numEdges)';
        varEdges = accumarray(varIdx, edgeIdx, [N 1], @(x) {x}, {zeros(0,1)});

        [sortCheck, sortOrder] = sort(checkIdx);
        checkEdges = accumarray(sortCheck, sortOrder, [M 1], @(x) {x}, {zeros(0,1)});

        cached_M = M;
        cached_N = N;
        cached_nnz = current_nnz;
        cached_checkEdges = checkEdges;
        cached_varEdges = varEdges;
        cached_varIdx = varIdx;
        cached_numEdges = numEdges;
    else
        checkEdges = cached_checkEdges;
        varEdges = cached_varEdges;
        varIdx = cached_varIdx;
        numEdges = cached_numEdges;
    end

    % Variable-to-check and check-to-variable message vectors (vectorized)
    v2c = llrChannel(varIdx);
    c2v = zeros(1, numEdges);

    decodedBits = double(llrChannel < 0);
    converged = false;
    iterUsed = 0;

    recordTelemetry = (nargout > 3);
    if recordTelemetry
        llr_history = zeros(maxIter, N);
        syndrome_history = zeros(maxIter, 1);
    end

    for it = 1:maxIter
        iterUsed = it;

        % --- check-node update (min-sum) ---
        for m = 1:M
            eIdx = checkEdges{m};
            deg = length(eIdx);
            if deg == 0, continue; end

            v = v2c(eIdx);
            s = sign(v + (v == 0) * 1e-12);
            prodS = prod(s);

            absV = abs(v);
            [minVal1, minIdx1] = min(absV);
            absV(minIdx1) = inf;
            minVal2 = min(absV);

            for k = 1:deg
                if k == minIdx1
                    c2v(eIdx(k)) = (prodS * s(k)) * minVal2;
                else
                    c2v(eIdx(k)) = (prodS * s(k)) * minVal1;
                end
            end
        end

        % --- variable-node update (extrinsic) ---
        totalLLR = llrChannel;
        for n = 1:N
            eIdx = varEdges{n};
            deg = length(eIdx);
            if deg == 0, continue; end

            inc = c2v(eIdx);
            sumInc = sum(inc);
            totalLLR(n) = llrChannel(n) + sumInc;

            for k = 1:deg
                v2c(eIdx(k)) = totalLLR(n) - inc(k);
            end
        end

        decodedBits = double(totalLLR < 0);

        syndrome = mod(H * decodedBits', 2);
        unsatCount = sum(syndrome ~= 0);

        if recordTelemetry
            llr_history(it, :) = totalLLR;
            syndrome_history(it) = unsatCount;
        end

        if unsatCount == 0
            converged = true;
            break;
        end
    end

    if recordTelemetry
        telemetry.llr_history = llr_history(1:iterUsed, :);
        telemetry.syndrome_history = syndrome_history(1:iterUsed);
        telemetry.converged = converged;
        telemetry.iterUsed = iterUsed;
    end
end
