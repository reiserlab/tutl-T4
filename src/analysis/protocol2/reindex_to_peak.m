function aligned = reindex_to_peak(traces_11xN, peak_pos, step)
% REINDEX_TO_PEAK  Shift an 11-row trace matrix so peak_pos -> row 6.
%
%   ALIGNED = REINDEX_TO_PEAK(TRACES_11xN, PEAK_POS) reindexes the rows
%   of TRACES_11xN (size 11-by-N) so that row PEAK_POS in the input
%   appears at row 6 (the central row) of the output.  Out-of-range
%   source rows become NaN.
%
%   ALIGNED = REINDEX_TO_PEAK(..., STEP) takes STEP input rows per output
%   row (default 1).  Use STEP = 2 for sessions recorded at 1-px bar-flash
%   pitch so that every output row spans 2 arena px like the 2-px sessions
%   (only rows 4..8 are then filled; the rest are NaN).
%
%   Used by BATCH_ANALYZE_1DRF and the manuscript figure scripts to
%   align per-cell flash traces by the M6 spatial centroid.

    if nargin < 3, step = 1; end
    n_pos = 11;
    center = 6;
    N = size(traces_11xN, 2);
    aligned = NaN(n_pos, N);
    for i = 1:n_pos
        src = peak_pos + (i - center) * step;
        if src >= 1 && src <= n_pos
            aligned(i, :) = traces_11xN(src, :);
        end
    end

end
