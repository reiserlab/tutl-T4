function results = batch_analyze_1DRF(data_root, opts)
% BATCH_ANALYZE_1DRF  Process all 1DRF experiments and generate population plots.
%
%   RESULTS = BATCH_ANALYZE_1DRF(DATA_ROOT) processes every experiment
%   folder in DATA_ROOT, classifies cells as ON/OFF and control/ttl, then
%   generates population plots with median +/- MAD (default) or mean +/-
%   SEM (configurable via opts.stat_method).
%
%   RESULTS = BATCH_ANALYZE_1DRF(DATA_ROOT, OPTS) uses the options
%   structure to override default parameters.
%
%   INPUTS:
%     data_root - Path to the 1DRF data directory containing experiment
%                 folders (e.g. '/path/to/protocol2/data/1DRF')
%     opts      - (Optional) structure. Inherits all fields from
%                 analyze_single_experiment opts, plus:
%                   .on_threshold  - Frame value threshold for ON vs OFF
%                                    classification (default: 129)
%                   .save_figs     - Save population figures as PDF
%                                    (default: true)
%                   .save_dir      - Output directory for population PDFs
%                                    (default: <data_root>/population_results)
%                   .flash_baseline - Sample indices for bar flash baseline
%                                     (default: 1:5000)
%                   .stat_method   - 'median_mad' (default) or 'mean_sem'
%
%   OUTPUT:
%     results - Structure array (one element per cell) with fields:
%       .folder            - Experiment folder name
%       .date_str          - Experiment date
%       .strain            - Strain name from metadata
%       .frame             - Frame number from metadata
%       .is_on             - Logical, true if ON cell (Frame > on_threshold)
%       .is_ttl            - Logical, true if Strain contains 'ttl'
%       .group             - String: 'on_control', 'on_ttl', 'off_control',
%                            or 'off_ttl'
%       .max_v_aligned     - 16x2 PD-aligned [angles, responses] from
%                            find_PD_and_order_idx
%       .pd_direction      - Preferred direction in degrees
%       .pd_orientation    - Bar orientation at PD in degrees
%       .pd_flash_bl       - 11xN baseline-subtracted mean flash traces
%                            along PD-ND axis
%       .ortho_flash_bl    - 11xN baseline-subtracted mean flash traces
%                            along orthogonal axis
%       .pd_flash_baselines     - 11x1 per-position baseline values (mV)
%                                 used to recover absolute voltage from .pd_flash_bl
%       .ortho_flash_baselines  - 11x1 per-position baseline values (mV)
%       .centroid_m6_rounded        - integer M6 centroid (1..11) on PD axis
%       .ortho_centroid_m6_rounded  - integer M6 centroid (1..11) on ortho axis
%       .pd_flash_m6_aligned        - 11xN M6-aligned flash traces, PD axis
%       .ortho_flash_m6_aligned     - 11xN M6-aligned flash traces, ortho axis
%       .flash_pitch_px    - bar-flash position pitch in arena px (2, or 1 for
%                            the three Oct-2025 sessions); M6-aligned traces
%                            are always on the 2-px grid
%       .dsi_vector        - vector-sum direction selectivity index (manuscript)
%
%   FIGURES GENERATED:
%     For each of ON and OFF cells (6 figures total):
%       1. PD-aligned polar tuning curve (control=black, ttl=red, +/- spread)
%       2. 1x11 PD-ND bar flash (control=black, ttl=red, +/- spread)
%       3. 1x11 orthogonal bar flash (control=black, ttl=red, +/- spread)
%     Spread is MAD (default) or SEM, set via opts.stat_method.
%
%   ANALYSIS PIPELINE:
%     1. Discover experiment folders in data_root
%     2. For each experiment:
%        a. Load data and parse bar sweeps/flashes
%        b. Compute bar sweep responses, find PD, align to pi/2
%        c. Extract baseline-subtracted bar flash traces
%        d. Classify cell as ON/OFF and control/ttl
%     3. Group cells and generate population-averaged plots
%
%   EXAMPLE:
%     results = batch_analyze_1DRF('/path/to/protocol2/data/1DRF');
%
%   See also ANALYZE_SINGLE_EXPERIMENT, FIND_PD_AND_ORDER_IDX,
%            PLOT_POLAR_POPULATION, PLOT_FLASH_1X11_POPULATION,
%            VERIFY_LUT_DIRECTIONS, COMPUTE_BAR_SWEEP_RESPONSES,
%            FIND_PD_FROM_LUT
% ________________________________________________________________________

    %% Set defaults
    if nargin < 2, opts = struct(); end
    opts = set_batch_defaults(opts, data_root);

    %% Verify CircStat is on the path (used indirectly via find_PD_and_order_idx)
    if isempty(which('circ_vmpar'))
        error('batch_analyze_1DRF:CircStatMissing', ...
            ['CircStat toolbox not on MATLAB path. ' ...
             'Add CircStat2012a (or equivalent) to the path before calling. ' ...
             'See https://github.com/circstat/circstat-matlab.']);
    end

    %% Step 1: Load LUT (shared across experiments)
    S_lut = load(opts.lut_path, 'Tbl');
    Tbl = S_lut.Tbl;

    %% Step 2: Discover experiment folders
    d = dir(data_root);
    d = d([d.isdir]);
    d = d(~startsWith({d.name}, '.'));

    % Filter to folders that contain currentExp.mat
    valid = false(numel(d), 1);
    for i = 1:numel(d)
        valid(i) = isfile(fullfile(data_root, d(i).name, 'currentExp.mat'));
    end
    d = d(valid);
    n_exp = numel(d);

    fprintf('Found %d experiment folders in %s\n', n_exp, data_root);

    %% Step 3: Process each experiment
    results = struct([]);

    for exp_idx = 1:n_exp
        folder = d(exp_idx).name;
        exp_folder = fullfile(data_root, folder);
        fprintf('\n[%d/%d] Processing %s...\n', exp_idx, n_exp, folder);

        try
            r = process_single_cell(exp_folder, Tbl, opts);
            r.folder = folder;

            if isempty(results)
                results = r;
            else
                results(end + 1) = r; %#ok<AGROW>
            end

            fprintf('  -> %s | PD: %.0f° | %s | flash pitch %d px\n', ...
                r.group, r.pd_direction, r.strain, r.flash_pitch_px);

        catch ME
            fprintf('  ERROR: %s\n', ME.message);
            continue;
        end
    end

    fprintf('\n=== Processing Complete ===\n');
    fprintf('Total cells: %d\n', numel(results));

    % Count per group
    groups = {results.group};
    for g = ["on_control", "on_ttl", "off_control", "off_ttl"]
        fprintf('  %s: %d\n', g, sum(strcmp(groups, g)));
    end

    %% Step 4: Save results
    if opts.save_figs
        if ~isfolder(opts.save_dir)
            mkdir(opts.save_dir);
        end
        save_path = fullfile(opts.save_dir, 'batch_results.mat');
        save(save_path, 'results');
        fprintf('\nResults saved to: %s\n', save_path);
    end

    %% Step 5: Generate population plots
    fprintf('\nGenerating population plots...\n');

    for on_off_label = ["ON", "OFF"]
        if on_off_label == "ON"
            mask = [results.is_on];
        else
            mask = ~[results.is_on];
        end

        ctrl_mask = mask & ~[results.is_ttl];
        ttl_mask  = mask & [results.is_ttl];

        % --- Polar tuning curves ---
        aligned_ctrl = {results(ctrl_mask).max_v_aligned};
        aligned_ttl  = {results(ttl_mask).max_v_aligned};

        polar_opts.stat_method = opts.stat_method;
        polar_title = sprintf('%s Cells — PD-Aligned Polar Tuning', on_off_label);
        fig_polar = plot_polar_population(aligned_ctrl, aligned_ttl, ...
            polar_title, polar_opts);

        % --- PD-ND bar flash traces ---
        pd_traces_ctrl = {results(ctrl_mask).pd_flash_bl};
        pd_traces_ttl  = {results(ttl_mask).pd_flash_bl};

        flash_opts.y_limits    = opts.flash_ylim;
        flash_opts.plot_type   = 'pd_nd';
        flash_opts.stat_method = opts.stat_method;
        pd_flash_title = sprintf('%s Cells — PD-ND Bar Flash', on_off_label);
        fig_pd = plot_flash_1x11_population(...
            pd_traces_ctrl, pd_traces_ttl, pd_flash_title, flash_opts);

        % --- Orthogonal bar flash traces ---
        ortho_traces_ctrl = {results(ctrl_mask).ortho_flash_bl};
        ortho_traces_ttl  = {results(ttl_mask).ortho_flash_bl};

        flash_opts.fig_position = [50 100 1800 300];
        flash_opts.plot_type    = 'orthogonal';
        ortho_flash_title = sprintf('%s Cells — Orthogonal Bar Flash', on_off_label);
        fig_ortho = plot_flash_1x11_population(...
            ortho_traces_ctrl, ortho_traces_ttl, ortho_flash_title, flash_opts);

        % Save figures
        if opts.save_figs
            save_population_figures(opts.save_dir, on_off_label, ...
                opts.stat_method, fig_polar, fig_pd, fig_ortho);
        end
    end

    fprintf('\nDone.\n');

end


%% ========================= Core Processing ============================

function r = process_single_cell(exp_folder, Tbl, opts)
% PROCESS_SINGLE_CELL  Run the full per-cell pipeline and return results.

    % Save/restore working directory (load_protocol2_data uses cd)
    orig_dir = pwd;
    cleanup = onCleanup(@() cd(orig_dir));

    % Load data
    [date_str, ~, Log, ~, ~] = load_protocol2_data(exp_folder);

    f_data   = Log.ADC.Volts(1, :);
    v_data   = Log.ADC.Volts(2, :) * 10;
    median_v = median(v_data);

    ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
        'pattern_order', 'func_order', 'metadata');
    metadata = ce.metadata;

    % Classify cell
    r.date_str = date_str;
    r.strain   = metadata.Strain;
    r.frame    = metadata.Frame;
    r.is_on    = metadata.Frame > opts.on_threshold;
    r.is_ttl   = contains(metadata.Strain, 'ttl');

    if r.is_on && ~r.is_ttl
        r.group = 'on_control';
    elseif r.is_on && r.is_ttl
        r.group = 'on_ttl';
    elseif ~r.is_on && ~r.is_ttl
        r.group = 'off_control';
    else
        r.group = 'off_ttl';
    end

    % Parse bar sweep data and verify LUT
    bar_data = parse_bar_data(f_data, v_data);

    [lut_directions, lut_orientations, lut_patterns, lut_functions] = ...
        verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, opts.plot_order);

    % Compute bar sweep responses
    sweep_opts.baseline_range = opts.baseline_range;
    sweep_opts.stim_trim_end  = opts.stim_trim_end;
    sweep_opts.percentile     = opts.percentile;
    max_v = compute_bar_sweep_responses(bar_data, opts.plot_order, sweep_opts);

    % PD-align the tuning curve for population averaging
    % find_PD_and_order_idx expects data sorted by angle (0, 22.5, ..., 337.5)
    lut_dirs_ordered = lut_directions(opts.plot_order);
    [~, sort_idx] = sort(lut_dirs_ordered);
    max_v_sorted = max_v(sort_idx);
    max_v_polar = [max_v_sorted; max_v_sorted(1)];  % 17x1 for circular

    % median_voltage = 0 because max_v is already baseline-relative
    [d_aligned, ~, ~, ~, ~, ~, ~, ~] = find_PD_and_order_idx(max_v_polar, 0);
    r.max_v_aligned = d_aligned;

    % Direction-selectivity vector index (used by manuscript figures)
    [~, r.dsi_vector, ~, ~] = compute_bar_response_metrics(d_aligned);

    % Find PD and map to bar flash columns
    pd_info = find_pd_from_lut(max_v, lut_directions, lut_orientations, ...
        lut_patterns, lut_functions, opts.plot_order, Tbl, opts.pattern_offset);

    r.pd_direction  = pd_info.pd_direction;
    r.pd_orientation = pd_info.pd_orientation;

    % Parse bar flash data (prop_int = 0.5 matches the original protocol-2
    % gap_between_flashes setting; required since parse_bar_flash_data was
    % made parameterized).
    [data_slow_bf, ~, mean_slow_bf, ~] = parse_bar_flash_data(f_data, v_data, 0.5);

    % Extract baseline-subtracted mean flash traces (11 x N_timepoints)
    bl_samples = opts.flash_baseline;

    [r.pd_flash_bl,    r.pd_flash_baselines]    = extract_flash_traces(mean_slow_bf, ...
        pd_info.bar_flash_col,   pd_info.pos_order,       bl_samples);
    [r.ortho_flash_bl, r.ortho_flash_baselines] = extract_flash_traces(mean_slow_bf, ...
        pd_info.ortho_flash_col, pd_info.ortho_pos_order, bl_samples);

    % --- M6 alignment for manuscript figures (68%-area centroid -> row 6) ---
    % Bar-flash pitch is 2 px in all sessions except three Oct-23 sessions
    % recorded at 1 px. Every cell is placed on a common 2-px (2.5 deg)
    % grid: 1-px cells contribute every second stimulus to rows 4..8 only.
    r.flash_pitch_px = read_flash_pitch(exp_folder);
    step = 2 / r.flash_pitch_px;                 % 1 (2-px) or 2 (1-px)

    pd_peaks    = compute_pos_peak_amplitudes(mean_slow_bf, ...
        pd_info.bar_flash_col,   pd_info.pos_order,       bl_samples);
    ortho_peaks = compute_pos_peak_amplitudes(mean_slow_bf, ...
        pd_info.ortho_flash_col, pd_info.ortho_pos_order, bl_samples);

    m6_pd                       = compute_m6_centroid(max(pd_peaks, 0), 0.68);
    r.centroid_m6_rounded       = m6_pd.centroid_int;
    r.pd_flash_m6_aligned       = reindex_to_peak(r.pd_flash_bl, r.centroid_m6_rounded, step);

    m6_ortho                    = compute_m6_centroid(max(ortho_peaks, 0), 0.68);
    r.ortho_centroid_m6_rounded = m6_ortho.centroid_int;
    r.ortho_flash_m6_aligned    = reindex_to_peak(r.ortho_flash_bl, r.ortho_centroid_m6_rounded, step);

end


function pitch = read_flash_pitch(exp_folder)
% READ_FLASH_PITCH  Bar-flash position pitch (arena px) from the FLASHES pattern.
%
%   Frames 2 and 3 of Patterns/0011_*FLASHES*.mat are positions 1 and 2 of
%   the first (cardinal) bar orientation; the distance between their bar
%   centroids is the pitch (1 or 2 px). Frame 1 is the blank frame.

    d = dir(fullfile(exp_folder, 'Patterns', '*FLASHES*.mat'));
    S = load(fullfile(d(1).folder, d(1).name), 'pattern');
    c = zeros(1, 2);
    for k = 1:2
        fr = S.pattern.Pats(:, :, k + 1);
        [~, cols] = find(fr ~= mode(fr(:)));
        c(k) = mean(cols);
    end
    pitch = round(abs(c(2) - c(1)));
    assert(pitch == 1 || pitch == 2, ...
        'read_flash_pitch:UnexpectedPitch', 'Unexpected bar-flash pitch %g px', pitch);
end


function [traces, baselines] = extract_flash_traces(mean_slow_bf, flash_col, pos_order, bl_samples)
% EXTRACT_FLASH_TRACES  Get baseline-subtracted mean traces for one orientation.
%
%   [TRACES, BASELINES] = EXTRACT_FLASH_TRACES(...)
%   Returns an 11 x N matrix (positions x timepoints), ordered ND to PD.
%   BASELINES is an 11 x 1 vector of per-position baseline values (mV);
%   absolute voltage can be reconstructed as TRACES + BASELINES (with
%   broadcast).
%   All traces are truncated to the minimum length across positions to
%   handle minor length variations between flash stimuli.

    n_pos = 11;

    % First pass: find minimum trace length across all positions
    min_pts = Inf;
    for i = 1:n_pos
        ts = mean_slow_bf{i, flash_col};
        if ~isempty(ts)
            min_pts = min(min_pts, numel(ts));
        end
    end

    if isinf(min_pts)
        traces    = [];
        baselines = [];
        return;
    end

    traces    = NaN(n_pos, min_pts);
    baselines = NaN(n_pos, 1);

    for pos_idx = 1:n_pos
        flash_pos = pos_order(pos_idx);
        ts = mean_slow_bf{flash_pos, flash_col};
        if ~isempty(ts)
            ts_trunc = ts(1:min_pts);
            bl = mean(ts_trunc(bl_samples(bl_samples <= min_pts)));
            baselines(pos_idx) = bl;
            traces(pos_idx, :) = ts_trunc(:)' - bl;
        end
    end

end


function peak_amplitudes = compute_pos_peak_amplitudes(mean_slow_bf, flash_col, pos_order, bl_samples)
% COMPUTE_POS_PEAK_AMPLITUDES  Per-position 99.5%ile depolarization (1x11).
%
%   Used as the input to COMPUTE_M6_CENTROID for M6-axis alignment.
%   Each output element is the 99.5th percentile of (response_window
%   - baseline) at that flash position, over response samples 5001:6551
%   (200 ms onset window + 75 ms tail).

    n_pos = 11;
    peak_amplitudes = zeros(1, n_pos);
    resp_start = 5001;
    resp_end   = 6551;  % 5801 stim offset + 750 sample tail

    for pos_idx = 1:n_pos
        flash_pos = pos_order(pos_idx);
        ts = mean_slow_bf{flash_pos, flash_col};
        if isempty(ts), continue; end
        bl  = mean(ts(bl_samples(bl_samples <= numel(ts))));
        win = resp_start : min(resp_end, numel(ts));
        if ~isempty(win)
            peak_amplitudes(pos_idx) = prctile(ts(win) - bl, 99.5);
        end
    end

end


%% ========================= Default Options ============================

function opts = set_batch_defaults(opts, data_root)
% SET_BATCH_DEFAULTS  Fill in default values for batch analysis options.

    % LUT path: same directory as this script
    if ~isfield(opts, 'lut_path')
        script_dir = fileparts(mfilename('fullpath'));
        opts.lut_path = fullfile(script_dir, 'bar_lut.mat');
    end
    if ~isfield(opts, 'plot_order')
        opts.plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
    end
    if ~isfield(opts, 'baseline_range')
        opts.baseline_range = [1000 9000];
    end
    if ~isfield(opts, 'stim_trim_end')
        opts.stim_trim_end = 7000;
    end
    if ~isfield(opts, 'percentile')
        opts.percentile = 98;
    end
    if ~isfield(opts, 'flash_baseline')
        opts.flash_baseline = 1:5000;
    end
    if ~isfield(opts, 'flash_ylim')
        opts.flash_ylim = [-15 35];
    end
    if ~isfield(opts, 'pattern_offset')
        opts.pattern_offset = 2;
    end
    if ~isfield(opts, 'on_threshold')
        opts.on_threshold = 129;
    end
    if ~isfield(opts, 'save_figs')
        opts.save_figs = true;
    end
    if ~isfield(opts, 'save_dir')
        opts.save_dir = fullfile(data_root, 'population_results');
    end
    if ~isfield(opts, 'stat_method')
        opts.stat_method = 'median_mad'; %  'mean_sem'
    end

end


%% ========================= Figure Saving ==============================

function save_population_figures(save_dir, on_off_label, stat_method, ...
    fig_polar, fig_pd, fig_ortho)
% SAVE_POPULATION_FIGURES  Export population figures as 300 dpi PDFs.
%   Filenames include the stat_method (e.g. 'median_mad' or 'mean_sem').

    if ~isfolder(save_dir)
        mkdir(save_dir);
    end

    export_opts = {'ContentType', 'image', 'Resolution', 300};
    prefix = lower(char(on_off_label));

    exportgraphics(fig_polar, ...
        fullfile(save_dir, sprintf('%s_polar_population_%s.pdf', prefix, stat_method)), ...
        export_opts{:});
    exportgraphics(fig_pd, ...
        fullfile(save_dir, sprintf('%s_pd_nd_flash_population_%s.pdf', prefix, stat_method)), ...
        export_opts{:});
    exportgraphics(fig_ortho, ...
        fullfile(save_dir, sprintf('%s_ortho_flash_population_%s.pdf', prefix, stat_method)), ...
        export_opts{:});

    fprintf('  %s figures saved to: %s\n', on_off_label, save_dir);

end
