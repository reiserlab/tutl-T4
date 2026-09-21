function fig = generate_manuscript_fig_ds(speed_dps, opts)
% GENERATE_MANUSCRIPT_FIG_DS  Bar-sweep ring/polar/DSI/AR/Vm sub-figure.
%
%   FIG = GENERATE_MANUSCRIPT_FIG_DS(SPEED_DPS) builds the bottom-half
%   bar-sweep figure for the manuscript at the requested bar speed:
%
%       SPEED_DPS = 56  (main figure, 56 dps)
%       SPEED_DPS = 28  (supplementary figure, 28 dps)
%
%   Combined 48-cell dataset (23 early-batch + 25 late-batch).
%   Gaussian-convolution alignment (99.9th percentile).
%
%   Panels:
%       A: T4 (ON) ring-of-traces + polar (ctrl vs tutl-)
%       B: T5 (OFF) ring-of-traces + polar (ctrl vs tutl-)
%       (C/E for 56 dps, C for 28 dps): DSI + aspect-ratio box plots
%       (E/G for 56 dps, E for 28 dps): Vm pre vs during stim
%
%   FIG = GENERATE_MANUSCRIPT_FIG_DS(SPEED_DPS, OPTS) accepts an options
%   struct with fields:
%       .data_root   - data_root (default '/Users/reiserm/Documents/ttl_1DRF')
%       .skip_export - true to suppress PDF/PNG export (default false)
%       .stamp_path  - if non-empty, save pre-plot variables to this
%                      .mat file and return without plotting.
%
%   See also GENERATE_MANUSCRIPT_FIG_EF, GENERATE_MANUSCRIPT_FIG.

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'data_root'),   opts.data_root   = '/Users/reiserm/Documents/ttl_1DRF'; end
if ~isfield(opts, 'skip_export'), opts.skip_export = false; end
if ~isfield(opts, 'stamp_path'),  opts.stamp_path  = ''; end

data_root = opts.data_root;
out_dir   = fullfile(data_root, 'manuscript_figures');
if ~isfolder(out_dir) && ~opts.skip_export, mkdir(out_dir); end

%% ========================= User-configurable flags ==========================
SUBTRACT_BASELINE = false;            % v3: absolute voltage ring traces
POLAR_METRIC_MODE = 'baseline_subtracted';

% --- Speed dispatch ---
SPEED_DPS = speed_dps;
switch SPEED_DPS
    case 56
        ROW_OFFSET_LATE = 16;
        PRE_BF_SPEED    = 'medium';
        FUNC_PAIR       = [5 6];
        FWHM_SAMPLES    = 2500;
        % Panel letters for combined main figure (E,F,G,H placement in this fig)
        PANEL_LETTER_T4_RING = 'E';
        PANEL_LETTER_T5_RING = 'F';
        PANEL_LETTER_BOX     = 'G';
        PANEL_LETTER_VM      = 'H';
        SHOW_RING_ZOOM       = true;   % apply RING_XLIM zoom + thicker traces
        SHOW_POLAR_DIAG      = true;   % print polar_mean max diagnostics
        POLAR_RPAD_FACTOR    = 1.10;   % overlay slightly larger
        RING_SCALEBAR_T      = 500;    % scale bar 0.5 s
    case 28
        ROW_OFFSET_LATE = 0;
        PRE_BF_SPEED    = 'slow';
        FUNC_PAIR       = [3 4];
        FWHM_SAMPLES    = 5000;
        PANEL_LETTER_T4_RING = 'C';
        PANEL_LETTER_T5_RING = 'D';
        PANEL_LETTER_BOX     = 'E';
        PANEL_LETTER_VM      = 'F';
        SHOW_RING_ZOOM       = false;
        SHOW_POLAR_DIAG      = false;
        POLAR_RPAD_FACTOR    = 1.15;
        RING_SCALEBAR_T      = 1000;   % scale bar 1.0 s
    otherwise
        error('generate_manuscript_fig_ds:BadSpeed', ...
            'speed_dps must be 28 or 56, got %d', SPEED_DPS);
end

% Processing constants (must match explore_population_ring_of_traces.m)
BL_START       = 1000;
BL_END         = 9000;
STIM_START     = 9000;
STIM_TRIM_END  = 7000;
DISPLAY_HALF   = 14000;
DS_FACTOR      = 10;

% --- Gaussian alignment ---
POLAR_PERCENTILE = 99.9;
GAUSS_SIGMA  = FWHM_SAMPLES / (2 * sqrt(2 * log(2)));
GAUSS_HALF   = ceil(3 * GAUSS_SIGMA);
gauss_kernel = exp(-(-GAUSS_HALF:GAUSS_HALF).^2 / (2 * GAUSS_SIGMA^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);
SEARCH_HALF  = 6000;   % +/-600 ms constraint window

plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

REF_SAMPLE  = DISPLAY_HALF + 1;
DISPLAY_LEN = 2 * DISPLAY_HALF + 1;

% LUT
lut_path = which('bar_lut.mat');
if isempty(lut_path)
    lut_path = fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'analysis', ...
        'protocol2', 'bar_lut.mat');
end
S_lut = load(lut_path, 'Tbl');
Tbl = S_lut.Tbl;

%% Data loading (with cache)
cache_tag = sprintf('gauss_999_%ddps', SPEED_DPS);
if ~SUBTRACT_BASELINE
    cache_tag = [cache_tag '_abs'];
end
% CACHE_VERSION guards against silently reusing a cache built by older code
% (bump whenever extraction, alignment or per-cell labelling changes).
CACHE_VERSION = 'v2';   % v2: canonical peak angle labels (2026-09-20)
cache_file = fullfile(data_root, 'population_results', ...
    sprintf('ring_of_traces_cache_%s_%s.mat', cache_tag, CACHE_VERSION));

if isfile(cache_file)
    fprintf('Loading cached ring-of-traces data (%s)...\n', cache_tag);
    C = load(cache_file);
    all_cells         = C.all_cells;
    display_len_ds    = C.display_len_ds;
    ref_sample_ds     = C.ref_sample_ds;
    groups_final      = C.groups_final;
    pd_aligned_angles = C.pd_aligned_angles;
    fprintf('  Loaded %d cells from cache.\n', numel(all_cells));
else
    fprintf('No cache found — extracting from raw data...\n');

    % Load batch results for validated PD directions
    S_late_batch = load(fullfile(data_root, 'population_results', ...
        'batch_results.mat'), 'results');
    late_batch = S_late_batch.results;
    S_early_batch = load(fullfile(data_root, 'pre-bar-flash', ...
        'population_results', 'batch_results_pre_bf.mat'), 'results');
    early_batch = S_early_batch.results;

    % --- Discover experiments ---
    % Late dataset
    late_root = data_root;
    d_late = dir(late_root);
    d_late = d_late([d_late.isdir]);
    d_late = d_late(~startsWith({d_late.name}, '.'));
    late_valid = false(numel(d_late), 1);
    for i = 1:numel(d_late)
        late_valid(i) = isfile(fullfile(late_root, d_late(i).name, 'currentExp.mat'));
    end
    d_late = d_late(late_valid);

    % Early dataset
    early_root = fullfile(data_root, 'pre-bar-flash');
    early_list = discover_early_experiments(early_root);

    n_late  = numel(d_late);
    n_early = numel(early_list);
    fprintf('  Found %d late + %d early experiments\n', n_late, n_early);

    % --- Extract traces ---
    all_cells = struct([]);
    cell_idx = 0;

    % Process late dataset
    fprintf('  Loading late dataset...\n');
    for exp_idx = 1:n_late
        folder_name = d_late(exp_idx).name;
        exp_folder = fullfile(late_root, folder_name);
        try
            orig_dir = pwd;
            cleanup = onCleanup(@() cd(orig_dir));

            [date_str, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
            f_data = Log.ADC.Volts(1, :);
            v_data = Log.ADC.Volts(2, :) * 10;

            ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
                'pattern_order', 'func_order', 'metadata');
            bar_data = parse_bar_data(f_data, v_data);

            is_on  = ce.metadata.Frame > 129;
            is_ttl = contains(ce.metadata.Strain, 'ttl');

            [lut_directions, ~, ~, ~] = ...
                verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);

            batch_entry = find_batch_entry(late_batch, folder_name);
            if isempty(batch_entry), continue; end

            [traces_16, pd_shift, lut_dirs_ordered] = extract_and_align_traces( ...
                bar_data, plot_order, lut_directions, ...
                batch_entry.pd_direction, ROW_OFFSET_LATE);

            cell_idx = cell_idx + 1;
            all_cells(cell_idx).traces_aligned = traces_16;
            all_cells(cell_idx).pd_shift       = pd_shift;
            all_cells(cell_idx).is_on          = is_on;
            all_cells(cell_idx).is_ttl         = is_ttl;
            all_cells(cell_idx).date_str       = [date_str '_' folder_name(end-4:end)];
            all_cells(cell_idx).batch          = 'late';
            all_cells(cell_idx).lut_dirs       = lut_dirs_ordered;
            all_cells(cell_idx).max_v_aligned  = batch_entry.max_v_aligned;
            all_cells(cell_idx).dsi            = batch_entry.dsi_vector;

            fprintf('    [%d/%d] %s OK\n', exp_idx, n_late, folder_name);
        catch ME
            fprintf('    [%d/%d] %s ERROR: %s\n', exp_idx, n_late, folder_name, ME.message);
        end
    end

    % Process early dataset
    fprintf('  Loading early dataset...\n');
    for exp_idx = 1:n_early
        ei = early_list(exp_idx);
        try
            orig_dir = pwd;
            cleanup = onCleanup(@() cd(orig_dir));

            [date_str, ~, Log, ~, ~] = load_protocol2_data(ei.folder);
            f_data = Log.ADC.Volts(1, :);
            v_data = Log.ADC.Volts(2, :) * 10;

            ce = load(fullfile(ei.folder, 'currentExp.mat'), ...
                'pattern_order', 'func_order', 'metadata');
            bar_data = parse_bar_data_pre_bf(f_data, v_data, PRE_BF_SPEED);

            is_on  = strcmpi(ei.cell_type, 'ON');
            is_ttl = strcmpi(ei.treatment, 'ttl');
            is_off = ~is_on;

            [bar_data, ~] = correct_off_polarity_swap( ...
                bar_data, ce.pattern_order, ce.func_order, ...
                ei.date_str, is_off, FUNC_PAIR);

            [lut_directions, ~, ~, ~] = ...
                verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);

            batch_entry = find_batch_entry(early_batch, ei.folder);
            if isempty(batch_entry), continue; end

            [traces_16, pd_shift, lut_dirs_ordered] = extract_and_align_traces( ...
                bar_data, plot_order, lut_directions, ...
                batch_entry.pd_direction, 0);  % no offset — pre_bf returns 16 rows

            cell_idx = cell_idx + 1;
            all_cells(cell_idx).traces_aligned = traces_16;
            all_cells(cell_idx).pd_shift       = pd_shift;
            all_cells(cell_idx).is_on          = is_on;
            all_cells(cell_idx).is_ttl         = is_ttl;
            all_cells(cell_idx).date_str       = ei.date_str;
            all_cells(cell_idx).batch          = 'early';
            all_cells(cell_idx).lut_dirs       = lut_dirs_ordered;
            all_cells(cell_idx).max_v_aligned  = batch_entry.max_v_aligned;
            all_cells(cell_idx).dsi            = batch_entry.dsi_vector;

            fprintf('    [%d/%d] %s OK\n', exp_idx, n_early, ei.date_str);
        catch ME
            fprintf('    [%d/%d] %s ERROR: %s\n', exp_idx, n_early, ei.date_str, ME.message);
        end
    end

    n_cells = numel(all_cells);
    fprintf('  Loaded %d cells\n', n_cells);

    % --- Baseline subtraction (always performed for alignment) ---
    % Store per-trace baselines so we can add them back in absolute mode.
    cell_baselines = cell(n_cells, 16);
    for ci = 1:n_cells
        for di = 1:16
            tr = all_cells(ci).traces_aligned{di};
            if isempty(tr), continue; end
            bl = mean(tr(BL_START:min(BL_END, numel(tr))));
            cell_baselines{ci, di} = bl;
            all_cells(ci).traces_aligned{di} = tr - bl;
        end
    end

    % --- Recompute peak amplitudes using current POLAR_PERCENTILE ---
    % Angle labels are the canonical PD-aligned frame of traces_aligned (row 5 = PD,
    % rows 1/9 orthogonal, row 13 = ND; see extract_and_align_traces). The batch
    % max_v_aligned(:,1) column must NOT be used here: for PD = 292.5 deg its 0-deg
    % label wraps to 2*pi in find_PD_and_order_idx and sorts last, which made
    % compute_ar read rows 4/1/8 instead of 5/1/9 for six cells (fixed 2026-09-20).
    for ci = 1:n_cells
        angles = deg2rad((0:15)' * 22.5);
        peak_amps = NaN(16, 1);
        for di = 1:16
            tr = all_cells(ci).traces_aligned{di};
            if isempty(tr), continue; end
            stim_end = max(1, numel(tr) - STIM_TRIM_END);
            if STIM_START > stim_end, continue; end
            d_stim = tr(STIM_START:stim_end);
            peak_amps(di) = prctile(d_stim, POLAR_PERCENTILE);
        end
        all_cells(ci).peak_amps = [angles, peak_amps];
    end

    % =================================================================
    % Gaussian convolution alignment
    % =================================================================
    fprintf('  Using Gaussian convolution alignment (FWHM=%d, pctile=%.1f)...\n', ...
            FWHM_SAMPLES, POLAR_PERCENTILE);

        time_to_max   = NaN(n_cells, 16);
        align_valid   = false(n_cells, 16);
        aligned_traces = cell(n_cells, 16);

        for ci = 1:n_cells
            for di = 1:16
                tr = all_cells(ci).traces_aligned{di};
                if isempty(tr), continue; end
                bl_std = std(tr(BL_START:min(BL_END, numel(tr))));
                tr_conv = conv(tr, gauss_kernel, 'same');
                stim_end = max(1, numel(tr_conv) - STIM_TRIM_END);
                if STIM_START > stim_end, continue; end
                % Rough peak in full stim window
                [~, rough_idx] = max(tr_conv(STIM_START:stim_end));
                rough_time = STIM_START + rough_idx - 1;
                % Constrained fine peak within +/-SEARCH_HALF
                search_lo = max(1, rough_time - SEARCH_HALF);
                search_hi = min(numel(tr_conv), rough_time + SEARCH_HALF);
                [max_val, fine_idx] = max(tr_conv(search_lo:search_hi));
                peak_time = search_lo + fine_idx - 1;
                time_to_max(ci, di) = peak_time;
                if max_val > bl_std
                    align_valid(ci, di) = true;
                end
            end

            % Shift traces to align peaks
            for di = 1:16
                if align_valid(ci, di)
                    shift_time = time_to_max(ci, di);
                else
                    shift_time = borrow_nearest_shift(time_to_max(ci, :), align_valid(ci, :), di);
                    if isnan(shift_time)
                        aligned_traces{ci, di} = NaN(DISPLAY_LEN, 1);
                        continue;
                    end
                end
                tr = all_cells(ci).traces_aligned{di};
                if isempty(tr)
                    aligned_traces{ci, di} = NaN(DISPLAY_LEN, 1);
                    continue;
                end
                aligned_traces{ci, di} = shift_and_window(tr, shift_time, REF_SAMPLE, DISPLAY_LEN);
            end
        end

        % Downsample
        aligned_traces_ds = cell(size(aligned_traces));
        for ci = 1:n_cells
            for di = 1:16
                aligned_traces_ds{ci, di} = blockavg_downsample(aligned_traces{ci, di}, DS_FACTOR);
            end
        end
        display_len_ds = ceil(DISPLAY_LEN / DS_FACTOR);
        ref_sample_ds  = ceil(REF_SAMPLE / DS_FACTOR);

        % --- Restore absolute voltage if requested ---
        if ~SUBTRACT_BASELINE
            fprintf('  Adding baselines back to aligned traces (absolute mode)...\n');
            for ci = 1:n_cells
                for di = 1:16
                    bl = cell_baselines{ci, di};
                    if ~isempty(bl) && ~isempty(aligned_traces_ds{ci, di})
                        aligned_traces_ds{ci, di} = aligned_traces_ds{ci, di} + bl;
                    end
                end
            end

            if strcmp(POLAR_METRIC_MODE, 'per_cell_normalized')
                % Option B: per-cell normalization — subtract each cell's
                % mean peak amplitude across all 16 directions.
                % Shifts tuning curves to zero-centered while preserving shape.
                fprintf('  Normalizing peak_amps per-cell (subtracting mean across directions)...\n');
                for ci = 1:n_cells
                    pa = all_cells(ci).peak_amps;
                    if isempty(pa) || size(pa, 2) < 2, continue; end
                    cell_mean = mean(pa(:, 2), 'omitnan');
                    all_cells(ci).peak_amps(:, 2) = pa(:, 2) - cell_mean;
                end
            end
            % Option A ('baseline_subtracted'): peak_amps stay as-is
            % (already baseline-subtracted from the alignment step)
        end

        % Group cells and compute population stats (no xcorr step)
        on_ctrl_idx  = find([all_cells.is_on] & ~[all_cells.is_ttl]);
        on_ttl_idx   = find([all_cells.is_on] &  [all_cells.is_ttl]);
        off_ctrl_idx = find(~[all_cells.is_on] & ~[all_cells.is_ttl]);
        off_ttl_idx  = find(~[all_cells.is_on] &  [all_cells.is_ttl]);

        pd_aligned_angles = (0:15)' * 22.5;

        groups_final = struct( ...
            'name',     {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'}, ...
            'indices',  {on_ctrl_idx, on_ttl_idx, off_ctrl_idx, off_ttl_idx}, ...
            'color',    {[0 0 0], [1 0 0], [0 0 0], [1 0 0]} ...
        );

        for g = 1:4
            idx = groups_final(g).indices;
            n_g = numel(idx);
            groups_final(g).n = n_g;
            groups_final(g).mean_traces = cell(16, 1);
            groups_final(g).sem_traces  = cell(16, 1);

            polar_amps = NaN(16, n_g);
            for k = 1:n_g
                pa = all_cells(idx(k)).peak_amps;
                if ~isempty(pa) && size(pa, 2) == 2
                    polar_amps(:, k) = pa(:, 2);
                end
            end
            groups_final(g).polar_mean = mean(polar_amps, 2, 'omitnan');
            groups_final(g).polar_sem  = std(polar_amps, 0, 2, 'omitnan') / sqrt(n_g);

            for di = 1:16
                traces_mat = NaN(display_len_ds, n_g);
                for k = 1:n_g
                    tr = aligned_traces_ds{idx(k), di};
                    len = min(numel(tr), display_len_ds);
                    traces_mat(1:len, k) = tr(1:len);
                end
                groups_final(g).mean_traces{di} = mean(traces_mat, 2, 'omitnan');
                groups_final(g).sem_traces{di}  = std(traces_mat, 0, 2, 'omitnan') / sqrt(n_g);
            end
        end

    % (xcorr alignment path removed — Gaussian alignment only)

    % Save cache
    fprintf('  Saving cache to %s...\n', cache_file);
    save(cache_file, 'all_cells', ...
        'display_len_ds', 'ref_sample_ds', 'groups_final', 'pd_aligned_angles', ...
        '-v7.3');
    fprintf('  Cache saved.\n');
end

%% --- Apply per-cell normalization if requested (runs after cache load too) ---
if ~SUBTRACT_BASELINE && strcmp(POLAR_METRIC_MODE, 'per_cell_normalized')
    fprintf('Applying per-cell normalization to peak_amps for polar/DSI/AR...\n');
    for ci = 1:numel(all_cells)
        pa = all_cells(ci).peak_amps;
        if isempty(pa) || size(pa, 2) < 2, continue; end
        cell_mean = mean(pa(:, 2), 'omitnan');
        all_cells(ci).peak_amps(:, 2) = pa(:, 2) - cell_mean;
    end
    % Also update groups_final polar stats to reflect normalization
    for g = 1:numel(groups_final)
        idx = groups_final(g).indices;
        polar_amps = NaN(16, numel(idx));
        for k = 1:numel(idx)
            pa = all_cells(idx(k)).peak_amps;
            if ~isempty(pa) && size(pa, 2) == 2
                polar_amps(:, k) = pa(:, 2);
            end
        end
        groups_final(g).polar_mean = mean(polar_amps, 2, 'omitnan');
        groups_final(g).polar_sem  = std(polar_amps, 0, 2, 'omitnan') / sqrt(numel(idx));
    end
end

%% Build Panel C data from all_cells (uses current POLAR_PERCENTILE)
fprintf('Building Panel C statistics from all_cells (%.1fth pctile)...\n', POLAR_PERCENTILE);
combined = struct([]);
for ci = 1:numel(all_cells)
    c = all_cells(ci);
    s.is_on  = c.is_on;
    s.is_ttl = c.is_ttl;
    d = c.peak_amps;  % 16x2 [angles, peaks]
    s.dsi_pdnd = (d(5,2) - d(13,2)) / (d(5,2) + d(13,2));
    s.max_v_aligned = d;  % needed by compute_ar
    combined = [combined, s]; %#ok<AGROW>
end
fprintf('  Combined: %d cells for DSI + AR\n', numel(combined));

%% Build Panel D data: per-cell sweep voltage (pre-stim + during-stim)
fprintf('Extracting sweep voltages for Panel D...\n');
n_cells_total = numel(all_cells);
voltage_pre  = NaN(n_cells_total, 1);
voltage_stim = NaN(n_cells_total, 1);
orig_dir_v = pwd;

% Build folder lookup for early cells (date_str -> full path)
early_root_v = fullfile(data_root, 'pre-bar-flash');
early_list_v = discover_early_experiments(early_root_v);
early_folder_map = containers.Map();
for ei_idx = 1:numel(early_list_v)
    early_folder_map(early_list_v(ei_idx).date_str) = early_list_v(ei_idx).folder;
end

for ci = 1:n_cells_total
    c = all_cells(ci);
    try
        % Resolve experiment folder
        if strcmp(c.batch, 'late')
            folder_v = fullfile(data_root, c.date_str);
        else
            folder_v = early_folder_map(c.date_str);
        end
        [~, ~, Log_v, ~, ~] = load_protocol2_data(folder_v);
        cd(orig_dir_v);

        f_data_v = Log_v.ADC.Volts(1, :);
        v_data_v = Log_v.ADC.Volts(2, :) * 10;

        % Use the sweeps of THIS figure's speed for every cell (2026-09-20 fix;
        % previously late cells pooled 28 + 56 dps while early cells used 28 only).
        if strcmp(c.batch, 'late')
            bar_data_v = parse_bar_data(f_data_v, v_data_v);
            use_rows_v = ROW_OFFSET_LATE + (1:16);          % rows 1-16 = 28 dps, 17-32 = 56 dps
        else
            bar_data_v = parse_bar_data_pre_bf(f_data_v, v_data_v, PRE_BF_SPEED);
            use_rows_v = 1:16;
        end
        use_rows_v = use_rows_v(use_rows_v <= size(bar_data_v, 1));
        % Both parsers pad 9000 samples (900 ms) before bar onset and after bar offset.
        % Pre-sweep window: 800 ms before onset. During-sweep window: onset to offset
        % exactly (2026-09-20 fix; previously 7000 were trimmed, leaving 200 ms of blank
        % screen inside the "during sweep" window).
        pre_v = []; stim_v = [];
        for d = use_rows_v
            tr = bar_data_v{d, 4};
            if isempty(tr), continue; end
            tr = tr(:)'; tlen = numel(tr);
            pre_end = min(9000, tlen);
            pre_v = [pre_v, tr(1000:pre_end)]; %#ok<AGROW>
            stim_end = tlen - 9000;
            if stim_end > 9001
                stim_v = [stim_v, tr(9001:stim_end)]; %#ok<AGROW>
            end
        end
        voltage_pre(ci)  = median(pre_v);
        voltage_stim(ci) = median(stim_v);
    catch ME
        cd(orig_dir_v);
        fprintf('  [%2d/%d] voltage ERROR: %s\n', ci, n_cells_total, ME.message);
    end
end
fprintf('  Voltage extracted: %d / %d cells\n', sum(~isnan(voltage_pre)), n_cells_total);
% Diagnostics: group medians and rank-sum p (T4 ctrl vs tutl, T5 ctrl vs tutl, T4 ctrl vs T5 ctrl)
vg = {[all_cells.is_on] & ~[all_cells.is_ttl], [all_cells.is_on] & [all_cells.is_ttl], ...
      ~[all_cells.is_on] & ~[all_cells.is_ttl], ~[all_cells.is_on] & [all_cells.is_ttl]};
for vv = {voltage_pre, 'Vm pre-sweep   '; voltage_stim, 'Vm during sweep'}'
    x = vv{1}(:)';
    fprintf('  %s (%d dps): medians %s mV | p T4 %.4f, T5 %.4f, ctrl/ctrl %.4f\n', vv{2}, SPEED_DPS, ...
        mat2str(round(cellfun(@(m) median(x(m), 'omitnan'), vg), 2)), ...
        ranksum(x(vg{1}), x(vg{2})), ranksum(x(vg{3}), x(vg{4})), ranksum(x(vg{1}), x(vg{3})));
end

% Build voltage struct array matching draw_boxplot_panel expectations
voltage_combined = struct([]);
for ci = 1:n_cells_total
    c = all_cells(ci);
    s.is_on  = c.is_on;
    s.is_ttl = c.is_ttl;
    s.voltage_pre  = voltage_pre(ci);
    s.voltage_stim = voltage_stim(ci);
    voltage_combined = [voltage_combined, s]; %#ok<AGROW>
end

%% ===================== STAMP GATE (validation only) ====================
if ~isempty(opts.stamp_path)
    stamp = struct();
    stamp.all_cells       = all_cells;
    stamp.groups_final    = groups_final;
    stamp.combined        = combined;
    stamp.voltage_pre     = voltage_pre;
    stamp.voltage_stim    = voltage_stim;
    stamp.voltage_combined = voltage_combined;
    stamp.SPEED_DPS       = SPEED_DPS;
    stamp.ROW_OFFSET_LATE = ROW_OFFSET_LATE;
    stamp.PRE_BF_SPEED    = PRE_BF_SPEED;
    stamp.FUNC_PAIR       = FUNC_PAIR;
    stamp.FWHM_SAMPLES    = FWHM_SAMPLES;
    save(opts.stamp_path, '-struct', 'stamp', '-v7.3');
    fprintf('STAMP wrote %s\n', opts.stamp_path);
    fig = [];
    return
end

%% Create composite figure
fprintf('\n=== Creating manuscript figure ===\n');

FIG_W_CM = 21;  FIG_H_CM = 9;
fig = figure('Units', 'centimeters', 'Position', [2 2 FIG_W_CM FIG_H_CM], ...
    'Color', 'w', 'PaperUnits', 'centimeters', ...
    'PaperSize', [FIG_W_CM FIG_H_CM], 'PaperPosition', [0 0 FIG_W_CM FIG_H_CM]);
set(fig, 'DefaultAxesFontName', 'Helvetica', 'DefaultTextFontName', 'Helvetica');

% Ring layout for manuscript (tightened vs exploratory)
% Aspect-ratio correction for circle on canvas.
FIG_AR = FIG_W_CM / FIG_H_CM;  % width / height
RING = struct( ...
    'radius_x',    0.126, ...     % ring radius in normalized x-coords (~20% up from 0.105)
    'radius_y',    0.126 * FIG_AR, ... % corrected for aspect ratio → true circle
    'subW',        0.048, ...     % subplot width  (~20% up from 0.040)
    'subH',        0.090, ...     % subplot height (~20% up from 0.075)
    'polar_scale', 0.65, ...      % polar plot fills more of the ring interior
    'font_label',  5);

t_ms = ((1:display_len_ds) - ref_sample_ds) * DS_FACTOR * 0.1;

% --- Compute shared y-limits across ON and OFF panels ---
all_mean_vals = [];
for g = [1 2 3 4]  % all 4 groups
    for di = 1:16
        all_mean_vals = [all_mean_vals; groups_final(g).mean_traces{di}]; %#ok<AGROW>
    end
end
y_range = [min(all_mean_vals, [], 'omitnan'), max(all_mean_vals, [], 'omitnan')];
y_pad = 0.15 * diff(y_range);
shared_ylim = [y_range(1) - y_pad, y_range(2) + y_pad];

% --- Per-direction Wilcoxon stats (ctrl vs tutl-) ---
[on_stats, on_stats_table] = compute_direction_stats( ...
    all_cells, groups_final(1).indices, groups_final(2).indices, pd_aligned_angles);
[off_stats, off_stats_table] = compute_direction_stats( ...
    all_cells, groups_final(3).indices, groups_final(4).indices, pd_aligned_angles);

% --- Ring trace display range and line width ---
if SHOW_RING_ZOOM
    % Central 50% of time axis for zoom; thicker lines for visibility (56 dps)
    t_half = (t_ms(end) - t_ms(1)) / 4;
    t_center = (t_ms(1) + t_ms(end)) / 2;
    RING_XLIM = [t_center - t_half, t_center + t_half];
    RING_TRACE_LW = 0.5;
else
    % 28 dps: use full time axis, default trace line width
    RING_XLIM = [];
    RING_TRACE_LW = [];
end

if SHOW_POLAR_DIAG
    fprintf('T4 ctrl polar_mean max: %.1f mV\n', max(groups_final(1).polar_mean));
    fprintf('T4 tutl  polar_mean max: %.1f mV\n', max(groups_final(2).polar_mean));
    fprintf('T5 ctrl polar_mean max: %.1f mV\n', max(groups_final(3).polar_mean));
    fprintf('T5 tutl  polar_mean max: %.1f mV\n', max(groups_final(4).polar_mean));
end

% --- Panel A: T4 (ON) ring-of-traces ---
panelA_center = [0.17, 0.47];
draw_ring_panel(fig, panelA_center, RING, ...
    groups_final(1), groups_final(2), ...   % ON ctrl, ON tutl-
    pd_aligned_angles, all_cells, ...
    t_ms, display_len_ds, shared_ylim, DS_FACTOR, on_stats, ...
    [0 0 0], [1 0 0], RING_XLIM, RING_TRACE_LW, RING_SCALEBAR_T, true, POLAR_RPAD_FACTOR);

% --- Panel B: T5 (OFF) ring-of-traces ---
panelB_center = [0.48, 0.47];
draw_ring_panel(fig, panelB_center, RING, ...
    groups_final(3), groups_final(4), ...   % OFF ctrl, OFF tutl-
    pd_aligned_angles, all_cells, ...
    t_ms, display_len_ds, shared_ylim, DS_FACTOR, off_stats, ...
    [0.4 0.4 0.4], [0.8 0.2 0.2], RING_XLIM, RING_TRACE_LW, [], false, POLAR_RPAD_FACTOR);

% --- Right-side panels: C (top row) and D (bottom row) ---
% Each row has 2 panels side by side spanning the right-side area.
R_LEFT = 0.71;   % left edge of right-side area (shifted right for narrower panels)
R_W    = 0.11;   % width per panel (narrower)
R_GAP  = 0.03;   % gap between the two panels in a row

% --- Panel C (top row): DSI + Aspect Ratio ---
C_Y = 0.56;  C_H = 0.28;

ax_dsi = axes(fig, 'Position', [R_LEFT, C_Y, R_W, C_H]);
draw_boxplot_panel(ax_dsi, combined, 'dsi_pdnd', {'DS index', '(PD-ND)/(PD+ND)'}, ...
    [0 0.85], [0 0.25 0.5 0.75]);

ax_ar = axes(fig, 'Position', [R_LEFT + R_W + R_GAP, C_Y, R_W, C_H]);
draw_ar_panel(ax_ar, combined);

% --- Panel D (bottom row): Pre-stim + During-sweep voltage ---
D_Y = 0.12;  D_H = 0.28;

% Compute shared voltage y-limits
v_all = [voltage_pre; voltage_stim];
v_lo = floor(min(v_all(~isnan(v_all)))) - 1;
v_hi = ceil(max(v_all(~isnan(v_all)))) + 2;
v_ticks = v_lo:2:v_hi;

ax_vpre = axes(fig, 'Position', [R_LEFT, D_Y, R_W, D_H]);
draw_boxplot_panel(ax_vpre, voltage_combined, 'voltage_pre', 'V_m (mV)', ...
    [v_lo v_hi], v_ticks);
title(ax_vpre, 'Pre-stim', 'FontSize', 7, 'FontWeight', 'normal');

ax_vstim = axes(fig, 'Position', [R_LEFT + R_W + R_GAP, D_Y, R_W, D_H]);
draw_boxplot_panel(ax_vstim, voltage_combined, 'voltage_stim', 'V_m (mV)', ...
    [v_lo v_hi], v_ticks);
title(ax_vstim, 'During sweep', 'FontSize', 7, 'FontWeight', 'normal');

% --- Panel labels ---
annotation(fig, 'textbox', [0.00, 0.90, 0.05, 0.08], ...
    'String', PANEL_LETTER_T4_RING, 'FontSize', 12, 'FontWeight', 'bold', ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
annotation(fig, 'textbox', [0.31, 0.90, 0.05, 0.08], ...
    'String', PANEL_LETTER_T5_RING, 'FontSize', 12, 'FontWeight', 'bold', ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
annotation(fig, 'textbox', [R_LEFT - 0.06, C_Y + C_H + 0.06, 0.05, 0.08], ...
    'String', PANEL_LETTER_BOX, 'FontSize', 12, 'FontWeight', 'bold', ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
annotation(fig, 'textbox', [R_LEFT - 0.06, D_Y + D_H + 0.06, 0.05, 0.08], ...
    'String', PANEL_LETTER_VM, 'FontSize', 12, 'FontWeight', 'bold', ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% --- Cell-type labels above ring panels ---
annotation(fig, 'textbox', [0.04, 0.92, 0.26, 0.06], ...
    'String', 'T4 (ON)', 'FontSize', 9, 'FontWeight', 'bold', ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
annotation(fig, 'textbox', [0.35, 0.92, 0.26, 0.06], ...
    'String', 'T5 (OFF)', 'FontSize', 9, 'FontWeight', 'bold', ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
% Speed label between B and C
annotation(fig, 'textbox', [0.28, 0.44, 0.10, 0.04], ...
    'String', sprintf('%d\\circ/s', SPEED_DPS), 'FontSize', 7, ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', 'Interpreter', 'tex', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

%% Print and save per-direction stats summary
fprintf('\n=== Per-direction Wilcoxon rank-sum: T4 (ON) ctrl vs tutl- ===\n');
print_direction_stats_table(on_stats_table);
fprintf('\n=== Per-direction Wilcoxon rank-sum: T5 (OFF) ctrl vs tutl- ===\n');
print_direction_stats_table(off_stats_table);

% Save stats table to file (only when actually exporting)
if ~opts.skip_export
    stats_file = fullfile(out_dir, sprintf('fig_ds_direction_stats_%ddps.txt', SPEED_DPS));
    fid = fopen(stats_file, 'w');
    if fid ~= -1
        fprintf(fid, 'Per-direction Wilcoxon rank-sum tests (ctrl vs tutl-)\n');
        fprintf(fid, 'FDR correction: Benjamini-Hochberg at q=0.05\n\n');
        fprintf(fid, '=== T4 (ON): ctrl (n=%d) vs tutl- (n=%d) ===\n', ...
            numel(groups_final(1).indices), numel(groups_final(2).indices));
        write_direction_stats_table(fid, on_stats_table);
        fprintf(fid, '\n=== T5 (OFF): ctrl (n=%d) vs tutl- (n=%d) ===\n', ...
            numel(groups_final(3).indices), numel(groups_final(4).indices));
        write_direction_stats_table(fid, off_stats_table);
        fclose(fid);
        fprintf('Saved: %s\n', stats_file);
    else
        warning('Could not open %s for writing; skipping stats sidecar.', stats_file);
    end
end

%% Export
ts = datestr(now, 'yyyymmdd_HHMM');
if ~SUBTRACT_BASELINE
    abs_tag = sprintf('_absolute_%s', POLAR_METRIC_MODE);
else
    abs_tag = '';
end
pdf_file = fullfile(out_dir, sprintf('fig_ds_panels_ABCD_%s%s_v3_%s.pdf', cache_tag, abs_tag, ts));
png_file = fullfile(out_dir, sprintf('fig_ds_panels_ABCD_%s%s_v3_%s.png', cache_tag, abs_tag, ts));

if ~opts.skip_export
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
    exportgraphics(fig, png_file, 'Resolution', 300);
    fprintf('Saved: %s\n', pdf_file);
    fprintf('Saved: %s\n', png_file);
end
fprintf('=== generate_manuscript_fig_ds(%d) done ===\n', SPEED_DPS);

end  % main function

%% ========================= Local Functions ===============================

function draw_ring_panel(fig, center, ring, g_ctrl, g_ttl, ...
    pd_aligned_angles, all_cells, t_ms, display_len_ds, y_lim, DS_FACTOR, dir_stats, ...
    ctrl_color, ttl_color, x_range, trace_lw, scalebar_t, show_pd_lines, rpad_factor)
% DRAW_RING_PANEL  Draw a ring-of-traces panel at the specified figure position.
%   show_pd_lines (optional): if true, draw PD/OD/ND reference lines to ctrl data values.
%   rpad_factor   (optional): polar overlay padding factor (1.10 main / 1.15 supp).

    if nargin < 12, dir_stats = []; end
    if nargin < 13 || isempty(ctrl_color), ctrl_color = [0 0 0]; end
    if nargin < 14 || isempty(ttl_color),  ttl_color  = [1 0 0]; end
    if nargin < 15 || isempty(x_range),    x_range = [t_ms(1), t_ms(end)]; end
    if nargin < 16 || isempty(trace_lw),   trace_lw = 0.5; end
    if nargin < 17 || isempty(scalebar_t), scalebar_t = 1000; end
    if nargin < 18 || isempty(show_pd_lines), show_pd_lines = false; end
    if nargin < 19 || isempty(rpad_factor), rpad_factor = 1.10; end

    % --- 16 radial timeseries subplots ---
    ax_scalebar = [];  % for scale bar (8 o'clock position = di 15)
    for di = 1:16
        angle_deg = pd_aligned_angles(di);
        angle_rad = deg2rad(180 - angle_deg);

        x_pos = center(1) + ring.radius_x * cos(angle_rad);
        y_pos = center(2) + ring.radius_y * sin(angle_rad);
        ax = axes(fig, 'Position', [x_pos - ring.subW/2, y_pos - ring.subH/2, ...
            ring.subW, ring.subH]); %#ok<LAXES>
        hold(ax, 'on');

        % -65 mV reference line (thin, full width, plotted first so data is on top)
        ref_mv = -65;
        if y_lim(1) <= ref_mv && y_lim(2) >= ref_mv
            plot(ax, [t_ms(1), t_ms(end)], [ref_mv ref_mv], 'k-', 'LineWidth', 0.3);
        end

        % Ctrl mean (no SEM band) — plotted on top of reference line
        plot(ax, t_ms, g_ctrl.mean_traces{di}, 'Color', ctrl_color, 'LineWidth', trace_lw);
        % tutl- mean (no SEM band)
        plot(ax, t_ms, g_ttl.mean_traces{di}, 'Color', ttl_color, 'LineWidth', trace_lw);

        ylim(ax, y_lim);
        xlim(ax, x_range);

        axis(ax, 'off');

        if di == 15, ax_scalebar = ax; end  % 8 o'clock position (315 deg)
    end

    % --- Central polar plot ---
    % Use radius_x for width and radius_y for height so polar fills the ring
    cs_x = ring.radius_x * ring.polar_scale;
    cs_y = ring.radius_y * ring.polar_scale;
    polar_pos = [center(1) - cs_x, center(2) - cs_y, 2*cs_x, 2*cs_y];

    theta = all_cells(g_ctrl.indices(1)).peak_amps(:, 1);

    % Lighter fill colors (alpha-blended with white)
    ctrl_fill = min(1, ctrl_color * 0.3 + 0.7);
    ttl_fill  = min(1, ttl_color * 0.3 + 0.7);

    polar_opts = struct( ...
        'n_ctrl', g_ctrl.n, 'n_ttl', g_ttl.n, ...
        'ctrl_label', 'control', 'ttl_label', '{\ittutl-}', ...
        'ctrl_line_color', ctrl_color, 'ctrl_fill_color', ctrl_fill, ...
        'ttl_line_color', ttl_color, 'ttl_fill_color', ttl_fill, ...
        'dir_stats', dir_stats, ...
        'rpad_factor', rpad_factor);
    [axPolar, axFill] = plot_polar_with_patch(polar_pos, theta, ...
        g_ctrl.polar_mean, g_ctrl.polar_sem, ...
        g_ttl.polar_mean, g_ttl.polar_sem, polar_opts);

    % --- Outward arrows at upper-right corner of each ring trace ---
    fig_sz = get(fig, 'PaperSize');  % [W, H] in cm
    W = fig_sz(1); H = fig_sz(2);
    arr_cm = 0.36;  % desired arrow length in cm
    for di = 1:16
        angle_deg = pd_aligned_angles(di);
        angle_rad = deg2rad(180 - angle_deg);
        x_pos = center(1) + ring.radius_x * cos(angle_rad);
        y_pos = center(2) + ring.radius_y * sin(angle_rad);
        tail_x = x_pos + ring.subW/2 - 0.003;
        tail_y = y_pos + ring.subH/2 - 0.003;
        % Direction in cm: cos/sin give unit direction in physical space
        dx_cm = cos(angle_rad);
        dy_cm = sin(angle_rad);
        phys_len = sqrt(dx_cm^2 + dy_cm^2);  % = 1
        % Convert desired cm displacement to normalized figure coords
        head_x = tail_x + (arr_cm * dx_cm / phys_len) / W;
        head_y = tail_y + (arr_cm * dy_cm / phys_len) / H;
        head_x = max(0.01, min(0.99, head_x));
        head_y = max(0.01, min(0.99, head_y));
        tail_x = max(0.01, min(0.99, tail_x));
        tail_y = max(0.01, min(0.99, tail_y));
        annotation(fig, 'arrow', [tail_x, head_x], [tail_y, head_y], ...
            'Color', [0 0 0], 'HeadWidth', 6, 'HeadLength', 4, ...
            'HeadStyle', 'plain', 'LineWidth', 1.5);
    end

    % --- PD / OD / ND reference lines on polar plot (T4 only) ---
    if show_pd_lines
        % Get actual ctrl data values at PD (90°), OD (0°), ND (270°)
        % theta is in radians; PD-aligned so PD = pi/2, OD = 0, ND = 3*pi/2
        ctrl_mean = g_ctrl.polar_mean;
        [~, pd_idx] = min(abs(theta - pi/2));
        [~, od_idx] = min(abs(theta - 0));
        [~, nd_idx] = min(abs(theta - 3*pi/2));
        r_pd = ctrl_mean(pd_idx);
        r_od = ctrl_mean(od_idx);
        r_nd = ctrl_mean(nd_idx);
        ref_lw = axPolar.LineWidth * 2;  % thicker than axis lines
        hold(axFill, 'on');
        % PD: up (y = +r), OD: right (x = +r), ND: down (y = -r)
        plot(axFill, [0 0],     [0 r_pd], 'k-', 'LineWidth', ref_lw);
        plot(axFill, [0 r_od],  [0 0],    'k-', 'LineWidth', ref_lw);
        plot(axFill, [0 0],     [0 -r_nd],'k-', 'LineWidth', ref_lw);
        % Labels near tips
        text(axFill, 0,          r_pd*1.1,  'PD', 'FontSize', 5, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        text(axFill, r_od*1.1,   0,         'OD', 'FontSize', 5, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left',   'VerticalAlignment', 'middle');
        text(axFill, 0,          -r_nd*1.1, 'ND', 'FontSize', 5, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end

    % --- Legend: upper-left of ring area ---
    leg_x = center(1) - ring.radius_x - 0.01;
    leg_y = center(2) + ring.radius_y + 0.01;
    annotation(fig, 'textbox', [leg_x, leg_y - 0.01, 0.08, 0.025], ...
        'String', 'control', 'FontSize', 7, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica', 'EdgeColor', 'none', 'Color', ctrl_color, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    annotation(fig, 'textbox', [leg_x, leg_y - 0.06, 0.08, 0.025], ...
        'String', '{\ittutl-}', 'Interpreter', 'tex', ...
        'FontSize', 7, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica', 'EdgeColor', 'none', 'Color', ttl_color, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

    % --- Scale bar at 8 o'clock subplot ---
    add_scale_bar_on_axes(ax_scalebar, t_ms, y_lim, ring.font_label, scalebar_t);
end


function draw_boxplot_panel(ax, results, field, y_label, y_limits, y_ticks)
% DRAW_BOXPLOT_PANEL  Box+dot plot with Wilcoxon brackets for 4 groups.
%   Groups: T4 ctrl, T4 tutl-, T5 ctrl, T5 tutl-
%   y_limits: [ymin ymax], y_ticks: vector of tick positions

    group_masks = {
        [results.is_on] & ~[results.is_ttl], ...   % T4 ctrl
        [results.is_on] &  [results.is_ttl], ...   % T4 tutl-
        ~[results.is_on] & ~[results.is_ttl], ...  % T5 ctrl
        ~[results.is_on] &  [results.is_ttl]       % T5 tutl-
    };

    % Colors: ctrl=black, tutl-=red, T4 darker, T5 lighter
    colors = [0 0 0; 1 0 0; 0.4 0.4 0.4; 0.8 0.2 0.2];
    % Lighter versions for box fill (alpha-blended with white)
    box_colors = min(1, colors * 0.3 + 0.7);

    n_groups = 4;
    all_vals = [];
    all_grp_idx = [];
    group_data = cell(1, n_groups);

    for g = 1:n_groups
        mask = group_masks{g};
        vals = [results(mask).(field)];
        vals = vals(~isnan(vals));
        group_data{g} = vals;
        ng = numel(vals);
        all_vals = [all_vals, vals]; %#ok<AGROW>
        all_grp_idx = [all_grp_idx, repmat(g, 1, ng)]; %#ok<AGROW>
    end

    if isempty(all_vals), return; end

    hold(ax, 'on');

    % Draw color-matched box plots per group
    for g = 1:n_groups
        idx = all_grp_idx == g;
        if sum(idx) < 2, continue; end
        bc = boxchart(ax, all_grp_idx(idx), all_vals(idx), ...
            'BoxFaceColor', box_colors(g,:), 'BoxEdgeColor', colors(g,:), ...
            'WhiskerLineColor', colors(g,:), 'MarkerStyle', 'none', ...
            'BoxWidth', 0.5, 'LineWidth', 1.0);
    end

    % Jittered dots — larger, more transparent, more spread
    for g = 1:n_groups
        idx = all_grp_idx == g;
        vals = all_vals(idx);
        x = g + 0.25 * (rand(size(vals)) - 0.5);
        scatter(ax, x, vals, 20, colors(g, :), 'filled', 'MarkerFaceAlpha', 0.5);
    end

    % Axis formatting
    ylim(ax, y_limits);
    ax.Clipping = 'off';  % prevent brackets from being clipped
    set(ax, 'YTick', y_ticks, 'TickDir', 'out', 'FontSize', 7);
    xlim(ax, [0.5, n_groups + 0.5]);
    box(ax, 'off');
    % Floating axes: hide defaults, draw manual lines between first/last ticks
    ax.XColor = 'none';
    ax.YColor = 'none';
    % Manual x-axis line and ticks
    line(ax, [1 n_groups], [y_limits(1) y_limits(1)], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
    tick_arm_y = 0.015 * diff(y_limits);
    for g = 1:n_groups
        line(ax, [g g], [y_limits(1) y_limits(1)-tick_arm_y], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
    end
    % Manual y-axis line and ticks (first to last tick only)
    line(ax, [0.5 0.5], [y_ticks(1) y_ticks(end)], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
    tick_arm_x = 0.03 * n_groups;
    for ti = 1:numel(y_ticks)
        line(ax, [0.5 0.5-tick_arm_x], [y_ticks(ti) y_ticks(ti)], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
        text(ax, 0.5 - tick_arm_x - 0.05, y_ticks(ti), num2str(y_ticks(ti)), ...
            'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
    ylabel(ax, y_label, 'FontSize', 8, 'Color', 'k');

    % Two-row x-axis labels: treatment on tick row, cell type below
    set(ax, 'XTick', 1:n_groups, 'XTickLabel', []);
    tx_labels = {'ctrl', '{\ittutl-}', 'ctrl', '{\ittutl-}'};
    tx_colors = [0 0 0; 1 0 0; 0.4 0.4 0.4; 0.8 0.2 0.2];
    for g = 1:n_groups
        text(ax, g, y_limits(1) - 0.06*diff(y_limits), tx_labels{g}, ...
            'FontSize', 6, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', 'Interpreter', 'tex', ...
            'Color', tx_colors(g,:));
    end
    % Cell type group labels
    text(ax, 1.5, y_limits(1) - 0.16*diff(y_limits), 'T4', ...
        'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');
    text(ax, 3.5, y_limits(1) - 0.16*diff(y_limits), 'T5', ...
        'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');

    % Wilcoxon brackets: T4 ctrl vs T4 tutl-, T5 ctrl vs T5 tutl-
    % Bracket height: ~4% of data range; y0 near top of data
    bracket_dy = 0.04 * diff(y_limits);
    y0 = y_limits(2) - 0.02 * diff(y_limits);

    if numel(group_data{1}) >= 2 && numel(group_data{2}) >= 2
        p = ranksum(group_data{1}, group_data{2});
        draw_bracket(ax, 1, 2, y0, p);
    end
    if numel(group_data{3}) >= 2 && numel(group_data{4}) >= 2
        p = ranksum(group_data{3}, group_data{4});
        draw_bracket(ax, 3, 4, y0, p);
    end
    % T4 ctrl vs T5 ctrl bracket (wider, higher)
    if numel(group_data{1}) >= 2 && numel(group_data{3}) >= 2
        p = ranksum(group_data{1}, group_data{3});
        draw_bracket(ax, 1, 3, y0 + bracket_dy, p);
    end
end


function draw_ar_panel(ax, results)
% DRAW_AR_PANEL  Aspect ratio box+dot plot for 4 groups.

    group_masks = {
        [results.is_on] & ~[results.is_ttl], ...
        [results.is_on] &  [results.is_ttl], ...
        ~[results.is_on] & ~[results.is_ttl], ...
        ~[results.is_on] &  [results.is_ttl]
    };

    colors = [0 0 0; 1 0 0; 0.4 0.4 0.4; 0.8 0.2 0.2];
    box_colors = min(1, colors * 0.3 + 0.7);

    n = numel(results);
    ar_vals = NaN(1, n);
    for k = 1:n
        d = results(k).max_v_aligned;
        if isnumeric(d) && size(d, 1) == 16 && size(d, 2) == 2
            ar_vals(k) = compute_ar(d);
        end
    end

    n_groups = 4;
    all_vals = [];
    all_grp_idx = [];
    group_data = cell(1, n_groups);

    for g = 1:n_groups
        mask = group_masks{g};
        vals = ar_vals(mask);
        vals = vals(~isnan(vals));
        group_data{g} = vals;
        ng = numel(vals);
        all_vals = [all_vals, vals]; %#ok<AGROW>
        all_grp_idx = [all_grp_idx, repmat(g, 1, ng)]; %#ok<AGROW>
    end

    if isempty(all_vals), return; end

    y_limits = [0 4.4];
    y_ticks = 0:1:4;

    hold(ax, 'on');

    % AR=1 reference line
    yline(ax, 1, '-', 'Color', 'k', 'LineWidth', 0.5);

    % Color-matched box plots per group
    for g = 1:n_groups
        idx = all_grp_idx == g;
        if sum(idx) < 2, continue; end
        bc = boxchart(ax, all_grp_idx(idx), all_vals(idx), ...
            'BoxFaceColor', box_colors(g,:), 'BoxEdgeColor', colors(g,:), ...
            'WhiskerLineColor', colors(g,:), 'MarkerStyle', 'none', ...
            'BoxWidth', 0.5, 'LineWidth', 1.0);
    end

    % Jittered dots
    for g = 1:n_groups
        idx = all_grp_idx == g;
        vals = all_vals(idx);
        x = g + 0.25 * (rand(size(vals)) - 0.5);
        scatter(ax, x, vals, 20, colors(g, :), 'filled', 'MarkerFaceAlpha', 0.5);
    end

    ylim(ax, y_limits);
    ax.Clipping = 'off';  % prevent brackets from being clipped
    set(ax, 'YTick', y_ticks, 'TickDir', 'out', 'FontSize', 7);
    xlim(ax, [0.5, n_groups + 0.5]);
    box(ax, 'off');
    % Floating axes: hide defaults, draw manual lines between first/last ticks
    ax.XColor = 'none';
    ax.YColor = 'none';
    line(ax, [1 n_groups], [y_limits(1) y_limits(1)], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
    tick_arm_y = 0.015 * diff(y_limits);
    for g = 1:n_groups
        line(ax, [g g], [y_limits(1) y_limits(1)-tick_arm_y], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
    end
    line(ax, [0.5 0.5], [y_ticks(1) y_ticks(end)], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
    tick_arm_x = 0.03 * n_groups;
    for ti = 1:numel(y_ticks)
        line(ax, [0.5 0.5-tick_arm_x], [y_ticks(ti) y_ticks(ti)], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
        text(ax, 0.5 - tick_arm_x - 0.05, y_ticks(ti), num2str(y_ticks(ti)), ...
            'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
    ylabel(ax, {'aspect ratio', 'PD/mean(OD_1,OD_2)'}, 'FontSize', 8, 'Color', 'k');

    % Two-row x-axis labels (color-coded)
    set(ax, 'XTick', 1:n_groups, 'XTickLabel', []);
    tx_labels = {'ctrl', '{\ittutl-}', 'ctrl', '{\ittutl-}'};
    tx_colors = [0 0 0; 1 0 0; 0.4 0.4 0.4; 0.8 0.2 0.2];
    for g = 1:n_groups
        text(ax, g, y_limits(1) - 0.06*diff(y_limits), tx_labels{g}, ...
            'FontSize', 6, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', 'Interpreter', 'tex', ...
            'Color', tx_colors(g,:));
    end
    text(ax, 1.5, y_limits(1) - 0.16*diff(y_limits), 'T4', ...
        'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');
    text(ax, 3.5, y_limits(1) - 0.16*diff(y_limits), 'T5', ...
        'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');

    % Wilcoxon brackets
    % Bracket height: ~4% of data range; y0 near top of data
    bracket_dy = 0.04 * diff(y_limits);
    y0 = y_limits(2) - 0.02 * diff(y_limits);

    if numel(group_data{1}) >= 2 && numel(group_data{2}) >= 2
        p = ranksum(group_data{1}, group_data{2});
        draw_bracket(ax, 1, 2, y0, p);
    end
    if numel(group_data{3}) >= 2 && numel(group_data{4}) >= 2
        p = ranksum(group_data{3}, group_data{4});
        draw_bracket(ax, 3, 4, y0, p);
    end
    % T4 ctrl vs T5 ctrl bracket (wider, higher)
    if numel(group_data{1}) >= 2 && numel(group_data{3}) >= 2
        p = ranksum(group_data{1}, group_data{3});
        draw_bracket(ax, 1, 3, y0 + bracket_dy, p);
    end
    ylim(ax, y_limits);
end


function [axPolar, axFill] = plot_polar_with_patch(ax_position, ...
    theta, center_ctrl, spread_ctrl, center_ttl, spread_ttl, opts)
% PLOT_POLAR_WITH_PATCH  Central polar plot with filled patch SEM shading.
%   Dual-axis: polaraxes for grid, Cartesian overlay for patch/line.

    if nargin < 7, opts = struct(); end
    if ~isfield(opts, 'ctrl_line_color'), opts.ctrl_line_color = [0 0 0]; end
    if ~isfield(opts, 'ctrl_fill_color'), opts.ctrl_fill_color = [0.80 0.80 0.80]; end
    if ~isfield(opts, 'ttl_line_color'),  opts.ttl_line_color  = [1 0 0]; end
    if ~isfield(opts, 'ttl_fill_color'),  opts.ttl_fill_color  = [1 0.70 0.70]; end
    if ~isfield(opts, 'alpha'),           opts.alpha           = 0.40; end
    if ~isfield(opts, 'ctrl_label'),      opts.ctrl_label      = 'control'; end
    if ~isfield(opts, 'ttl_label'),       opts.ttl_label       = '{\ittutl-}'; end
    if ~isfield(opts, 'dir_stats'),      opts.dir_stats       = []; end
    if ~isfield(opts, 'rpad_factor'),    opts.rpad_factor     = 1.10; end

    axPolar = polaraxes('Position', ax_position);
    hold(axPolar, 'on');

    all_upper = [];
    if ~isempty(center_ctrl)
        all_upper = [all_upper; center_ctrl(:) + spread_ctrl(:)];
    end
    if ~isempty(center_ttl)
        all_upper = [all_upper; center_ttl(:) + spread_ttl(:)];
    end
    % rlim at 30 mV — the boundary circle IS the 30 mV ring
    % RTick at 15 only; the rlim boundary draws the outer (30 mV) circle
    rlim(axPolar, [0, 30]);
    axPolar.RTick = [15];

    axFill = axes('Position', axPolar.Position, 'Color', 'none', ...
                  'XColor', 'none', 'YColor', 'none', 'HitTest', 'off');
    axis(axFill, 'equal');
    hold(axFill, 'on');

    rL = rlim(axPolar);
    rpad = rL(2) * opts.rpad_factor;  % overlay padding (1.10 main / 1.15 supp)
    set(axFill, 'XLim', [-rpad rpad], 'YLim', [-rpad rpad]);

    % Explicit black circle at 30 mV
    th_circ = linspace(0, 2*pi, 200);
    plot(axFill, 30*cos(th_circ), 30*sin(th_circ), 'k-', 'LineWidth', 0.4);


    hCtrl = gobjects(1, 1);
    if ~isempty(center_ctrl) && any(~isnan(center_ctrl))
        hCtrl = draw_polar_patch(axFill, theta, center_ctrl, spread_ctrl, ...
            opts.ctrl_line_color, opts.ctrl_fill_color, opts.alpha, 0.75);
    end

    hTtl = gobjects(1, 1);
    if ~isempty(center_ttl) && any(~isnan(center_ttl))
        hTtl = draw_polar_patch(axFill, theta, center_ttl, spread_ttl, ...
            opts.ttl_line_color, opts.ttl_fill_color, opts.alpha, 0.75);
    end

    set(axFill, 'XLim', [-rpad rpad], 'YLim', [-rpad rpad]);
    set(axFill, 'Position', axPolar.Position);

    axPolar.Color = 'none';
    uistack(axFill, 'top');

    % Legend is now drawn in draw_ring_panel as figure-level annotations

    axPolar.ThetaZeroLocation = 'right';
    axPolar.ThetaDir = 'counterclockwise';
    axPolar.ThetaTick = 0:22.5:337.5;
    axPolar.ThetaTickLabel = {};
    axPolar.RTickLabel = {};  % remove 0/10/20 radial labels
    axPolar.FontSize = 5;
    axPolar.LineWidth = 0.5;

    % "30 mV" label at 3 o'clock (0 deg = right side)
    rL = rlim(axPolar);
    text(axFill, rL(2)*1.05, 0, '30 mV', 'FontSize', 5, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

    % --- Per-direction significance asterisks (raw p, no FDR) ---
    if ~isempty(opts.dir_stats) && numel(opts.dir_stats) == 16
        rL = rlim(axPolar);
        r_ast = rL(2) * 1.20;  % clearly outside the outer ring
        for di = 1:16
            p_raw = opts.dir_stats(di).p_raw;
            if isnan(p_raw) || p_raw >= 0.05, continue; end
            if p_raw < 0.001
                ast_str = '***';
            elseif p_raw < 0.01
                ast_str = '**';
            else
                ast_str = '*';
            end
            th_ast = opts.dir_stats(di).angle_rad;
            [x_ast, y_ast] = pol2cart(th_ast, r_ast);
            text(axFill, x_ast, y_ast, ast_str, 'FontSize', 6, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', 'Color', 'k');
        end
    end
end


function hLine = draw_polar_patch(axFill, theta, centerVals, bandVals, ...
    lineColor, fillColor, alphaVal, lw)
% DRAW_POLAR_PATCH  Draw center line + spread patch on Cartesian overlay.

    th = theta(:);
    m  = centerVals(:);
    b  = bandVals(:);

    upper = max(m + b, 0);
    lower = max(m - b, 0);

    th_closed    = [th; th(1) + 2*pi];
    upper_closed = [upper; upper(1)];
    lower_closed = [lower; lower(1)];
    m_closed     = [m; m(1)];

    th_poly = [th_closed; flipud(th_closed)];
    r_poly  = [upper_closed; flipud(lower_closed)];
    [x_poly, y_poly] = pol2cart(th_poly, r_poly);

    patch('XData', x_poly, 'YData', y_poly, ...
          'FaceColor', fillColor, 'FaceAlpha', alphaVal, ...
          'EdgeColor', 'none', 'Parent', axFill, 'HitTest', 'off');

    [x_line, y_line] = pol2cart(th_closed, m_closed);
    hLine = plot(axFill, x_line, y_line, '-', 'Color', lineColor, 'LineWidth', lw);
end


function add_scale_bar_on_axes(ax, t_ms, y_lim, font_size, bar_t)
% ADD_SCALE_BAR_ON_AXES  Draw time + voltage L-shaped scale bars on left side.
    if isempty(ax) || ~isvalid(ax), return; end
    if nargin < 4, font_size = 5; end
    if nargin < 5 || isempty(bar_t), bar_t = 1000; end  % default 1 s

    bar_v = 15;    % 15 mV
    v_range = diff(y_lim);

    % Time label
    if bar_t >= 1000
        t_label = sprintf('%g s', bar_t / 1000);
    else
        t_label = sprintf('%g ms', bar_t);
    end

    % Position at bottom-left of subplot
    xl = xlim(ax);
    x_left   = xl(1) + 0.02 * diff(xl);
    y_bottom = y_lim(1) + 0.05 * v_range;
    x_right  = x_left + bar_t;
    y_top    = y_bottom + bar_v;

    % L-shaped bars
    plot(ax, [x_left, x_right], [y_bottom, y_bottom], 'k-', 'LineWidth', 1, ...
        'Clipping', 'off');
    plot(ax, [x_left, x_left], [y_bottom, y_top], 'k-', 'LineWidth', 1, ...
        'Clipping', 'off');

    text(ax, x_left, y_bottom - 0.08*v_range, t_label, ...
        'HorizontalAlignment', 'left', 'FontSize', font_size, 'Clipping', 'off');

    % "15 mV" — outside the vertical bar, starting at base
    text(ax, x_left - 0.03*diff(xl), y_bottom, '15 mV', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'FontSize', font_size, 'Rotation', 90, 'Clipping', 'off');
end


function draw_bracket(ax, x1, x2, y, p)
% DRAW_BRACKET  Draw a bracket between two x positions with p-value.
    yl = ylim(ax);
    arm = 0.02 * diff(yl);  % 2% of axis range
    line(ax, [x1 x1 x2 x2], [y - arm, y, y, y - arm], ...
        'Color', 'k', 'LineWidth', 0.25);
    if p < 0.001
        p_str = '***';
    elseif p < 0.01
        p_str = '**';
    elseif p < 0.05
        p_str = '*';
    else
        p_str = 'n.s.';
    end
    text(ax, mean([x1 x2]), y, p_str, 'FontSize', 6, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end


function ar = compute_ar(d_aligned)
% COMPUTE_AR  Aspect ratio from 16x2 PD-aligned tuning data.
    angles = d_aligned(:, 1);
    resps  = d_aligned(:, 2);

    [~, pd_idx]     = min(abs(angles - pi/2));
    [~, ortho1_idx] = min(abs(angles - 0));
    [~, ortho2_idx] = min(abs(angles - pi));

    pd_resp    = resps(pd_idx);
    ortho_mean = mean([resps(ortho1_idx), resps(ortho2_idx)]);

    if ortho_mean > 0
        ar = pd_resp / ortho_mean;
    else
        ar = NaN;
    end
end


function [traces_16, pd_shift, lut_dirs_ordered] = extract_and_align_traces( ...
    bar_data, plot_order, lut_directions, pd_direction, row_offset)
% EXTRACT_AND_ALIGN_TRACES  Extract 16 mean traces and PD-align by circular shift.
%   row_offset: 0 for 28 dps (rows 1-16), 16 for 56 dps, 32 for 168 dps.

    if nargin < 5, row_offset = 0; end

    n_reps = size(bar_data, 2) - 1;

    traces_subplot_order = cell(16, 1);
    for si = 1:16
        data_row = plot_order(si) + row_offset;
        traces_subplot_order{si} = bar_data{data_row, n_reps + 1};
    end

    lut_dirs_subplot = lut_directions(plot_order);
    [sorted_angles, sort_idx] = sort(lut_dirs_subplot(:));
    traces_sorted = traces_subplot_order(sort_idx);

    angle_diffs = abs(mod(sorted_angles - pd_direction + 180, 360) - 180);
    [~, pd_sorted_idx] = min(angle_diffs);

    pd_shift = 5 - pd_sorted_idx;
    traces_16 = circshift(traces_sorted, pd_shift);
    lut_dirs_ordered = circshift(sorted_angles, pd_shift);
end


function shifted = shift_and_window(tr, peak_time, ref_sample, display_len)
% SHIFT_AND_WINDOW  Shift trace so peak_time -> ref_sample, window to display_len.
    shifted = NaN(display_len, 1);
    offset = round(peak_time) - ref_sample;
    for out_idx = 1:display_len
        src_idx = out_idx + offset;
        if src_idx >= 1 && src_idx <= numel(tr)
            shifted(out_idx) = tr(src_idx);
        end
    end
end


function shift_time = borrow_nearest_shift(time_to_max_row, valid_row, target_di)
% BORROW_NEAREST_SHIFT  Find nearest valid direction's alignment time (circular).
    shift_time = NaN;
    for dist = 1:8
        cw  = mod(target_di - 1 + dist, 16) + 1;
        ccw = mod(target_di - 1 - dist, 16) + 1;
        if valid_row(cw)
            shift_time = time_to_max_row(cw);
            return;
        end
        if valid_row(ccw)
            shift_time = time_to_max_row(ccw);
            return;
        end
    end
end


function exp_list = discover_early_experiments(data_root)
% DISCOVER_EARLY_EXPERIMENTS  Walk {control,ttl}/{ON,OFF}/ tree.
    exp_list = struct('folder', {}, 'date_str', {}, ...
        'treatment', {}, 'cell_type', {});

    for treatment = ["control", "ttl"]
        for cell_type = ["ON", "OFF"]
            subdir = fullfile(data_root, treatment, cell_type);
            if ~isfolder(subdir), continue; end

            dd = dir(subdir);
            dd = dd([dd.isdir]);
            dd = dd(~startsWith({dd.name}, '.'));

            for k = 1:numel(dd)
                exp_folder = fullfile(subdir, dd(k).name);
                if ~isfile(fullfile(exp_folder, 'currentExp.mat'))
                    continue;
                end
                ei.folder    = exp_folder;
                ei.date_str  = dd(k).name;
                ei.treatment = char(treatment);
                ei.cell_type = char(cell_type);
                exp_list(end + 1) = ei; %#ok<AGROW>
            end
        end
    end
end


function entry = find_batch_entry(batch_results, folder_name_or_path)
% FIND_BATCH_ENTRY  Look up a batch results entry by folder name or path.
    entry = [];
    [~, query_base] = fileparts(folder_name_or_path);
    if isempty(query_base)
        query_base = folder_name_or_path;
    end

    for i = 1:numel(batch_results)
        batch_folder = batch_results(i).folder;
        [~, batch_base] = fileparts(batch_folder);
        if strcmp(batch_base, query_base) || ...
           endsWith(batch_folder, folder_name_or_path) || ...
           strcmp(batch_folder, folder_name_or_path)
            entry = batch_results(i);
            return;
        end
    end
end


function [stats, tbl] = compute_direction_stats(all_cells, ctrl_idx, ttl_idx, pd_aligned_angles)
% COMPUTE_DIRECTION_STATS  Wilcoxon rank-sum at each of 16 PD-aligned directions.
%   Returns stats(16) struct array and tbl (16-row table) for printing.
%   FDR correction: Benjamini-Hochberg at q=0.05.

    n_ctrl = numel(ctrl_idx);
    n_ttl  = numel(ttl_idx);
    n_dirs = 16;

    stats = struct('angle_deg', cell(n_dirs,1), 'angle_rad', cell(n_dirs,1), ...
        'n_ctrl', cell(n_dirs,1), 'n_ttl', cell(n_dirs,1), ...
        'median_ctrl', cell(n_dirs,1), 'median_ttl', cell(n_dirs,1), ...
        'mean_ctrl', cell(n_dirs,1), 'mean_ttl', cell(n_dirs,1), ...
        'p_raw', cell(n_dirs,1), 'p_fdr', cell(n_dirs,1), ...
        'sig_raw', cell(n_dirs,1), 'sig_fdr', cell(n_dirs,1));

    p_raw_all = NaN(n_dirs, 1);

    for di = 1:n_dirs
        % Gather per-cell peak amplitudes at this direction
        ctrl_vals = NaN(n_ctrl, 1);
        for k = 1:n_ctrl
            pa = all_cells(ctrl_idx(k)).peak_amps;
            if ~isempty(pa) && size(pa, 1) >= di
                ctrl_vals(k) = pa(di, 2);
            end
        end
        ttl_vals = NaN(n_ttl, 1);
        for k = 1:n_ttl
            pa = all_cells(ttl_idx(k)).peak_amps;
            if ~isempty(pa) && size(pa, 1) >= di
                ttl_vals(k) = pa(di, 2);
            end
        end

        ctrl_vals = ctrl_vals(~isnan(ctrl_vals));
        ttl_vals  = ttl_vals(~isnan(ttl_vals));

        stats(di).angle_deg   = pd_aligned_angles(di);
        stats(di).angle_rad   = deg2rad(pd_aligned_angles(di));
        stats(di).n_ctrl      = numel(ctrl_vals);
        stats(di).n_ttl       = numel(ttl_vals);
        stats(di).median_ctrl = median(ctrl_vals);
        stats(di).median_ttl  = median(ttl_vals);
        stats(di).mean_ctrl   = mean(ctrl_vals);
        stats(di).mean_ttl    = mean(ttl_vals);

        if numel(ctrl_vals) >= 2 && numel(ttl_vals) >= 2
            p_raw_all(di) = ranksum(ctrl_vals, ttl_vals);
        else
            p_raw_all(di) = NaN;
        end
        stats(di).p_raw = p_raw_all(di);
        stats(di).sig_raw = ~isnan(p_raw_all(di)) && p_raw_all(di) < 0.05;
    end

    % Benjamini-Hochberg FDR correction
    valid = ~isnan(p_raw_all);
    p_valid = p_raw_all(valid);
    n_valid = numel(p_valid);
    valid_idx = find(valid);

    [p_sorted, sort_order] = sort(p_valid);
    p_fdr_sorted = p_sorted;
    for i = 1:n_valid
        p_fdr_sorted(i) = p_sorted(i) * n_valid / i;
    end
    % Enforce monotonicity (from last to first)
    for i = n_valid-1:-1:1
        p_fdr_sorted(i) = min(p_fdr_sorted(i), p_fdr_sorted(i+1));
    end
    p_fdr_sorted = min(p_fdr_sorted, 1);

    % Map back
    p_fdr_all = NaN(n_dirs, 1);
    p_fdr_unsorted = NaN(n_valid, 1);
    p_fdr_unsorted(sort_order) = p_fdr_sorted;
    for i = 1:n_valid
        p_fdr_all(valid_idx(i)) = p_fdr_unsorted(i);
    end

    for di = 1:n_dirs
        stats(di).p_fdr   = p_fdr_all(di);
        stats(di).sig_fdr  = ~isnan(p_fdr_all(di)) && p_fdr_all(di) < 0.05;
    end

    % Build table for printing
    tbl = struct('angle_deg', {stats.angle_deg}, ...
        'n_ctrl', {stats.n_ctrl}, 'n_ttl', {stats.n_ttl}, ...
        'median_ctrl', {stats.median_ctrl}, 'median_ttl', {stats.median_ttl}, ...
        'mean_ctrl', {stats.mean_ctrl}, 'mean_ttl', {stats.mean_ttl}, ...
        'p_raw', {stats.p_raw}, 'p_fdr', {stats.p_fdr}, ...
        'sig_fdr', {stats.sig_fdr});
end


function print_direction_stats_table(tbl)
% PRINT_DIRECTION_STATS_TABLE  Print per-direction stats to console.
    fprintf('%-8s  %5s  %5s  %8s  %8s  %8s  %8s  %8s  %8s  %4s\n', ...
        'Dir(deg)', 'nCtrl', 'nTTL', 'MdnCtrl', 'MdnTTL', 'MeanCtrl', 'MeanTTL', ...
        'p(raw)', 'p(FDR)', 'Sig');
    fprintf('%s\n', repmat('-', 1, 85));
    n = numel([tbl.angle_deg]);
    angles   = [tbl.angle_deg];
    n_ctrls  = [tbl.n_ctrl];
    n_ttls   = [tbl.n_ttl];
    mdn_c    = [tbl.median_ctrl];
    mdn_t    = [tbl.median_ttl];
    mn_c     = [tbl.mean_ctrl];
    mn_t     = [tbl.mean_ttl];
    p_raws   = [tbl.p_raw];
    p_fdrs   = [tbl.p_fdr];
    sigs     = [tbl.sig_fdr];
    for i = 1:n
        sig_str = '';
        if sigs(i)
            if p_fdrs(i) < 0.001, sig_str = '***';
            elseif p_fdrs(i) < 0.01, sig_str = '**';
            else, sig_str = '*';
            end
        end
        fprintf('%7.1f   %5d  %5d  %8.2f  %8.2f  %8.2f  %8.2f  %8.4f  %8.4f  %4s\n', ...
            angles(i), n_ctrls(i), n_ttls(i), mdn_c(i), mdn_t(i), ...
            mn_c(i), mn_t(i), p_raws(i), p_fdrs(i), sig_str);
    end
    n_sig = sum(sigs);
    fprintf('\n  Significant (FDR<0.05): %d / %d directions\n', n_sig, n);
end


function write_direction_stats_table(fid, tbl)
% WRITE_DIRECTION_STATS_TABLE  Write per-direction stats to file.
    fprintf(fid, '%-8s  %5s  %5s  %8s  %8s  %8s  %8s  %8s  %8s  %4s\n', ...
        'Dir(deg)', 'nCtrl', 'nTTL', 'MdnCtrl', 'MdnTTL', 'MeanCtrl', 'MeanTTL', ...
        'p(raw)', 'p(FDR)', 'Sig');
    fprintf(fid, '%s\n', repmat('-', 1, 85));
    n = numel([tbl.angle_deg]);
    angles   = [tbl.angle_deg];
    n_ctrls  = [tbl.n_ctrl];
    n_ttls   = [tbl.n_ttl];
    mdn_c    = [tbl.median_ctrl];
    mdn_t    = [tbl.median_ttl];
    mn_c     = [tbl.mean_ctrl];
    mn_t     = [tbl.mean_ttl];
    p_raws   = [tbl.p_raw];
    p_fdrs   = [tbl.p_fdr];
    sigs     = [tbl.sig_fdr];
    for i = 1:n
        sig_str = '';
        if sigs(i)
            if p_fdrs(i) < 0.001, sig_str = '***';
            elseif p_fdrs(i) < 0.01, sig_str = '**';
            else, sig_str = '*';
            end
        end
        fprintf(fid, '%7.1f   %5d  %5d  %8.2f  %8.2f  %8.2f  %8.2f  %8.4f  %8.4f  %4s\n', ...
            angles(i), n_ctrls(i), n_ttls(i), mdn_c(i), mdn_t(i), ...
            mn_c(i), mn_t(i), p_raws(i), p_fdrs(i), sig_str);
    end
    n_sig = sum(sigs);
    fprintf(fid, '\n  Significant (FDR<0.05): %d / %d directions\n', n_sig, n);
end


function r_out = harmonize_results(r_in, batch_label)
% HARMONIZE_RESULTS  Extract common fields into a uniform struct array.
    n = numel(r_in);
    r_out = struct([]);

    for k = 1:n
        s.date_str        = r_in(k).date_str;
        s.is_on           = r_in(k).is_on;
        s.is_ttl          = r_in(k).is_ttl;
        s.group           = r_in(k).group;
        s.max_v_aligned   = r_in(k).max_v_aligned;
        s.dsi_vector      = r_in(k).dsi_vector;
        s.dsi_pdnd        = r_in(k).dsi_pdnd;
        s.dir_tuning_fwhm = r_in(k).dir_tuning_fwhm;
        s.batch           = batch_label;

        if isempty(r_out)
            r_out = s;
        else
            r_out(end + 1) = s; %#ok<AGROW>
        end
    end
end


function y = blockavg_downsample(x, n)
% BLOCKAVG_DOWNSAMPLE  Anti-aliased downsampling via block averaging.
%   Local helper: reduces sample rate of vector X by factor N by averaging
%   N consecutive samples per output sample.

    if n == 1, y = x; return; end
    was_row = isrow(x);
    if was_row, x = x(:); end
    len    = size(x, 1);
    n_full = floor(len / n);
    n_out  = ceil(len / n);
    y = zeros(n_out, size(x, 2));
    if n_full > 0
        y(1:n_full, :) = mean(reshape(x(1:n_full*n, :), n, n_full, []), 1);
    end
    if n_out > n_full
        y(end, :) = mean(x(n_full*n+1:end, :), 1);
    end
    if was_row, y = y'; end
end
