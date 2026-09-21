function results = batch_analyze_pre_bar_flash(data_root, opts)
% BATCH_ANALYZE_PRE_BAR_FLASH  Batch-process pre-bar-flash experiments.
%
%   RESULTS = BATCH_ANALYZE_PRE_BAR_FLASH(DATA_ROOT) processes all experiments
%   under the pre-bar-flash directory tree. Extracts directional tuning from
%   28 dps bar sweeps and computes per-cell DS metrics.
%
%   This pipeline is a simplified version of batch_analyze_1DRF for the summer
%   2025 protocol which has bar sweeps at 2 speeds (28 + 56 dps) but NO bar
%   flashes, NO 168 dps speed, and NO RF mapping data.
%
%   DATA_ROOT should point to the pre-bar-flash directory, e.g.:
%     /Users/reiserm/Documents/ttl_1DRF/pre-bar-flash/
%
%   Expected directory structure:
%     DATA_ROOT/{control,ttl}/{ON,OFF}/YYYY_MM_DD_HH_MM/
%
%   Cell classification is done by DIRECTORY PATH (not metadata), because some
%   post-Oct-15 experiments have inconsistent metadata fields.
%
%   RESULTS = BATCH_ANALYZE_PRE_BAR_FLASH(DATA_ROOT, OPTS) uses options:
%     opts.save_results - Save batch_results_pre_bf.mat (default: true)
%     opts.save_dir     - Output directory (default: <data_root>/population_results)
%     opts.lut_path     - Path to bar_lut.mat
%     opts.plot_order   - 1x16 bar ordering (default: standard)
%     opts.baseline_range - Baseline sample range (default: [1000 9000])
%     opts.stim_trim_end  - Trim from end of stim (default: 7000)
%     opts.percentile     - Peak detection percentile (default: 98)
%
%   OUTPUT:
%     results - Structure array (one per cell) with fields:
%       .folder            - Full path to experiment folder
%       .date_str          - Experiment date string
%       .strain            - Strain from metadata
%       .frame             - Frame number from metadata
%       .pd_direction_vs   - continuous vector-sum preferred direction (deg), same as
%                            .pd_direction here (this batch never snapped)
%       .metadata_mismatch - true if Strain / Frame in metadata disagree with
%                            the directory-based classification (warning issued)
%       .is_on             - true for ON/T4 cells
%       .is_ttl            - true for TTL cells
%       .group             - 'on_control', 'on_ttl', 'off_control', 'off_ttl'
%       .max_v_aligned     - 16x2 PD-aligned [angles, responses]
%       .pd_direction      - Preferred direction (degrees)
%       .dsi_vector        - DSI via vector sum (0-1)
%       .dsi_pdnd          - DSI via PD-ND method (-1 to 1)
%       .dir_tuning_fwhm   - Direction tuning FWHM (degrees)
%       .dir_tuning_cv     - Circular variance (0-1)
%       .dir_tuning_kappa  - Von Mises concentration
%       .tuning_aspect_ratio - PD / mean(ortho)
%       .sym_ratio         - Symmetry ratio (0-1)
%       .polarity_corrected - true if dark-bar polarity swap was applied
%
%   See also BATCH_ANALYZE_1DRF, PARSE_BAR_DATA_PRE_BF,
%            COMPUTE_BAR_SWEEP_RESPONSES, FIND_PD_AND_ORDER_IDX

    %% Set defaults
    if nargin < 2, opts = struct(); end
    opts = set_defaults(opts, data_root);

    %% Verify CircStat is on the path (used indirectly via find_PD_and_order_idx)
    if isempty(which('circ_vmpar'))
        error('batch_analyze_pre_bar_flash:CircStatMissing', ...
            ['CircStat toolbox not on MATLAB path. ' ...
             'Add CircStat2012a (or equivalent) to the path before calling. ' ...
             'See https://github.com/circstat/circstat-matlab.']);
    end

    %% Step 1: Load LUT
    S_lut = load(opts.lut_path, 'Tbl');
    Tbl = S_lut.Tbl;

    %% Step 2: Discover experiment folders
    exp_list = discover_experiments(data_root);
    n_exp = numel(exp_list);
    fprintf('Found %d experiments in %s\n\n', n_exp, data_root);

    %% Step 3: Process each experiment
    results = struct([]);

    for exp_idx = 1:n_exp
        ei = exp_list(exp_idx);
        fprintf('[%d/%d] %s  (%s/%s) ...', ...
            exp_idx, n_exp, ei.date_str, ei.treatment, ei.cell_type);

        try
            r = process_single_cell(ei.folder, Tbl, opts, ei);
            r.folder = ei.folder;

            if isempty(results)
                results = r;
            else
                results(end + 1) = r; %#ok<AGROW>
            end

            corr_tag = '';
            if r.polarity_corrected
                corr_tag = ' *polarity-corrected*';
            end
            fprintf(' PD=%.0f° DSI=%.2f AR=%.2f [%s]%s\n', ...
                r.pd_direction, r.dsi_vector, r.tuning_aspect_ratio, ...
                r.group, corr_tag);

        catch ME
            fprintf(' ERROR: %s\n', ME.message);
            continue;
        end
    end

    fprintf('\n=== Processing Complete ===\n');
    fprintf('Total cells: %d\n', numel(results));

    % Per-group counts
    if ~isempty(results)
        groups = {results.group};
        for g = ["on_control", "on_ttl", "off_control", "off_ttl"]
            fprintf('  %s: %d\n', g, sum(strcmp(groups, g)));
        end
    end

    %% Step 4: Save results
    if opts.save_results && ~isempty(results)
        if ~isfolder(opts.save_dir)
            mkdir(opts.save_dir);
        end
        save_path = fullfile(opts.save_dir, 'batch_results_pre_bf.mat');
        save(save_path, 'results');
        fprintf('\nResults saved to: %s\n', save_path);
    end

    %% Step 5: Print summary table
    if ~isempty(results)
        print_summary_table(results);
    end

end


%% ========================= Experiment Discovery ==========================

function exp_list = discover_experiments(data_root)
% DISCOVER_EXPERIMENTS  Walk {control,ttl}/{ON,OFF}/ tree and collect folders.
%   Classification is by directory path, not metadata.

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

                % Must contain currentExp.mat
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


%% ========================= Core Processing ============================

function r = process_single_cell(exp_folder, Tbl, opts, ei)
% PROCESS_SINGLE_CELL  Run the per-cell bar sweep pipeline.

    % Save/restore working directory (load_protocol2_data uses cd)
    orig_dir = pwd;
    cleanup = onCleanup(@() cd(orig_dir));

    % Load data
    [date_str, ~, Log, ~, ~] = load_protocol2_data(exp_folder);

    f_data = Log.ADC.Volts(1, :);
    v_data = Log.ADC.Volts(2, :) * 10;

    ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
        'pattern_order', 'func_order', 'metadata');
    metadata = ce.metadata;

    % --- Classify by DIRECTORY PATH (not metadata) ---
    r.date_str = date_str;
    r.strain   = metadata.Strain;
    r.frame    = metadata.Frame;

    r.is_on  = strcmpi(ei.cell_type, 'ON');
    r.is_ttl = strcmpi(ei.treatment, 'ttl');

    if r.is_on && ~r.is_ttl
        r.group = 'on_control';
    elseif r.is_on && r.is_ttl
        r.group = 'on_ttl';
    elseif ~r.is_on && ~r.is_ttl
        r.group = 'off_control';
    else
        r.group = 'off_ttl';
    end

    % --- Cross-check the directory-based classification against metadata ---
    %   (a) genotype: when Strain explicitly names 'ttl' or 'control' it must
    %       agree with the ttl/ vs control/ folder (summer sessions carry the
    %       uninformative '42F06_T4T5' and are not checked);
    %   (b) polarity: the RF frame from the localisation pattern lies in the
    %       bright half (> 129, the threshold used by batch_analyze_1DRF) iff
    %       the folder is under ON/.
    %   Warn only; the flag is stored so mismatches can be listed from the results.
    strain_l = lower(string(metadata.Strain));
    strain_informative = contains(strain_l, 'ttl') || contains(strain_l, 'control');
    geno_mismatch = strain_informative && (contains(strain_l, 'ttl') ~= r.is_ttl);
    pol_mismatch  = isnumeric(metadata.Frame) && ((metadata.Frame > 129) ~= r.is_on);
    r.metadata_mismatch = geno_mismatch || pol_mismatch;
    if r.metadata_mismatch
        warning('batch_analyze_pre_bar_flash:MetadataMismatch', ...
            '%s filed as %s but metadata says Strain=%s, Frame=%g -- check the folder', ...
            date_str, r.group, string(metadata.Strain), metadata.Frame);
    end

    % --- Parse bar sweep data (pre-bar-flash protocol) ---
    bar_data = parse_bar_data_pre_bf(f_data, v_data);

    % --- Correct dark-bar polarity swap for pre-Oct-15 OFF cells ---
    %   Patterns 9/10 used bar_0 (dark bar) instead of bar_15 (bright bar),
    %   causing func=3/func=4 directions to be swapped for OFF cells only.
    %   See correct_off_polarity_swap.m for full details.
    [bar_data, r.polarity_corrected] = correct_off_polarity_swap( ...
        bar_data, ce.pattern_order, ce.func_order, ei.date_str, ~r.is_on);

    % --- Verify LUT directions ---
    [lut_directions, ~, ~, ~] = ...
        verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, opts.plot_order);

    % --- Compute bar sweep responses ---
    sweep_opts.baseline_range = opts.baseline_range;
    sweep_opts.stim_trim_end  = opts.stim_trim_end;
    sweep_opts.percentile     = opts.percentile;
    max_v = compute_bar_sweep_responses(bar_data, opts.plot_order, sweep_opts);

    % --- PD-align the tuning curve ---
    % Sort by LUT direction (0, 22.5, ..., 337.5) for find_PD_and_order_idx
    lut_dirs_ordered = lut_directions(opts.plot_order);
    [~, sort_idx] = sort(lut_dirs_ordered);
    max_v_sorted = max_v(sort_idx);
    max_v_polar = [max_v_sorted; max_v_sorted(1)];  % 17x1 for circular

    [d_aligned, ~, ~, angle_rad, dir_fwhm, dir_cv, ~, dir_kappa] = ...
        find_PD_and_order_idx(max_v_polar, 0);
    r.max_v_aligned = d_aligned;

    % PD direction in degrees
    r.pd_direction = rad2deg(angle_rad);
    if r.pd_direction < 0
        r.pd_direction = r.pd_direction + 360;
    end
    r.pd_direction_vs = r.pd_direction;   % continuous vector-sum PD (same field name as the late batch)

    % --- Direction selectivity metrics ---
    [r.sym_ratio, r.dsi_vector, r.dsi_pdnd, ~] = ...
        compute_bar_response_metrics(d_aligned);
    r.dir_tuning_fwhm  = dir_fwhm;
    r.dir_tuning_cv    = dir_cv;
    r.dir_tuning_kappa = dir_kappa;

    % --- Tuning aspect ratio ---
    r.tuning_aspect_ratio = compute_tuning_aspect_ratio(d_aligned);

end


%% ========================= Tuning Aspect Ratio ========================

function ar = compute_tuning_aspect_ratio(d_aligned)
% COMPUTE_TUNING_ASPECT_RATIO  PD response / mean of two orthogonal responses.
%   Matches the implementation in batch_analyze_1DRF.m.

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


%% ========================= Default Options ============================

function opts = set_defaults(opts, data_root)
% SET_DEFAULTS  Fill in default parameter values.

    if ~isfield(opts, 'lut_path')
        % LUT is stored alongside batch_analyze_1DRF.m in the protocol2 dir
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
    if ~isfield(opts, 'save_results')
        opts.save_results = true;
    end
    if ~isfield(opts, 'save_dir')
        opts.save_dir = fullfile(data_root, 'population_results');
    end

end


%% ========================= Summary Table ==============================

function print_summary_table(results)
% PRINT_SUMMARY_TABLE  Print a per-cell metrics table to console.

    fprintf('\n');
    fprintf('%-4s  %-20s  %-12s  %6s  %6s  %6s  %6s  %6s  %6s\n', ...
        '#', 'Date', 'Group', 'PD', 'DSI_v', 'DSI_pn', 'FWHM', 'AR', 'SymR');
    fprintf('%s\n', repmat('-', 1, 92));

    for k = 1:numel(results)
        r = results(k);
        fprintf('%-4d  %-20s  %-12s  %6.1f  %6.3f  %6.3f  %6.1f  %6.2f  %6.3f\n', ...
            k, r.date_str, r.group, r.pd_direction, ...
            r.dsi_vector, r.dsi_pdnd, r.dir_tuning_fwhm, ...
            r.tuning_aspect_ratio, r.sym_ratio);
    end

    fprintf('\n');

end
