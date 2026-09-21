function fig_ef = generate_manuscript_fig_ef(axis_mode, opts)
% GENERATE_MANUSCRIPT_FIG_EF  Flash-row sub-figure (M6-aligned PD or ortho axis).
%
%   FIG_EF = GENERATE_MANUSCRIPT_FIG_EF(AXIS_MODE) builds an 18 x 7 cm
%   figure with 1x11 bar-flash trace tiles (T4 row + T5 row) and
%   right-hand-column amplitude summaries.  AXIS_MODE selects which
%   M6-aligned flash field of BATCH_RESULTS.MAT to plot:
%
%       'pd'    - pd_flash_m6_aligned, PD-on-left flip applied
%       'ortho' - ortho_flash_m6_aligned, no flip
%
%   FIG_EF = GENERATE_MANUSCRIPT_FIG_EF(AXIS_MODE, OPTS) accepts an
%   options struct with fields:
%       .data_root      - data_root (default '/Users/reiserm/Documents/ttl_1DRF')
%       .skip_export    - true to suppress PDF/PNG export (default false)
%       .stamp_path     - if non-empty, save pre-plot variables to this
%                         .mat file and return without plotting.  Used
%                         by the local validation harness.
%       .use_raw_traces - true for absolute voltage (default false)
%       .show_fwhm      - draw FWHM bars on the depolarization panels
%                         (default false = manuscript)
%       .pool_test      - 'percell' (default = manuscript): each cell's mean
%                         over the three positions, rank-sum on cells;
%                         'pooled': three-position pooled rank-sum on all
%                         cell x position values (treats within-cell
%                         neighbours as independent; kept for comparison)
%       .fwhm_method    - 'interp' (default): linear interpolation of the
%                         half-max crossings, undefined (NaN) if a flank is
%                         truncated; 'gauss': FWHM of a Gaussian fitted to the
%                         amplitude-by-position profile (handles truncation)
%
%   Variable naming convention inside this function: per-cell data
%   matrices are c_dep / t_dep / c_hyp / t_hyp (T4 / ON cells) and
%   off_c_* / off_t_* (T5 / OFF cells), regardless of axis. This matches
%   the stamp gate naming so the validation harness is axis-agnostic.
%
%   See also GENERATE_MANUSCRIPT_FIG_DS, GENERATE_MANUSCRIPT_FIG.

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'data_root'),       opts.data_root      = '/Users/reiserm/Documents/ttl_1DRF'; end
if ~isfield(opts, 'skip_export'),     opts.skip_export    = false; end
if ~isfield(opts, 'stamp_path'),      opts.stamp_path     = ''; end
if ~isfield(opts, 'use_raw_traces'),  opts.use_raw_traces = false; end
if ~isfield(opts, 'show_fwhm'),       opts.show_fwhm      = false; end      % manuscript setting
if ~isfield(opts, 'pool_test'),       opts.pool_test      = 'percell'; end  % manuscript setting
if ~isfield(opts, 'fwhm_method'),     opts.fwhm_method    = 'interp'; end

%% ===================== Axis dispatch =====================================
switch lower(axis_mode)
    case 'pd'
        flash_field    = 'pd_flash_m6_aligned';
        baseline_field = 'pd_flash_baselines';
        centroid_field = 'centroid_m6_rounded';
        do_flip        = true;
        cartoon_orient = 'vertical';
        panel_letter_traces = 'C';
        panel_letter_amp    = 'D';
        section_label  = 'PD axis';
        out_tag        = 'fig_ds_panels_EF';
        diag_label     = 'PD';
    case 'ortho'
        flash_field    = 'ortho_flash_m6_aligned';
        baseline_field = 'ortho_flash_baselines';
        centroid_field = 'ortho_centroid_m6_rounded';
        do_flip        = false;
        cartoon_orient = 'horizontal';
        panel_letter_traces = 'A';
        panel_letter_amp    = 'B';
        section_label  = 'Orthogonal axis';
        out_tag        = 'fig_supp_panels_EF_ortho';
        diag_label     = 'Ort';
    otherwise
        error('generate_manuscript_fig_ef:BadAxisMode', ...
            'axis_mode must be ''pd'' or ''ortho'', got %s', axis_mode);
end

USE_RAW_TRACES = opts.use_raw_traces;
DS_FACTOR      = 10;
data_root      = opts.data_root;
out_dir        = fullfile(data_root, 'manuscript_figures');
if ~isfolder(out_dir) && ~opts.skip_export, mkdir(out_dir); end

%% ===================== Load batch results ===============================
res_file = fullfile(data_root, 'population_results', 'batch_results.mat');
fprintf('Loading: %s\n', res_file);
S = load(res_file, 'results');
results = S.results;
fprintf('Loaded %d cells.\n', numel(results));
% Refuse batch files built before the bar-flash pitch correction (they place the
% three 1-px sessions on the wrong grid). Rebuild with scripts/build_batch_results.m.
if ~isfield(results, 'flash_pitch_px') || any(cellfun(@isempty, {results.flash_pitch_px}))
    error('generate_manuscript_fig_ef:OutdatedBatchFile', ...
        ['%s lacks results(k).flash_pitch_px: it was built before the bar-flash pitch ' ...
         'correction. Rebuild it with scripts/build_batch_results.m (batch_analyze_1DRF).'], res_file);
end

%% ===================== Constants =========================================
FONT_NAME  = 'Helvetica';
FONT_AX    = 6;
FONT_LABEL = 7;
FONT_TITLE = 8;
FONT_STAT  = 5.5;

% T4 (ON) colours
COL_T4 = struct('ctrl_line', [0 0 0],           'ttl_line', [1 0 0], ...
                'ctrl_fill', [0.80 0.80 0.80],  'ttl_fill', [1 0.70 0.70], ...
                'alpha', 0.35, 'stim_line', [0.2 0.7 0.2]);

% T5 (OFF) colours
COL_T5 = struct('ctrl_line', [0.4 0.4 0.4],     'ttl_line', [0.8 0.2 0.2], ...
                'ctrl_fill', [0.70 0.70 0.70],  'ttl_fill', [0.90 0.60 0.60], ...
                'alpha', 0.35, 'stim_line', [0.2 0.7 0.2]);

STAT_COLOR = [0 0 0];
ALPHA_AMP  = 0.20;

LINE_W_TRACE = 0.25;
LINE_W_AMP   = 0.25;
MARKER_SZ    = 2;

% Timing / amplitude constants (10 kHz)
STIM_ONSET    = 5001;
STIM_OFFSET   = 5801;
DEP_WINDOW    = [STIM_ONSET, STIM_ONSET + 2000 - 1];
HYP_WINDOW    = [STIM_ONSET + 1000, Inf];
DEP_PCTILE    = 99.9;
HYP_PCTILE    = 0.1;
REJECT_THRESH = 0.5;

TRACE_START = 3701;

% Downsampled timing constants
STIM_ONSET_DS  = ceil(STIM_ONSET  / DS_FACTOR);
STIM_OFFSET_DS = ceil(STIM_OFFSET / DS_FACTOR);
TRACE_START_DS = ceil(TRACE_START / DS_FACTOR);

if USE_RAW_TRACES
    Y_LIM        = [-75 -25];
    YLIM_HYP     = [-75 -50];
    YLIM_AMP_DEP = [-75 -25];
else
    Y_LIM        = [-5 25];
    YLIM_HYP     = [-5  0];
    YLIM_AMP_DEP = [0 25];
end

STAT_METHOD = 'mean_sem';
N_POS       = 11;
TTL_TEX     = '{\ittutl-}';

%% ===================== Group masks =======================================
on_ctrl  = [results.is_on] & ~[results.is_ttl];
on_ttl   = [results.is_on] &  [results.is_ttl];
off_ctrl = ~[results.is_on] & ~[results.is_ttl];
off_ttl  = ~[results.is_on] &  [results.is_ttl];

fprintf('T4 (ON):  ctrl=%d  tutl-=%d\n', sum(on_ctrl), sum(on_ttl));
fprintf('T5 (OFF): ctrl=%d  tutl-=%d\n', sum(off_ctrl), sum(off_ttl));

%% ===================== Reconstruct raw traces if requested ==============
if USE_RAW_TRACES
    fprintf('Reconstructing absolute-voltage traces from baselines...\n');
    for k = 1:numel(results)
        r = results(k);
        bl_aligned = reindex_to_peak(r.(baseline_field), r.(centroid_field), 2 / r.flash_pitch_px);
        results(k).(flash_field) = r.(flash_field) + bl_aligned;
    end
    Y_LABEL = 'mV';
else
    Y_LABEL = '\DeltamV';
end

%% ===================== Trace statistics ==================================
on_ctrl_tr  = {results(on_ctrl).(flash_field)};
on_ttl_tr   = {results(on_ttl).(flash_field)};
off_ctrl_tr = {results(off_ctrl).(flash_field)};
off_ttl_tr  = {results(off_ttl).(flash_field)};

[c_st,    t_st]    = compute_all_stats(on_ctrl_tr,  on_ttl_tr,  N_POS, STAT_METHOD, DS_FACTOR);
[off_c_st, off_t_st] = compute_all_stats(off_ctrl_tr, off_ttl_tr, N_POS, STAT_METHOD, DS_FACTOR);
if do_flip
    % Flip column order: PD on left, ND on right
    c_st     = c_st(11:-1:1);     t_st     = t_st(11:-1:1);
    off_c_st = off_c_st(11:-1:1); off_t_st = off_t_st(11:-1:1);
end

%% ===================== Amplitude matrices ================================
positions = -5:5;  % PD mode: PD at -5 (left), ND at +5 (right). Ortho mode: as-is.

[c_dep,    c_hyp]    = extract_robust_amps(results(on_ctrl),  flash_field, DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[t_dep,    t_hyp]    = extract_robust_amps(results(on_ttl),   flash_field, DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[off_c_dep, off_c_hyp] = extract_robust_amps(results(off_ctrl), flash_field, DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[off_t_dep, off_t_hyp] = extract_robust_amps(results(off_ttl),  flash_field, DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
if do_flip
    c_dep     = c_dep(:,11:-1:1);     c_hyp     = c_hyp(:,11:-1:1);
    t_dep     = t_dep(:,11:-1:1);     t_hyp     = t_hyp(:,11:-1:1);
    off_c_dep = off_c_dep(:,11:-1:1); off_c_hyp = off_c_hyp(:,11:-1:1);
    off_t_dep = off_t_dep(:,11:-1:1); off_t_hyp = off_t_hyp(:,11:-1:1);
end

%% ===================== Pooled rank-sum ===================================
[dep_pp,     pool_centers] = compute_pooled_ranksum(c_dep,     t_dep,     opts.pool_test);
[hyp_pp,     ~]            = compute_pooled_ranksum(c_hyp,     t_hyp,     opts.pool_test);
[off_dep_pp, ~]            = compute_pooled_ranksum(off_c_dep, off_t_dep, opts.pool_test);
[off_hyp_pp, ~]            = compute_pooled_ranksum(off_c_hyp, off_t_hyp, opts.pool_test);
fprintf('\nPooled test (%s), p per pool centre %s:\n', opts.pool_test, mat2str(pool_centers));
fprintf('  %s T4 dep:  %s\n', diag_label, mat2str(round(dep_pp, 3)));
fprintf('  %s T4 hyp:  %s\n', diag_label, mat2str(round(hyp_pp, 3)));
fprintf('  %s T5 dep:  %s\n', diag_label, mat2str(round(off_dep_pp, 3)));
fprintf('  %s T5 hyp:  %s\n', diag_label, mat2str(round(off_hyp_pp, 3)));
fprintf('Cells per grid position (ctrl/tutl), positions %s:\n', mat2str(positions));
fprintf('  T4: %s / %s\n', mat2str(sum(~isnan(c_dep), 1)),     mat2str(sum(~isnan(t_dep), 1)));
fprintf('  T5: %s / %s\n', mat2str(sum(~isnan(off_c_dep), 1)), mat2str(sum(~isnan(off_t_dep), 1)));

%% ===================== STAMP GATE (validation only) =====================
if ~isempty(opts.stamp_path)
    stamp = struct( ...
        'c_dep', c_dep, 't_dep', t_dep, 'c_hyp', c_hyp, 't_hyp', t_hyp, ...
        'off_c_dep', off_c_dep, 'off_t_dep', off_t_dep, ...
        'off_c_hyp', off_c_hyp, 'off_t_hyp', off_t_hyp, ...
        'c_st', c_st, 't_st', t_st, 'off_c_st', off_c_st, 'off_t_st', off_t_st, ...
        'pool_centers', pool_centers, ...
        'dep_pp', dep_pp, 'hyp_pp', hyp_pp, ...
        'off_dep_pp', off_dep_pp, 'off_hyp_pp', off_hyp_pp);
    save(opts.stamp_path, '-struct', 'stamp');
    fprintf('STAMP wrote %s\n', opts.stamp_path);
    fig_ef = [];
    return
end

%% ===================== Figure ============================================
FIG_W = 18;  FIG_H = 7;
fig_ef = figure('Units', 'centimeters', 'Position', [2 2 FIG_W FIG_H], ...
    'PaperUnits', 'centimeters', 'PaperSize', [FIG_W FIG_H], ...
    'PaperPosition', [0 0 FIG_W FIG_H], 'Color', 'w');
set(fig_ef, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

%% ===================== Layout geometry ===================================
TRACE_L  = 0.05;
TRACE_R  = 0.62;
TILE_GAP = 0.001;
TILE_W   = (TRACE_R - TRACE_L - (N_POS-1)*TILE_GAP) / N_POS;

AMP_T4_L = 0.70;
AMP_T5_L = 0.85;
AMP_W    = 0.11;

ROW_H   = 0.30;
ROW_GAP = 0.04;
E_T4_y  = 0.58;
E_T5_y  = E_T4_y - ROW_H - ROW_GAP;

PD_top   = E_T4_y + ROW_H;  PD_bot = E_T5_y;
PD_h     = PD_top - PD_bot;
PD_DEP_H = PD_h * 0.64;     PD_HYP_H = PD_h * 0.26;
PD_DEP_Y = PD_top - PD_DEP_H;
PD_HYP_Y = PD_bot;

%% ===================== Position labels ===================================
pos_labels = arrayfun(@(x) sprintf('%+d', x), 5:-1:-5, 'UniformOutput', false);
pos_labels{6} = '0';

%% ===================== TRACE panels (left) ===============================
draw_trace_row(fig_ef, TRACE_L, E_T4_y, TILE_W, ROW_H, TILE_GAP, ...
    c_st, t_st, N_POS, Y_LIM, COL_T4, LINE_W_TRACE, ...
    pos_labels, true, FONT_AX, FONT_LABEL, STIM_ONSET_DS, STIM_OFFSET_DS, TRACE_START_DS, 'E_T4', Y_LABEL);
add_row_label_colored(fig_ef, 0.01, E_T4_y + ROW_H/2, 'T4', COL_T4.ttl_line, FONT_TITLE);
annotation(fig_ef, 'textbox', [TRACE_L, E_T4_y + ROW_H - 0.04, 0.06, 0.03], ...
    'String', 'control', 'FontSize', 6, 'FontWeight', 'bold', ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', 'Color', COL_T4.ctrl_line);
annotation(fig_ef, 'textbox', [TRACE_L + 0.06, E_T4_y + ROW_H - 0.04, 0.06, 0.03], ...
    'String', '{\ittutl-}', 'Interpreter', 'tex', 'FontSize', 6, 'FontWeight', 'bold', ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', 'Color', COL_T4.ttl_line);

draw_trace_row(fig_ef, TRACE_L, E_T5_y, TILE_W, ROW_H, TILE_GAP, ...
    off_c_st, off_t_st, N_POS, Y_LIM, COL_T5, LINE_W_TRACE, ...
    {}, false, FONT_AX, FONT_LABEL, STIM_ONSET_DS, STIM_OFFSET_DS, TRACE_START_DS, 'E_T5', Y_LABEL);
add_row_label_colored(fig_ef, 0.01, E_T5_y + ROW_H/2, 'T5', COL_T5.ttl_line, FONT_TITLE);
annotation(fig_ef, 'textbox', [TRACE_L, E_T5_y + ROW_H - 0.04, 0.06, 0.03], ...
    'String', 'control', 'FontSize', 6, 'FontWeight', 'bold', ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', 'Color', COL_T5.ctrl_line);
annotation(fig_ef, 'textbox', [TRACE_L + 0.06, E_T5_y + ROW_H - 0.04, 0.06, 0.03], ...
    'String', '{\ittutl-}', 'Interpreter', 'tex', 'FontSize', 6, 'FontWeight', 'bold', ...
    'FontName', 'Helvetica', 'EdgeColor', 'none', 'Color', COL_T5.ttl_line);

CARTOON_Y_T4 = 0.945;
CARTOON_Y_T5 = 0.895;
draw_stim_cartoon_strip(fig_ef, TRACE_L, CARTOON_Y_T4, TILE_W, 0, TILE_GAP, N_POS, 'bright', cartoon_orient);
draw_stim_cartoon_strip(fig_ef, TRACE_L, CARTOON_Y_T5, TILE_W, 0, TILE_GAP, N_POS, 'dark',   cartoon_orient);

annotation(fig_ef, 'textbox', [TRACE_L, E_T5_y - 0.06, 0.15, 0.04], ...
    'String', 'proximal / trailing', 'FontSize', FONT_LABEL, ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
annotation(fig_ef, 'textbox', [TRACE_R - 0.15, E_T5_y - 0.06, 0.15, 0.04], ...
    'String', 'distal / leading', 'FontSize', FONT_LABEL, ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

%% ===================== AMPLITUDE panels (right) ==========================

% T4 (ON) depolarization
ax = axes(fig_ef, 'Position', [AMP_T4_L, PD_DEP_Y, AMP_W, PD_DEP_H]);
dep_st_t4 = draw_amp_line(ax, positions, c_dep, t_dep, ...
    COL_T4.ctrl_line, COL_T4.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T4.ctrl_fill, COL_T4.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, dep_pp, dep_st_t4, STAT_COLOR, FONT_STAT);
format_amp(ax, YLIM_AMP_DEP, FONT_AX, false, false);
ylabel(ax, 'mV', 'FontSize', FONT_LABEL);
if opts.show_fwhm
    add_fwhm_bars(ax, positions, dep_st_t4, COL_T4.ctrl_line, COL_T4.ttl_line, ...
        YLIM_AMP_DEP, FONT_STAT, TTL_TEX, c_dep, t_dep, [diag_label ' T4'], opts.fwhm_method);
end

% T5 (OFF) depolarization
ax = axes(fig_ef, 'Position', [AMP_T5_L, PD_DEP_Y, AMP_W, PD_DEP_H]);
dep_st_t5 = draw_amp_line(ax, positions, off_c_dep, off_t_dep, ...
    COL_T5.ctrl_line, COL_T5.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T5.ctrl_fill, COL_T5.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, off_dep_pp, dep_st_t5, STAT_COLOR, FONT_STAT);
format_amp(ax, YLIM_AMP_DEP, FONT_AX, false, true);
if opts.show_fwhm
    add_fwhm_bars(ax, positions, dep_st_t5, COL_T5.ctrl_line, COL_T5.ttl_line, ...
        YLIM_AMP_DEP, FONT_STAT, TTL_TEX, off_c_dep, off_t_dep, [diag_label ' T5'], opts.fwhm_method);
end

% T4 (ON) hyperpolarization
ax = axes(fig_ef, 'Position', [AMP_T4_L, PD_HYP_Y, AMP_W, PD_HYP_H]);
hyp_st_t4 = draw_amp_line(ax, positions, c_hyp, t_hyp, ...
    COL_T4.ctrl_line, COL_T4.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T4.ctrl_fill, COL_T4.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, hyp_pp, hyp_st_t4, STAT_COLOR, FONT_STAT);
format_amp(ax, YLIM_HYP, FONT_AX, false, false);
ylabel(ax, 'mV', 'FontSize', FONT_LABEL);

% T5 (OFF) hyperpolarization
ax = axes(fig_ef, 'Position', [AMP_T5_L, PD_HYP_Y, AMP_W, PD_HYP_H]);
hyp_st_t5 = draw_amp_line(ax, positions, off_c_hyp, off_t_hyp, ...
    COL_T5.ctrl_line, COL_T5.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T5.ctrl_fill, COL_T5.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, off_hyp_pp, hyp_st_t5, STAT_COLOR, FONT_STAT);
format_amp(ax, YLIM_HYP, FONT_AX, true, true);

%% ===================== Diagnostic ========================================
fprintf('\nAsterisk summary (p<0.05 count):\n');
fprintf('  %s T4 dep:  %d sig\n', diag_label, sum(dep_pp     < 0.05));
fprintf('  %s T5 dep:  %d sig\n', diag_label, sum(off_dep_pp < 0.05));
fprintf('  %s T4 hyp:  %d sig\n', diag_label, sum(hyp_pp     < 0.05));
fprintf('  %s T5 hyp:  %d sig\n', diag_label, sum(off_hyp_pp < 0.05));

%% ===================== Labels & annotations ==============================
lbl_fsize = 12;

annotation(fig_ef, 'textbox', [0.00, E_T4_y + ROW_H - 0.02, 0.04, 0.06], ...
    'String', panel_letter_traces, 'FontSize', lbl_fsize, 'FontWeight', 'bold', ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

annotation(fig_ef, 'textbox', [AMP_T4_L - 0.06, PD_DEP_Y + PD_DEP_H - 0.02, 0.04, 0.06], ...
    'String', panel_letter_amp, 'FontSize', lbl_fsize, 'FontWeight', 'bold', ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

annotation(fig_ef, 'textbox', [TRACE_L, E_T4_y + ROW_H + 0.005, 0.30, 0.03], ...
    'String', section_label, 'FontSize', FONT_TITLE, ...
    'FontWeight', 'bold', 'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

annotation(fig_ef, 'textbox', [AMP_T4_L, PD_DEP_Y + PD_DEP_H + 0.005, AMP_W, 0.03], ...
    'String', 'T4 (ON)', 'FontSize', FONT_TITLE, ...
    'FontWeight', 'bold', 'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
annotation(fig_ef, 'textbox', [AMP_T5_L, PD_DEP_Y + PD_DEP_H + 0.005, AMP_W, 0.03], ...
    'String', 'T5 (OFF)', 'FontSize', FONT_TITLE, ...
    'FontWeight', 'bold', 'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%% ===================== Scale bar =========================================
ax_sb = findobj(fig_ef, 'Type', 'axes', 'Tag', 'E_T5_1');
if ~isempty(ax_sb)
    add_horiz_scale_bar(ax_sb(1), Y_LIM, FONT_AX, TRACE_START_DS, DS_FACTOR);
end

%% ===================== Export ============================================
ts = datestr(now, 'yyyymmdd_HHMM');
if USE_RAW_TRACES, raw_tag = '_absolute'; else, raw_tag = ''; end
if opts.show_fwhm, raw_tag = [raw_tag '_withFWHM']; end                 % tags mark departures from the manuscript defaults
if strcmpi(opts.pool_test, 'pooled'), raw_tag = [raw_tag '_pooledTest']; end
if strcmpi(opts.fwhm_method, 'gauss'), raw_tag = [raw_tag '_gaussFWHM']; end
pdf_file = fullfile(out_dir, sprintf('%s%s_%s.pdf', out_tag, raw_tag, ts));
png_file = fullfile(out_dir, sprintf('%s%s_%s.png', out_tag, raw_tag, ts));

if ~opts.skip_export
    exportgraphics(fig_ef, pdf_file, 'ContentType', 'vector');
    exportgraphics(fig_ef, png_file, 'Resolution', 300);
    fprintf('Saved: %s\n', pdf_file);
    fprintf('Saved: %s\n', png_file);
end
fprintf('=== generate_manuscript_fig_ef(%s) done ===\n', axis_mode);

end


%% =========================================================================
%%                          LOCAL FUNCTIONS
%% =========================================================================

function [ctrl_stats, ttl_stats] = compute_all_stats(traces_ctrl, traces_ttl, n_pos, stat_method, ds_factor)
    if nargin < 5, ds_factor = 1; end
    ctrl_stats = struct('center', cell(1, n_pos), 'spread', cell(1, n_pos), 'n', num2cell(zeros(1, n_pos)));
    ttl_stats  = struct('center', cell(1, n_pos), 'spread', cell(1, n_pos), 'n', num2cell(zeros(1, n_pos)));
    for pi = 1:n_pos
        [ctrl_stats(pi).center, ctrl_stats(pi).spread, ctrl_stats(pi).n] = stack_and_stat(traces_ctrl, pi, stat_method, ds_factor);
        [ttl_stats(pi).center,  ttl_stats(pi).spread,  ttl_stats(pi).n]  = stack_and_stat(traces_ttl,  pi, stat_method, ds_factor);
    end
end


function [center, spread, n] = stack_and_stat(traces_cell, pos_idx, stat_method, ds_factor)
    if nargin < 4, ds_factor = 1; end
    center = []; spread = []; n = 0;
    if isempty(traces_cell), return; end
    min_len = Inf;
    for k = 1:numel(traces_cell)
        mat = traces_cell{k};
        if ~isempty(mat) && pos_idx <= size(mat, 1)
            min_len = min(min_len, size(mat, 2));
        end
    end
    if isinf(min_len), return; end
    all_tr = [];
    for k = 1:numel(traces_cell)
        mat = traces_cell{k};
        if ~isempty(mat) && pos_idx <= size(mat, 1)
            row = mat(pos_idx, 1:min_len);
            if ~all(isnan(row)), all_tr = [all_tr; row]; end %#ok<AGROW>
        end
    end
    if isempty(all_tr), return; end
    if ds_factor > 1
        all_tr_ds = zeros(size(all_tr, 1), ceil(size(all_tr, 2) / ds_factor));
        for kk = 1:size(all_tr, 1)
            all_tr_ds(kk, :) = blockavg_downsample(all_tr(kk, :), ds_factor);
        end
        all_tr = all_tr_ds;
    end
    n = size(all_tr, 1);
    if strcmpi(stat_method, 'mean_sem')
        center = mean(all_tr, 1, 'omitnan');
        spread = std(all_tr, 0, 1, 'omitnan') ./ sqrt(n);
    else
        center = median(all_tr, 1, 'omitnan');
        spread = mad(all_tr, 1, 1);
    end
end


function draw_trace_row(fig, x_left, y_bottom, tile_w, tile_h, tile_gap, ...
    ctrl_stats, ttl_stats, n_pos, y_lim, col, line_w, ...
    ~, ~, font_size, label_font, stim_onset, stim_offset, trace_start, row_tag, y_label)

    if nargin < 21, y_label = '\DeltamV'; end

    for pi = 1:n_pos
        x = x_left + (pi-1) * (tile_w + tile_gap);
        ax = axes(fig, 'Position', [x, y_bottom, tile_w, tile_h]); %#ok<LAXES>
        if ~isempty(row_tag), ax.Tag = sprintf('%s_%d', row_tag, pi); end
        hold(ax, 'on');

        ctr_c = ctrl_stats(pi).center;  spr_c = ctrl_stats(pi).spread;
        ctr_t = ttl_stats(pi).center;   spr_t = ttl_stats(pi).spread;

        x_len = 0;
        if ~isempty(ctr_c), x_len = max(x_len, numel(ctr_c)); end
        if ~isempty(ctr_t), x_len = max(x_len, numel(ctr_t)); end
        if x_len == 0
            ylim(ax, y_lim); set(ax, 'XTick', [], 'YTick', []); box(ax, 'off'); continue;
        end

        stim_rect_bot = y_lim(1) + 0.06 * diff(y_lim);
        stim_rect_top = y_lim(1) + 0.14 * diff(y_lim);
        fill(ax, [stim_onset, stim_offset, stim_offset, stim_onset], ...
            [stim_rect_bot, stim_rect_bot, stim_rect_top, stim_rect_top], ...
            col.stim_line, 'EdgeColor', 'none', 'FaceAlpha', 1.0);

        xv = 1:x_len;
        if ~isempty(ctr_c) && ctrl_stats(pi).n >= 2
            fill(ax, [xv, fliplr(xv)], [ctr_c + spr_c, fliplr(ctr_c - spr_c)], ...
                col.ctrl_fill, 'FaceAlpha', col.alpha, 'EdgeColor', 'none');
            plot(ax, xv, ctr_c, '-', 'Color', col.ctrl_line, 'LineWidth', line_w);
        end
        if ~isempty(ctr_t) && ttl_stats(pi).n >= 2
            fill(ax, [xv, fliplr(xv)], [ctr_t + spr_t, fliplr(ctr_t - spr_t)], ...
                col.ttl_fill, 'FaceAlpha', col.alpha, 'EdgeColor', 'none');
            plot(ax, xv, ctr_t, '-', 'Color', col.ttl_line, 'LineWidth', line_w);
        end

        ylim(ax, y_lim);
        xlim(ax, [trace_start, x_len]);
        set(ax, 'XTick', []);
        if pi == 1
            ylabel(ax, y_label, 'FontSize', label_font);
            set(ax, 'FontSize', font_size, 'TickDir', 'out', 'TickLength', [0.02 0.02]);
            set(ax, 'YTick', -5:5:25);
        else
            set(ax, 'YTick', []);
        end
        box(ax, 'off');
        ax.XColor = 'none';
        ax.YColor = 'none';
        if pi == 1
            for yt = -5:5:25
                text(ax, trace_start - 0.02*(x_len - trace_start), yt, num2str(yt), ...
                    'FontSize', font_size, 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'middle');
            end
            line(ax, [trace_start trace_start], [-5 25], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
            for yt = -5:5:25
                line(ax, [trace_start trace_start - 0.01*(x_len-trace_start)], [yt yt], ...
                    'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
            end
        end
    end

    for pi = 1:(n_pos - 1)
        sep_x = x_left + pi * (tile_w + tile_gap) - tile_gap/2;
        annotation(fig, 'line', [sep_x sep_x], [y_bottom, y_bottom + tile_h], ...
            'Color', 'w', 'LineWidth', 2);
    end
end


function draw_stim_cartoon_strip(fig, x_left, y_bottom, tile_w, ~, tile_gap, n_pos, bar_type, bar_orient)
% bar_type: 'bright' (T4/ON) or 'dark' (T5/OFF)
% bar_orient: 'vertical' (PD axis) or 'horizontal' (ortho axis)
    GREEN_MID  = [0.15 0.55 0.15];
    GREEN_HIGH = [0.4  0.95 0.2];
    BAR_OFF    = [0.0  0.0  0.0];
    if strcmp(bar_type, 'bright'), bar_color = GREEN_HIGH; else, bar_color = BAR_OFF; end

    n_arena  = 15;
    bar_w    = 4;
    bar_frac = bar_w / n_arena;
    fig_sz = get(fig, 'PaperSize');
    sq_w   = tile_w * 0.30;
    sq_h   = sq_w * fig_sz(1) / fig_sz(2);

    for pi = 1:n_pos
        tile_x = x_left + (pi - 1) * (tile_w + tile_gap);
        sq_x   = tile_x + (tile_w - sq_w) / 2;
        ax     = axes(fig, 'Position', [sq_x, y_bottom, sq_w, sq_h]); %#ok<LAXES>
        hold(ax, 'on');
        patch(ax, [0 1 1 0], [0 0 1 1], GREEN_MID, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.3);
        bar_center = ((bar_w/2) + (pi-1) * (n_arena - bar_w) / 10) / n_arena;
        if strcmp(bar_orient, 'horizontal')
            by = bar_center - bar_frac/2;
            patch(ax, [0 1 1 0], [by by by+bar_frac by+bar_frac], bar_color, 'EdgeColor', 'none');
        else
            bx = bar_center - bar_frac/2;
            patch(ax, [bx bx+bar_frac bx+bar_frac bx], [0 0 1 1], bar_color, 'EdgeColor', 'none');
        end
        xlim(ax, [0 1]); ylim(ax, [0 1]);
        axis(ax, 'off');
    end
end


function add_row_label_colored(fig, x, y, cell_type, ~, font_size)
    annotation(fig, 'textbox', [x - 0.005, y - 0.06, 0.025, 0.12], ...
        'String', cell_type, ...
        'FontSize', font_size + 2, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica', 'EdgeColor', 'none', 'Color', [0 0 0], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Rotation', 90);
end


function add_horiz_scale_bar(ax, y_lim, font_size, trace_start, ds_factor)
    if nargin < 5, ds_factor = 1; end
    v_range = diff(y_lim);
    bar_t = 5000 / ds_factor;
    x_left   = trace_start + 0.05 * (xlim(ax) * [0;1] - trace_start);
    y_bottom = y_lim(1) + 0.08 * v_range;
    plot(ax, [x_left, x_left + bar_t], [y_bottom, y_bottom], 'k-', ...
        'LineWidth', 0.8, 'Clipping', 'off');
    text(ax, x_left, y_bottom - 0.06*v_range, '0.5 s', ...
        'HorizontalAlignment', 'left', 'FontSize', font_size, 'Clipping', 'off');
end


function [dep_mat, hyp_mat] = extract_robust_amps(results_sub, trace_field, ...
    dep_win, hyp_win, dep_pct, hyp_pct, ~)
    n = numel(results_sub);
    dep_mat = NaN(n, 11); hyp_mat = NaN(n, 11);
    for k = 1:n
        if ~isfield(results_sub(k), trace_field), continue; end
        traces = results_sub(k).(trace_field);
        if isempty(traces), continue; end
        n_samp = size(traces, 2);
        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end
            de = min(dep_win(2), n_samp);
            if ~isinf(de)
                dep_mat(k, pos) = prctile(row(dep_win(1):de), dep_pct);
            else
                dep_mat(k, pos) = prctile(row(dep_win(1):end), dep_pct);
            end
            he = hyp_win(2); if isinf(he), he = n_samp; end
            he = min(he, n_samp);
            hval = prctile(row(hyp_win(1):he), hyp_pct);
            hyp_mat(k, pos) = min(hval, 0);
        end
    end
end


function stats = draw_amp_line(ax, positions, ctrl_data, ttl_data, ...
    col_c, col_t, line_w, marker_sz, ~, col_cf, col_tf, alpha_fill)
    hold(ax, 'on');
    n_pos = numel(positions);
    stats = struct('pos', num2cell(positions), 'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, 'sem_ctrl', NaN, 'sem_ttl', NaN);

    ctrl_mean = NaN(1, n_pos); ctrl_sem = NaN(1, n_pos);
    ttl_mean  = NaN(1, n_pos); ttl_sem  = NaN(1, n_pos);

    for i = 1:n_pos
        cv = ctrl_data(:, i); cv = cv(~isnan(cv));
        tv = ttl_data(:, i);  tv = tv(~isnan(tv));
        stats(i).n_ctrl = numel(cv); stats(i).n_ttl = numel(tv);
        if numel(cv) >= 2
            ctrl_mean(i) = mean(cv); ctrl_sem(i) = std(cv)/sqrt(numel(cv));
        elseif numel(cv) == 1
            ctrl_mean(i) = cv; ctrl_sem(i) = 0;
        end
        stats(i).mean_ctrl = ctrl_mean(i); stats(i).sem_ctrl = ctrl_sem(i);
        if numel(tv) >= 2
            ttl_mean(i) = mean(tv); ttl_sem(i) = std(tv)/sqrt(numel(tv));
        elseif numel(tv) == 1
            ttl_mean(i) = tv; ttl_sem(i) = 0;
        end
        stats(i).mean_ttl = ttl_mean(i); stats(i).sem_ttl = ttl_sem(i);
    end

    valid_c = ~isnan(ctrl_mean); valid_t = ~isnan(ttl_mean);
    if any(valid_c)
        xc = positions(valid_c);
        fill(ax, [xc, fliplr(xc)], ...
            [ctrl_mean(valid_c)+ctrl_sem(valid_c), fliplr(ctrl_mean(valid_c)-ctrl_sem(valid_c))], ...
            col_cf, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    end
    if any(valid_t)
        xt = positions(valid_t);
        fill(ax, [xt, fliplr(xt)], ...
            [ttl_mean(valid_t)+ttl_sem(valid_t), fliplr(ttl_mean(valid_t)-ttl_sem(valid_t))], ...
            col_tf, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    end

    if any(valid_c)
        plot(ax, positions(valid_c), ctrl_mean(valid_c), '-o', ...
            'Color', col_c, 'MarkerFaceColor', col_c, 'MarkerSize', marker_sz, 'LineWidth', line_w);
    end
    if any(valid_t)
        plot(ax, positions(valid_t), ttl_mean(valid_t), '-o', ...
            'Color', col_t, 'MarkerFaceColor', col_t, 'MarkerSize', marker_sz, 'LineWidth', line_w);
    end
end


function format_amp(ax, y_lim, font_ax, show_xlabel, hide_yticklabels)
    xline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    xlim(ax, [-5.5 5.5]); ylim(ax, y_lim);
    set(ax, 'FontSize', font_ax, 'TickDir', 'out', 'LineWidth', 0.4, 'XTick', -4:2:4);
    box(ax, 'off');
    if show_xlabel
        xlabel(ax, 'Position', 'FontSize', font_ax);
    else
        set(ax, 'XTickLabel', []);
    end
    if hide_yticklabels, set(ax, 'YTickLabel', []); end
end


function [pool_pvals, pool_centers] = compute_pooled_ranksum(ctrl_mat, ttl_mat, pool_test)
% 'pooled'  : all cell x position values of the three positions (manuscript)
% 'percell' : each cell's mean over the three positions, one value per cell
    if nargin < 3, pool_test = 'pooled'; end
    pool_centers = -4:4;
    pool_pvals = NaN(1, 9);
    for p = 1:9
        cols = p:p+2;
        c = ctrl_mat(:, cols);
        t = ttl_mat(:, cols);
        if strcmpi(pool_test, 'percell')
            c = mean(c, 2, 'omitnan');
            t = mean(t, 2, 'omitnan');
        end
        c = c(:); c = c(~isnan(c));
        t = t(:); t = t(~isnan(t));
        if numel(c) >= 2 && numel(t) >= 2
            pool_pvals(p) = ranksum(c, t);
        end
    end
end


function draw_pooled_asterisks(ax, pool_centers, pool_pvals, stats_11, col_stat, font_stat)
    pos11 = -5:5;
    yl = ylim(ax);
    is_hyp = yl(2) <= 0;
    for p = 1:numel(pool_centers)
        if isnan(pool_pvals(p)) || pool_pvals(p) >= 0.05, continue; end
        x = pool_centers(p);
        idx = find(pos11 == x, 1);
        if isempty(idx), continue; end
        y_c = stats_11(idx).mean_ctrl;  sem_c = stats_11(idx).sem_ctrl;
        y_t = stats_11(idx).mean_ttl;   sem_t = stats_11(idx).sem_ttl;
        if isnan(y_c) && isnan(y_t), continue; end
        if is_hyp
            y_bot  = min([y_c - sem_c, y_t - sem_t], [], 'omitnan');
            y_star = y_bot - 0.08 * diff(yl);
            va = 'top';
        else
            y_top  = max([y_c + sem_c, y_t + sem_t], [], 'omitnan');
            y_star = y_top + 0.06 * diff(yl);
            y_star = min(y_star, yl(2) - 0.02 * diff(yl));
            va = 'bottom';
        end
        if pool_pvals(p) < 0.001, str = '***';
        elseif pool_pvals(p) < 0.01, str = '**';
        else, str = '*'; end
        text(ax, x, y_star, str, 'FontSize', font_stat, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', va, ...
            'Color', col_stat, 'FontWeight', 'bold', 'Clipping', 'off');
    end
end


function [fw, left_x, right_x] = compute_fwhm_positions(positions, mean_vals)
% Linear-interpolated half-max crossings. If a flank never falls to
% half-max within the mapped positions the FWHM is undefined (truncated)
% and NaN is returned, so truncated cells drop out of the per-cell test.
    fw = NaN; left_x = NaN; right_x = NaN;
    valid = ~isnan(mean_vals);
    if sum(valid) < 3, return; end
    pos_v = positions(valid); val_v = mean_vals(valid);
    [pk, pk_idx] = max(val_v);
    if pk <= 0, return; end
    half_max = pk / 2;
    for j = pk_idx:-1:2
        if val_v(j-1) <= half_max
            frac = (half_max - val_v(j-1)) / (val_v(j) - val_v(j-1));
            left_x = pos_v(j-1) + frac * (pos_v(j) - pos_v(j-1)); break;
        end
    end
    if isnan(left_x), return; end
    for j = pk_idx:numel(val_v)-1
        if val_v(j+1) <= half_max
            frac = (half_max - val_v(j+1)) / (val_v(j) - val_v(j+1));
            right_x = pos_v(j+1) - frac * (pos_v(j+1) - pos_v(j)); break;
        end
    end
    if isnan(right_x), left_x = NaN; return; end
    fw = right_x - left_x;
end


function [fw, left_x, right_x] = compute_fwhm_gauss(positions, vals)
% FWHM of a single Gaussian A*exp(-(x-mu)^2/(2 sigma^2)) fitted (least squares,
% fminsearch) to the amplitude-by-position profile; needs >= 4 valid points.
% Returns NaN if the fit runs to an implausible width (sigma > 10 or < 0.3 grid
% positions, i.e. FWHM > ~59 deg or < ~1.8 deg).
    fw = NaN; left_x = NaN; right_x = NaN;
    ok = ~isnan(vals);
    if sum(ok) < 4, return; end
    x = positions(ok); y = vals(ok);
    [pk, i_pk] = max(y);
    if pk <= 0, return; end
    cost = @(p) sum((y - p(1) * exp(-(x - p(2)).^2 / (2 * p(3)^2))).^2);
    p = fminsearch(cost, [pk, x(i_pk), 2], optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6));
    sigma = abs(p(3));
    if sigma > 10 || sigma < 0.3, return; end
    fw = 2 * sqrt(2 * log(2)) * sigma;
    left_x = p(2) - fw / 2; right_x = p(2) + fw / 2;
end


function [fw, left_x, right_x] = compute_fwhm_any(positions, vals, method)
    if strcmpi(method, 'gauss')
        [fw, left_x, right_x] = compute_fwhm_gauss(positions, vals);
    else
        [fw, left_x, right_x] = compute_fwhm_positions(positions, vals);
    end
end


function fwhm_vec = compute_percell_fwhm(positions, amp_mat, method)
    if nargin < 3, method = 'interp'; end
    n = size(amp_mat, 1); fwhm_vec = NaN(n, 1);
    for k = 1:n
        [fw, ~, ~] = compute_fwhm_any(positions, amp_mat(k, :), method);
        fwhm_vec(k) = fw;
    end
end


function add_fwhm_bars(ax, positions, stats, col_c, col_t, y_lim, font_stat, ~, ...
    ctrl_dep_mat, ttl_dep_mat, label, method)
    if nargin < 11, label = ''; end
    if nargin < 12, method = 'interp'; end
    DEG_PER_POS = 2.5;   % 2-px grid, 1.25 deg/px
    ctrl_mn = [stats.mean_ctrl]; ttl_mn = [stats.mean_ttl];
    [fw_c, lx_c, rx_c] = compute_fwhm_any(positions, ctrl_mn, method);
    [fw_t, lx_t, rx_t] = compute_fwhm_any(positions, ttl_mn, method);
    fwhm_y_base = y_lim(1) + 0.5; fwhm_y_gap = 1.0;
    if ~isnan(fw_c)
        y_c = fwhm_y_base;
        plot(ax, [lx_c rx_c], [y_c y_c], '-', 'Color', col_c, 'LineWidth', 2.5);
        plot(ax, [lx_c lx_c], y_c+[-0.3 0.3], '-', 'Color', col_c, 'LineWidth', 1.0);
        plot(ax, [rx_c rx_c], y_c+[-0.3 0.3], '-', 'Color', col_c, 'LineWidth', 1.0);
    end
    if ~isnan(fw_t)
        y_t = fwhm_y_base + fwhm_y_gap;
        plot(ax, [lx_t rx_t], [y_t y_t], '-', 'Color', col_t, 'LineWidth', 2.5);
        plot(ax, [lx_t lx_t], y_t+[-0.3 0.3], '-', 'Color', col_t, 'LineWidth', 1.0);
        plot(ax, [rx_t rx_t], y_t+[-0.3 0.3], '-', 'Color', col_t, 'LineWidth', 1.0);
    end
    cv_fw = compute_percell_fwhm(positions, ctrl_dep_mat, method);
    tv_fw = compute_percell_fwhm(positions, ttl_dep_mat, method);
    cv_fw = cv_fw(~isnan(cv_fw)); tv_fw = tv_fw(~isnan(tv_fw));
    p_str = '';
    p_fw = NaN;
    if numel(cv_fw) >= 2 && numel(tv_fw) >= 2
        p_fw = ranksum(cv_fw, tv_fw);
        if p_fw < 0.001,     p_str = ' ***';
        elseif p_fw < 0.01,  p_str = ' **';
        elseif p_fw < 0.05,  p_str = ' *';
        end
    end
    text(ax, 0, fwhm_y_base + fwhm_y_gap + 1.5, ...
        ['FWHM' p_str], ...
        'FontSize', font_stat, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    fprintf(['FWHM (%s) %s: group-mean curve ctrl %.1f deg, tutl %.1f deg; per-cell ' ...
             'median ctrl %.1f deg (n=%d), tutl %.1f deg (n=%d), rank-sum p = %.3f ' ...
             '(cells without a defined width excluded)\n'], method, label, ...
        fw_c * DEG_PER_POS, fw_t * DEG_PER_POS, ...
        median(cv_fw) * DEG_PER_POS, numel(cv_fw), ...
        median(tv_fw) * DEG_PER_POS, numel(tv_fw), p_fw);
end


function y_ds = blockavg_downsample(y, factor)
    n = numel(y);
    n_out = ceil(n / factor);
    y_ds  = zeros(1, n_out);
    for bi = 1:n_out
        i1 = (bi-1)*factor + 1;
        i2 = min(bi*factor, n);
        y_ds(bi) = mean(y(i1:i2), 'omitnan');
    end
end
