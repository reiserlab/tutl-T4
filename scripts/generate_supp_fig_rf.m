function out = generate_supp_fig_rf(opts)
% GENERATE_SUPP_FIG_RF  Figure S5A-C: receptive fields from 4-px square flashes.
%
%   OUT = GENERATE_SUPP_FIG_RF(OPTS) builds, for every recorded T4 / T5 cell (late batch
%   in batch_results.mat plus the early, sweep-only batch), the signed 4-px square-flash
%   map, fits a 2-D Gaussian with axes fixed to the cell's preferred direction, aligns
%   every cell on its fitted centre in the PD frame (u along PD, proximal side negative;
%   v orthogonal) and exports one figure (18 x 10 cm):
%
%     row 1  group-mean signed maps (T4 ctrl, T4 tutl-, T5 ctrl, T5 tutl-) with the
%            group-mean half-maximum footprint, OD line, PD arrow and RF centre   (S5A)
%     row 2  T4 and T5 half-maximum footprints, ctrl vs tutl-, mean +- SEM       (S5B)
%            per-cell FWHM along PD-ND and orthogonal, proximal and distal
%            half-widths, box + dots, rank-sum brackets                            (S5C)
%
%   Per-square value ("clean" signed map): each square's 3-repeat mean trace (from
%   PARSE_FLASH_DATA, 8th output) is referenced to its own pre-flash baseline (mean of
%   the 100 ms before onset); early window = mean voltage change over 0..260 ms, late
%   window = mean(150..500 ms) kept only if negative (t = 0 is the first sample of the flash
%   frame, parse_flash_data sample 1002); the component with the larger magnitude is
%   kept; |value| < 1 mV is set to 0. The 2-D Gaussian fit and all RF numbers use the
%   depolarisation map (98th percentile, PARSE_FLASH_DATA max_data) with an offset term.
%   PD for the fit = continuous vector-sum direction of the bar-sweep responses
%   (pd_direction_vs in both batch files; not snapped to the stimulus grid). Footprints: per-cell half-maximum contour built from the split-fit
%   proximal / distal widths and the 2-D-fit orthogonal width, centred on the 2-D-fit
%   centre; line = pointwise mean across cells, band = pointwise SEM. A cell enters every
%   panel only if both fits converged (n in the titles = contributing cells).
%   Statistics: two-sided Wilcoxon rank-sum on per-cell values, uncorrected.
%
%   Outputs: manuscript_figures/fig_S5_rf_squares_<timestamp>.pdf/.png and a .mat with
%   every per-cell value (OUT).
%
%   OPTS fields (optional): .data_root ('/Users/reiserm/Documents/ttl_1DRF'),
%   .skip_export (false).
%
%   See also PARSE_FLASH_DATA, GENERATE_MANUSCRIPT_FIG_EF.

if nargin < 1, opts = struct(); end
if ~isfield(opts, 'data_root'),   opts.data_root   = '/Users/reiserm/Documents/ttl_1DRF'; end
if ~isfield(opts, 'skip_export'), opts.skip_export = false; end
data_root = opts.data_root;
out_dir   = fullfile(data_root, 'manuscript_figures');
if ~isfolder(out_dir) && ~opts.skip_export, mkdir(out_dir); end

DEG_PER_PX = 1.25;                          % arena pixel pitch (deg)
FWHM_C     = 2 * sqrt(2 * log(2));          % FWHM / sigma
CQ         = sqrt(2 * log(2));              % half-maximum radius / sigma
U          = -14:14;                        % aligned PD-frame grid (px), 1-px resolution
CLEAN_DEAD = 1.0;                           % clean maps: |value| below this (mV) set to 0
XL = [-17.25 10.75]; YL = [-12 12];         % map window (deg), one square further proximal than distal
COLS  = {[0 0 0], [1 0 0], [0.4 0.4 0.4], [0.8 0.2 0.2]};
GNAME = {'T4 ctrl', 'T4 tutl-', 'T5 ctrl', 'T5 tutl-'};

%% ===================== Cell list =========================================
S = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
R = S.results;
C = struct('path', {}, 'name', {}, 'group', {}, 'is_on', {}, 'is_ttl', {}, 'pd_direction', {}, 'filter_idx', {});
assert(isfield(R, 'pd_direction_vs'), 'generate_supp_fig_rf:OutdatedBatchFile', 'batch_results.mat lacks pd_direction_vs; rebuild with build_batch_results.m');
for k = 1:numel(R)
    C(end+1) = struct('path', fullfile(data_root, R(k).folder), 'name', R(k).folder, 'group', R(k).group, ...
        'is_on', R(k).is_on, 'is_ttl', R(k).is_ttl, 'pd_direction', R(k).pd_direction_vs, 'filter_idx', true); %#ok<AGROW>
end
early_root = fullfile(data_root, 'pre-bar-flash');
Se = load(fullfile(early_root, 'population_results', 'batch_results_pre_bf.mat'), 'results');
Re = Se.results;
assert(isfield(Re, 'pd_direction_vs'), 'generate_supp_fig_rf:OutdatedBatchFile', 'batch_results_pre_bf.mat lacks pd_direction_vs; rebuild with build_batch_results.m');
for k = 1:numel(Re)
    [~, nm] = fileparts(Re(k).folder);
    if Re(k).is_ttl, s1 = 'ttl'; else, s1 = 'control'; end
    if Re(k).is_on,  s2 = 'ON';  else, s2 = 'OFF'; end
    p = fullfile(early_root, s1, s2, nm);
    if ~isfolder(p), warning('early session folder not found: %s', p); continue; end
    C(end+1) = struct('path', p, 'name', char(nm), 'group', char(Re(k).group), 'is_on', Re(k).is_on, ...
        'is_ttl', Re(k).is_ttl, 'pd_direction', Re(k).pd_direction_vs, 'filter_idx', false); %#ok<AGROW>
end
n_cells = numel(C);
fprintf('%d cells\n', n_cells);

%% ===================== Per-cell maps, fits, alignment ====================
nU = numel(U); [uu, vv] = meshgrid(U, U);
fitp = NaN(n_cells, 6); r2 = NaN(n_cells, 1);
fwhm_pd = NaN(n_cells, 1); fwhm_or = NaN(n_cells, 1);
hw_prox = NaN(n_cells, 1); hw_dist = NaN(n_cells, 1); split_p = NaN(n_cells, 5);
Mal  = NaN(nU, nU, n_cells);                % aligned depolarisation maps, offset removed (rows v, cols u)
MalC = NaN(nU, nU, n_cells);                % aligned clean signed maps
ok_cell = false(n_cells, 1);
orig_dir = pwd; cleanup = onCleanup(@() cd(orig_dir)); %#ok<NASGU>
for k = 1:n_cells
    c = C(k);
    fprintf('[%2d/%d] %s %s\n', k, n_cells, c.name, c.group);
    try
        [~, ~, Log] = load_protocol2_data(c.path); cd(orig_dir);
    catch ME
        cd(orig_dir); warning('  load failed: %s', ME.message); continue;
    end
    assert(isscalar(c.pd_direction) && isfinite(c.pd_direction), 'generate_supp_fig_rf:BadPD', '%s: pd_direction_vs is not a finite scalar', c.name);
    assert(size(Log.ADC.Volts, 1) >= 2 && all(isfinite(Log.ADC.Volts(2, :))), 'generate_supp_fig_rf:BadRecording', '%s: voltage channel missing or non-finite', c.name);
    f_data = Log.ADC.Volts(1, :); v_data = Log.ADC.Volts(2, :) * 10;   % channel 1 frame number, channel 2 Vm (x10 gain); 10 kHz
    if c.is_on, on_off = "on"; else, on_off = "off"; end
    try
        [~, ~, ~, ~, ~, mx, ~, mt] = parse_flash_data(f_data, v_data, on_off, "slow", 4, tempdir, c.filter_idx);
    catch ME
        warning('  parse_flash_data failed: %s', ME.message); continue;
    end
    [X, Y] = square_centroids(c.path, on_off);          % arena px of every square (x right, y up)
    Sc = clean_signed_map(mt, CLEAN_DEAD);

    % 2-D Gaussian, axes fixed to the PD frame: p = [A x0 y0 sigma_PD sigma_orth B].
    % PD = continuous vector-sum direction of the bar-sweep responses (pd_direction_vs,
    % stored by both batch scripts), not the direction snapped to the 22.5-deg stimulus grid.
    phi = deg2rad(c.pd_direction);
    [p, r2(k), fit_ok] = fit_gauss2d(X, Y, mx, phi);
    fitp(k, :) = p;
    fwhm_pd(k) = FWHM_C * p(4) * DEG_PER_PX;
    fwhm_or(k) = FWHM_C * p(5) * DEG_PER_PX;

    % align on the fitted centre in the PD frame (u < 0 = proximal)
    xs = p(2) - uu * cos(phi) + vv * sin(phi);
    ys = p(3) - uu * sin(phi) - vv * cos(phi);
    Mal(:, :, k)  = interp2(X, Y, mx - p(6), xs, ys, 'linear', NaN);
    MalC(:, :, k) = interp2(X, Y, Sc,        xs, ys, 'linear', NaN);

    % proximal vs distal half-widths: split Gaussian along the PD line through the centre
    [hw_prox(k), hw_dist(k), split_p(k, :), split_ok] = fit_split_gauss(U * DEG_PER_PX, Mal(U == 0, :, k));

    % a cell counts only if both fits converged to finite parameters; every panel and
    % statistic below uses the same ok_cell set, so displayed n = cells contributing
    ok_cell(k) = fit_ok && split_ok && all(isfinite(p)) && isfinite(hw_prox(k)) && isfinite(hw_dist(k));
    if ~ok_cell(k), warning('  %s: fit did not converge; cell excluded from all panels', c.name); end
end
fprintf('%d of %d cells parsed and fitted\n', sum(ok_cell), n_cells);

%% ===================== Group statistics ==================================
is_on = [C.is_on]'; is_ttl = [C.is_ttl]';
masks = {ok_cell & is_on & ~is_ttl, ok_cell & is_on & is_ttl, ok_cell & ~is_on & ~is_ttl, ok_cell & ~is_on & is_ttl};
measures = { ...
    fwhm_pd, 'FWHM along PD-ND (deg)',      [0 24], 0:6:24; ...
    fwhm_or, 'FWHM along orthogonal (deg)', [0 24], 0:6:24; ...
    hw_prox, 'proximal half-width (deg)',   [0 14], 0:2:14; ...
    hw_dist, 'distal half-width (deg)',     [0 14], 0:2:14};
fprintf('\n=== 4-px square-flash RF, %d cells (%s); median per group, rank-sum ctrl vs tutl- ===\n', sum(ok_cell), ...
    strjoin(arrayfun(@(g) sprintf('%s n = %d', GNAME{g}, sum(masks{g})), 1:4, 'uni', 0), ', '));
p_meas = NaN(size(measures, 1), 2);                                   % rank-sum p, ctrl vs tutl-, [T4 T5]
for s = 1:size(measures, 1)
    v = measures{s, 1}; fprintf('%-30s', measures{s, 2});
    for t = 1:2
        a = v(masks{2*t-1}); b = v(masks{2*t});
        p_meas(s, t) = ranksum(a, b);
        fprintf('   %s %6.2f vs %6.2f  p = %.3f', GNAME{2*t-1}(1:2), median(a), median(b), p_meas(s, t));
    end
    fprintf('\n');
end

%% ===================== Figure ============================================
th   = linspace(0, 2 * pi, 361);
Udeg = U * DEG_PER_PX;
cl_map = pop_clim(MalC, masks);
FP = cell(1, 4); MM = cell(1, 4);                                      % per-group: contours (cells x angles), masked mean map
for g = 1:4
    FP{g} = cell_footprints(th, split_p, fitp, masks{g}, DEG_PER_PX, CQ);
    MM{g} = group_mean_map(MalC(:, :, masks{g}));
end
jitter_rs = RandStream('mt19937ar', 'Seed', 0);                        % reproducible dot jitter, caller's RNG untouched
fig = new_fig(18, 10);
ry = [0.60 0.13];
% row 1: group-mean maps (S5A)
for g = 1:4
    ax = axes(fig, 'Position', [0.05 + (g-1) * 0.215, ry(1), 0.18, 0.34]); hold(ax, 'on'); %#ok<LAXES>
    draw_group_map(ax, MM{g}, cl_map, Udeg, XL, YL, mean(FP{g}, 1), th);
    title(ax, sprintf('%s mean (n = %d)', GNAME{g}, sum(masks{g})), 'FontWeight', 'normal', 'FontSize', 6.5);
    if g == 1, ylabel(ax, 'orthogonal (deg)'); end
    xlabel(ax, 'PD-ND axis (deg)');
end
cb = colorbar(ax, 'eastoutside'); cb.Position = [0.925 ry(1) + 0.03 0.012 0.28];
cb.Ticks = unique([cl_map(1), 0:5:cl_map(2), cl_map(2)]);                          % both limits labelled
annotation(fig, 'textbox', [0.91 ry(1) + 0.31 0.04 0.03], 'String', 'mV', 'FontSize', 6.5, 'EdgeColor', 'none', ...
    'FontName', 'Helvetica', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
% row 2: footprints (S5B) and per-cell widths (S5C)
for t = 1:2
    ax = axes(fig, 'Position', [0.05 + (t-1) * 0.20, ry(2), 0.16, 0.32]); %#ok<LAXES>
    draw_footprint_panel(ax, t, FP, masks, th, COLS, GNAME, t == 1);
end
box_x = [0.50 0.63 0.76 0.89];
for b = 1:size(measures, 1)
    ax = axes(fig, 'Position', [box_x(b), ry(2), 0.085, 0.32]); %#ok<LAXES>
    draw_box_panel(ax, measures{b, 1}, masks, p_meas(b, :), measures{b, 2}, measures{b, 3}, measures{b, 4}, COLS, jitter_rs);
end

%% ===================== Export / output ===================================
ts_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
if ~opts.skip_export
    pdf_file = fullfile(out_dir, sprintf('fig_S5_rf_squares_%s.pdf', ts_str));
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
    exportgraphics(fig, strrep(pdf_file, '.pdf', '.png'), 'Resolution', 300);
    fprintf('Saved: %s (+ .png)\n', pdf_file);
end
out = struct('cells', C, 'ok', ok_cell, 'fit_A_x0_y0_sPD_sOrth_B', fitp, 'R2', r2, 'fwhm_pd_deg', fwhm_pd, ...
    'fwhm_orth_deg', fwhm_or, 'halfwidth_prox_deg', hw_prox, 'halfwidth_dist_deg', hw_dist, ...
    'split_fit_A_u0_sProx_sDist_B', split_p, 'maps_aligned_depol', Mal, 'maps_aligned_signed_clean', MalC, 'U_px', U);
if ~opts.skip_export, save(fullfile(out_dir, sprintf('fig_S5_rf_squares_%s.mat', ts_str)), 'out'); end
fprintf('=== generate_supp_fig_rf done ===\n');
end


%% ========================================================================
%  Analysis helpers
%% ========================================================================
function [X, Y] = square_centroids(exp_folder, on_off)
% Arena-pixel centroid (x = column, y = -row so that up is positive) of the 4-px square
% shown in each cell of PARSE_FLASH_DATA's 14 x 14 map. Map index (row, col) of a flash
% with frame value F follows parse_flash_data: fnum = F - 1; rows = n - mod(fnum - off, n);
% cols = floor((fnum - off) / n) + 1, with off = n_flashes for ON and 0 for OFF.
% Frame value F is stored in Pats(:, :, F + 1) (frame 1 of the pattern is blank).
    n = 14; n_fl = n * n;
    d = dir(fullfile(exp_folder, 'Patterns', '*4px_square*.mat'));
    assert(numel(d) == 1, 'square_centroids:PatternFile', '%s: expected one 4px_square pattern file, found %d', exp_folder, numel(d));
    S = load(fullfile(d(1).folder, d(1).name), 'pattern');
    P = S.pattern.Pats;
    assert(size(P, 3) >= 2 * n_fl + 1, 'square_centroids:PatternFrames', '%s: 4px_square pattern has %d frames, expected >= %d', exp_folder, size(P, 3), 2 * n_fl + 1);
    if on_off == "on", off = n_fl; else, off = 0; end
    X = NaN(n, n); Y = NaN(n, n);
    for row = 1:n
        for col = 1:n
            fnum = off + (col - 1) * n + (n - row);
            fr = double(P(:, :, fnum + 2));
            bg = mode(fr(:)); [rr, cc] = find(fr ~= bg);
            assert(numel(rr) == 16, 'square_centroids:SquareSize', '%s: frame %d lights %d pixels, expected 16', exp_folder, fnum + 1, numel(rr));
            X(row, col) = mean(cc); Y(row, col) = -mean(rr);
        end
    end
end


function S = clean_signed_map(M, dead)
% Signed per-square value from the 3-repeat mean traces returned by parse_flash_data
% (rows x cols x samples; 10 kHz). parse_flash_data extracts each trace from 1000 samples
% before the last pre-flash frame sample, so the FIRST sample of the flash frame is 1002
% (t = 0); the trace spans -100.1..+599.8 ms; flash 160 ms.
%   baseline          = mean(-100..0 ms)   = samples    1..1001
%   early window      = mean(0..260 ms)    = samples 1002..3601, signed (mean voltage change)
%   late window       = mean(150..500 ms)  = samples 2502..6001, kept only if negative
% The component with the larger magnitude is kept; |value| < dead (mV) -> 0.
    ON = 1002;                                      % first sample of the flash frame
    base = mean(M(:, :, 1:ON-1), 3);
    S = mean(M(:, :, ON:ON+2599), 3) - base;
    hyp = min(mean(M(:, :, ON+1500:ON+4999), 3) - base, 0);
    use_h = abs(hyp) > abs(S); S(use_h) = hyp(use_h);
    S(abs(S) < dead) = 0;
end


function [p, r2, ok] = fit_gauss2d(X, Y, Z, phi)
% 2-D Gaussian with offset whose axes are fixed to the PD frame (angle phi):
% p = [A x0 y0 sigma_PD sigma_orth B]; bounded least squares (sigma 1..12 px, amplitude
% 0..80 mV, centre within the mapped area), best of 3 starting widths (3, 5, 8 px).
% OK is false if no start converged (lsqcurvefit exit flag <= 0).
    xy = [X(:), Y(:)]; z = Z(:);
    [zmax, imax] = max(z);
    lb = [0,   min(X(:)), min(Y(:)), 1,  1,  min(z) - 5];
    ub = [80,  max(X(:)), max(Y(:)), 12, 12, max(z)];
    o = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxFunctionEvaluations', 4000);
    f = @(p, xy) gauss2d_eval(p, xy, phi);
    best = Inf; p = NaN(1, 6); ok = false;
    for s0 = [3 5 8]
        p0 = [zmax - median(z), xy(imax, 1), xy(imax, 2), s0, s0, median(z)];
        [pp, res, ~, flag] = lsqcurvefit(f, p0, xy, z, lb, ub, o);
        if flag > 0 && res < best, best = res; p = pp; ok = true; end
    end
    r2 = 1 - best / sum((z - mean(z)).^2);
end


function z = gauss2d_eval(p, xy, phi)
% p = [A x0 y0 sigma_PD sigma_orth B]; xy = [x y] columns; phi = PD direction (rad).
    dx = xy(:, 1) - p(2);
    dy = xy(:, 2) - p(3);
    a  =  dx * cos(phi) + dy * sin(phi);        % along PD
    b  = -dx * sin(phi) + dy * cos(phi);        % orthogonal
    z  = p(1) * exp(-(a.^2 / (2 * p(4)^2) + b.^2 / (2 * p(5)^2))) + p(6);
end


function [hw_prox, hw_dist, p, ok] = fit_split_gauss(u_deg, prof)
% Split Gaussian along the PD line through the aligned centre (v = 0):
% A*exp(-(u-u0)^2/(2 s^2)) + B with s = s_prox for u < u0 (proximal side, u negative) and
% s = s_dist for u >= u0; free centre u0; bounded least squares (sigma 1.5..15 deg,
% amplitude 0..80 mV), best of 2 starting widths (4, 7 deg). Returns the half-widths at
% half-maximum (1.177 sigma, deg) and p = [A u0 s_prox s_dist B]. OK is false when the
% profile has fewer than 7 finite points, no positive peak, or no start converged.
    hw_prox = NaN; hw_dist = NaN; p = NaN(1, 5); ok = false;
    fin = isfinite(prof); u = u_deg(fin); y = prof(fin);
    if sum(fin) < 7, return; end
    [pk, ipk] = max(y);
    if pk <= 0, return; end
    f = @(p, u) p(1) * exp(-(u - p(2)).^2 ./ (2 * (p(3) * (u < p(2)) + p(4) * (u >= p(2))).^2)) + p(5);
    lb = [0,   min(u), 1.5, 1.5, min(y) - 5];
    ub = [80,  max(u), 15,  15,  max(y)];
    o = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxFunctionEvaluations', 4000);
    best = Inf;
    for s0 = [4 7]
        [pp, res, ~, flag] = lsqcurvefit(f, [pk, u(ipk), s0, s0, 0], u, y, lb, ub, o);
        if flag > 0 && res < best, best = res; p = pp; ok = true; end
    end
    c = sqrt(2 * log(2));
    hw_prox = c * p(3); hw_dist = c * p(4);
end


function R = cell_footprints(th, split_p, fitp, mask, deg_per_px, cq)
% Half-maximum contour radius of each cell in MASK (rows) at the angles TH (columns),
% from its own split-fit sigma_prox / sigma_dist (deg) and 2-D-fit sigma_orth (px -> deg).
    idx = find(mask); R = NaN(numel(idx), numel(th));
    for k = 1:numel(idx)
        R(k, :) = footprint_radius(th, split_p(idx(k), 3), split_p(idx(k), 4), fitp(idx(k), 5) * deg_per_px, cq);
    end
    assert(all(isfinite(R(:))), 'generate_supp_fig_rf:Footprint', 'non-finite footprint for an included cell (mask / ok_cell inconsistency)');
end


function r = footprint_radius(th, s_prox, s_dist, s_orth, cq)
% Half-maximum radius of the asymmetric 2-D Gaussian in direction th (PD frame;
% u = cos th along PD with proximal = negative u, v = sin th orthogonal).
    su = s_dist * ones(size(th)); su(cos(th) < 0) = s_prox;
    r = cq ./ sqrt(cos(th).^2 ./ su.^2 + sin(th).^2 ./ s_orth.^2);
end


%% ========================================================================
%  Plotting helpers
%% ========================================================================
function fig = new_fig(w, h)
    fig = figure('Units', 'centimeters', 'Position', [2 2 w h], 'PaperUnits', 'centimeters', ...
        'PaperSize', [w h], 'PaperPosition', [0 0 w h], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', 'Helvetica', 'DefaultTextFontName', 'Helvetica', 'DefaultAxesFontSize', 7);
end


function cm = redblue_local(cl)
% blue -> white -> red with white at 0 for the colour range cl = [lo hi], lo < 0 < hi.
    n = 256; nb = round(n * -cl(1) / (cl(2) - cl(1))); nr = n - nb;
    tb = linspace(0, 1, nb)'; tr = linspace(1, 0, nr)';
    cm = [[tb, tb, ones(nb, 1)]; [ones(nr, 1), tr, tr]];
end


function cl = pop_clim(Mx, masks)
% Colour range for group-mean maps: white at 0, limits rounded to 2 mV over the four groups.
    lo = 0; hi = 0;
    for g = 1:4, mm = mean(Mx(:, :, masks{g}), 3, 'omitnan'); lo = min(lo, min(mm(:))); hi = max(hi, max(mm(:))); end
    cl = [-2 * max(1, ceil(-lo / 2)), 2 * ceil(hi / 2)];
end


function mm = group_mean_map(stack)
% Group mean of aligned signed maps; pixels with data from fewer than half the cells blank.
    mm = mean(stack, 3, 'omitnan');
    mm(sum(isfinite(stack), 3) < ceil(size(stack, 3) / 2)) = NaN;
end


function draw_group_map(ax, mm, cl, Udeg, XL, YL, rm, th)
% Masked group-mean signed map MM in the PD frame (proximal left), with the mean
% half-maximum footprint RM, OD line, PD arrow (from the distal edge of the footprint
% through the centre to 9 deg beyond the proximal edge) and the RF centre.
    imagesc(ax, Udeg, Udeg, mm, 'AlphaData', isfinite(mm)); set(ax, 'YDir', 'normal');
    colormap(ax, redblue_local(cl)); clim(ax, cl);
    plot(ax, rm .* cos(th), rm .* sin(th), 'k-', 'LineWidth', 0.8);
    r_prox = interp1(th, rm, pi); r_orth = interp1(th, rm, pi / 2); r_dist = interp1(th, rm, 0);
    plot(ax, [0 0], r_orth * [-1 1], 'k-', 'LineWidth', 1.0);
    text(ax, 0.7, r_orth - 1.5, 'OD', 'FontSize', 6.5, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    quiver(ax, r_dist, 0, -(r_dist + r_prox + 9), 0, 0, 'Color', 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.22);
    text(ax, -r_prox - 4.5, -1.8, 'PD', 'FontSize', 6.5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot(ax, 0, 0, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 0.8);
    text(ax, XL(1) + 0.5, YL(2) - 0.5, 'proximal', 'FontSize', 5.5, 'VerticalAlignment', 'top');
    text(ax, XL(2) - 0.5, YL(2) - 0.5, 'distal',   'FontSize', 5.5, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    axis(ax, 'image'); xlim(ax, XL); ylim(ax, YL); set(ax, 'TickDir', 'out', 'XTick', -15:5:10, 'YTick', -10:5:10);
end


function draw_footprint_panel(ax, t, FP, masks, th, cols, gname, show_ylabel)
% Half-maximum footprints of cell type t (1 = T4, 2 = T5), ctrl vs tutl-, centred, PD frame
% (proximal left): line = pointwise mean of the per-cell contours FP{g} (cells x angles),
% band = pointwise +- SEM.
    hold(ax, 'on');
    for gi = [2*t-1, 2*t]
        R = FP{gi};
        rm = mean(R, 1); se = std(R, 0, 1) / sqrt(size(R, 1));
        ri = rm - se; ro = rm + se;
        fill(ax, [ri .* cos(th), fliplr(ro .* cos(th))], [ri .* sin(th), fliplr(ro .* sin(th))], cols{gi}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot(ax, rm .* cos(th), rm .* sin(th), '-', 'Color', cols{gi}, 'LineWidth', 1.6);
    end
    plot(ax, 0, 0, 'k+', 'MarkerSize', 6); xline(ax, 0, ':', 'Color', [0.6 0.6 0.6]); yline(ax, 0, ':', 'Color', [0.6 0.6 0.6]);
    axis(ax, 'equal'); xlim(ax, [-13 13]); ylim(ax, [-13 13]); set(ax, 'TickDir', 'out', 'XTick', -10:5:10, 'YTick', -10:5:10); box(ax, 'off');
    title(ax, sprintf('%s: half-maximum footprint', gname{2*t-1}(1:2)), 'FontWeight', 'normal', 'FontSize', 6.5);
    xlabel(ax, 'PD-ND axis (deg)'); if show_ylabel, ylabel(ax, 'orthogonal (deg)'); end
    text(ax, -12.5, 12.5, sprintf('control n = %d', sum(masks{2*t-1})), 'Color', cols{2*t-1}, 'FontSize', 6, 'VerticalAlignment', 'top');
    text(ax, -12.5, 11,   sprintf('tutl- n = %d',   sum(masks{2*t})),   'Color', cols{2*t},   'FontSize', 6, 'VerticalAlignment', 'top');
    text(ax, -12.5, -12.5, 'line: group mean; band: +- SEM', 'FontSize', 5.5, 'VerticalAlignment', 'bottom');
end


function draw_box_panel(ax, vals, masks, p_ct, y_label, y_limits, y_ticks, cols, rs)
% Box + dot panel for the four groups (T4 ctrl, T4 tutl-, T5 ctrl, T5 tutl-), following
% draw_boxplot_panel in generate_manuscript_fig_ds.m: colour-matched boxchart (median, IQR,
% whiskers), jittered dots, floating axes, two-row x labels, and the ctrl-vs-tutl- bracket
% per cell type with the rank-sum p in P_CT = [p_T4 p_T5] (computed once with the printout).
    colors = vertcat(cols{:});
    box_colors = min(1, colors * 0.3 + 0.7);
    hold(ax, 'on');
    for g = 1:4
        v = vals(masks{g});
        boxchart(ax, g * ones(size(v)), v, 'BoxFaceColor', box_colors(g, :), 'BoxEdgeColor', colors(g, :), ...
            'WhiskerLineColor', colors(g, :), 'MarkerStyle', 'none', 'BoxWidth', 0.5, 'LineWidth', 1.0);
        x = g + 0.25 * (rand(rs, size(v)) - 0.5);
        scatter(ax, x, v, 20, colors(g, :), 'filled', 'MarkerFaceAlpha', 0.5);
    end
    ylim(ax, y_limits); xlim(ax, [0.5 4.5]); ax.Clipping = 'off';
    set(ax, 'YTick', y_ticks, 'TickDir', 'out', 'FontSize', 7, 'XTick', 1:4, 'XTickLabel', []); box(ax, 'off');
    ax.XColor = 'none'; ax.YColor = 'none';
    line(ax, [1 4], y_limits(1) * [1 1], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
    tick_arm_y = 0.015 * diff(y_limits);
    for g = 1:4, line(ax, [g g], [y_limits(1), y_limits(1) - tick_arm_y], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off'); end
    line(ax, [0.5 0.5], [y_ticks(1) y_ticks(end)], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
    tick_arm_x = 0.12;
    for ti = 1:numel(y_ticks)
        line(ax, [0.5, 0.5 - tick_arm_x], y_ticks(ti) * [1 1], 'Color', 'k', 'LineWidth', 0.4, 'Clipping', 'off');
        text(ax, 0.5 - tick_arm_x - 0.05, y_ticks(ti), num2str(y_ticks(ti)), 'FontSize', 7, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
    ylabel(ax, y_label, 'FontSize', 8, 'Color', 'k');
    tx = {'ctrl', '{\ittutl-}', 'ctrl', '{\ittutl-}'};
    for g = 1:4
        text(ax, g, y_limits(1) - 0.06 * diff(y_limits), tx{g}, 'FontSize', 6, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', 'Interpreter', 'tex', 'Color', colors(g, :));
    end
    text(ax, 1.5, y_limits(1) - 0.16 * diff(y_limits), 'T4', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(ax, 3.5, y_limits(1) - 0.16 * diff(y_limits), 'T5', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    y0 = y_limits(2) - 0.02 * diff(y_limits);
    draw_bracket(ax, 1, 2, y0, p_ct(1));
    draw_bracket(ax, 3, 4, y0, p_ct(2));
end


function draw_bracket(ax, x1, x2, y, p)
% As in generate_manuscript_fig_ds.m: bracket with asterisks.
    yl = ylim(ax); arm = 0.02 * diff(yl);
    line(ax, [x1 x1 x2 x2], [y - arm, y, y, y - arm], 'Color', 'k', 'LineWidth', 0.25);
    if p < 0.001, s = '***'; elseif p < 0.01, s = '**'; elseif p < 0.05, s = '*'; else, s = 'n.s.'; end
    text(ax, mean([x1 x2]), y, s, 'FontSize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
