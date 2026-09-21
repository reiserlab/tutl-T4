function fig_combined = generate_manuscript_fig(mode, opts)
% GENERATE_MANUSCRIPT_FIG  Combined manuscript figure (main or supp).
%
%   FIG_COMBINED = GENERATE_MANUSCRIPT_FIG(MODE) builds a single 18 x 16 cm
%   figure for the T4/T5 *tutl-* manuscript by calling the flash sub-figure
%   (top, ~42%) and the bar-sweep sub-figure (bottom, ~53%) and copying
%   their axes + annotations into a unified canvas.
%
%       MODE = 'main'  -> PD-axis flash (top) + 56 dps bar sweep (bottom)
%       MODE = 'supp'  -> orthogonal flash (top) + 28 dps bar sweep (bottom)
%
%   FIG_COMBINED = GENERATE_MANUSCRIPT_FIG(MODE, OPTS) accepts an options
%   struct with fields:
%       .data_root   - data_root (default '/Users/reiserm/Documents/ttl_1DRF')
%       .skip_export - true to suppress PDF/PNG export (default false)
%       .show_fwhm   - FWHM bars on the flash depolarization panels (default true)
%       .pool_test   - 'pooled' (default) or 'percell'; see GENERATE_MANUSCRIPT_FIG_EF
%
%   See also GENERATE_MANUSCRIPT_FIG_EF, GENERATE_MANUSCRIPT_FIG_DS.

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'data_root'),   opts.data_root   = '/Users/reiserm/Documents/ttl_1DRF'; end
if ~isfield(opts, 'skip_export'), opts.skip_export = false; end
if ~isfield(opts, 'show_fwhm'),   opts.show_fwhm   = true; end
if ~isfield(opts, 'pool_test'),   opts.pool_test   = 'pooled'; end
if ~isfield(opts, 'fwhm_method'), opts.fwhm_method = 'interp'; end

switch lower(mode)
    case 'main'
        axis_mode = 'pd';
        speed_dps = 56;
        out_tag   = 'fig_main_ABCDEF_v3';
    case 'supp'
        axis_mode = 'ortho';
        speed_dps = 28;
        out_tag   = 'fig_supp_ABCDEF_v3';
    otherwise
        error('generate_manuscript_fig:BadMode', ...
            'mode must be ''main'' or ''supp'', got %s', mode);
end

data_root = opts.data_root;
out_dir   = fullfile(data_root, 'manuscript_figures');
if ~isfolder(out_dir) && ~opts.skip_export, mkdir(out_dir); end

%% ===================== Generate sub-figures =============================
sub_opts = struct('data_root', data_root, 'skip_export', true, ...
    'show_fwhm', opts.show_fwhm, 'pool_test', opts.pool_test, 'fwhm_method', opts.fwhm_method);
variant_tag = '';
if ~opts.show_fwhm, variant_tag = [variant_tag '_noFWHM']; end
if strcmpi(opts.pool_test, 'percell'), variant_tag = [variant_tag '_percell']; end
if strcmpi(opts.fwhm_method, 'gauss'), variant_tag = [variant_tag '_gaussFWHM']; end
out_tag = [out_tag variant_tag];
fprintf('\n========== Generating flash sub-figure (axis_mode=%s) ==========\n', axis_mode);
fig_ef   = generate_manuscript_fig_ef(axis_mode, sub_opts);
fprintf('\n========== Generating bar-sweep sub-figure (speed=%d) ==========\n', speed_dps);
fig_efgh = generate_manuscript_fig_ds(speed_dps, sub_opts);

%% ===================== Composite canvas =================================
fprintf('\n========== Combining into %s figure ==========\n', mode);
FIG_W = 18;  FIG_H = 16;
fig_combined = figure('Units', 'centimeters', 'Position', [2 2 FIG_W FIG_H], ...
    'PaperUnits', 'centimeters', 'PaperSize', [FIG_W FIG_H], ...
    'PaperPosition', [0 0 FIG_W FIG_H], 'Color', 'w');
set(fig_combined, 'DefaultAxesFontName', 'Helvetica', 'DefaultTextFontName', 'Helvetica');

dummy_ann = annotation(fig_combined, 'textbox', [0 0 0.01 0.01], 'String', '');
ann_layer = dummy_ann.Parent;
delete(dummy_ann);

% --- Top portion: flash (CD or AB) ---
CD_S = 0.42;  CD_O = 0.56;
copy_axes_and_annotations(fig_ef,   fig_combined, ann_layer, CD_S,   CD_O);

% --- Bottom portion: bar sweep (EFGH or CDEF) ---
EFGH_S = 0.53; EFGH_O = 0.01;
copy_axes_and_annotations(fig_efgh, fig_combined, ann_layer, EFGH_S, EFGH_O);

close(fig_ef);
close(fig_efgh);

%% ===================== Export ===========================================
ts = datestr(now, 'yyyymmdd_HHMM');
pdf_file = fullfile(out_dir, sprintf('%s_%s.pdf', out_tag, ts));
png_file = fullfile(out_dir, sprintf('%s_%s.png', out_tag, ts));

if ~opts.skip_export
    exportgraphics(fig_combined, pdf_file, 'ContentType', 'vector');
    exportgraphics(fig_combined, png_file, 'Resolution', 300);
    fprintf('Saved: %s\n', pdf_file);
    fprintf('Saved: %s\n', png_file);
end
fprintf('=== generate_manuscript_fig(%s) done ===\n', mode);

end


function copy_axes_and_annotations(src_fig, dst_fig, ann_layer, y_scale, y_offset)
% COPY_AXES_AND_ANNOTATIONS  Copy axes + annotations from src to dst with y-rescale.
%   y_scale, y_offset: map src y in [0,1] -> [y_offset, y_offset + y_scale].

    src_axes = [findobj(src_fig, 'Type', 'axes'); findobj(src_fig, 'Type', 'polaraxes')];
    for i = 1:numel(src_axes)
        ax_new = copyobj(src_axes(i), dst_fig);
        pos = ax_new.Position;
        ax_new.Position = [pos(1), pos(2)*y_scale + y_offset, pos(3), pos(4)*y_scale];
    end

    src_anns = findall(src_fig, '-isa', 'matlab.graphics.shape.TextBox');
    src_anns = [src_anns; findall(src_fig, '-isa', 'matlab.graphics.shape.Arrow')];
    src_anns = [src_anns; findall(src_fig, '-isa', 'matlab.graphics.shape.Line')];
    for i = 1:numel(src_anns)
        ann_new = copyobj(src_anns(i), ann_layer);
        if isprop(ann_new, 'Position')
            pos = ann_new.Position;
            if numel(pos) >= 4
                ann_new.Position = [pos(1), pos(2)*y_scale + y_offset, pos(3), pos(4)*y_scale];
            elseif numel(pos) >= 2
                ann_new.Position = [pos(1), pos(2)*y_scale + y_offset];
            end
        end
    end

    fprintf('  Copied %d axes + %d annotations.\n', numel(src_axes), numel(src_anns));
end
