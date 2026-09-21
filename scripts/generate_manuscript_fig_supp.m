% GENERATE_MANUSCRIPT_FIG_SUPP  Build the supplementary manuscript figure.
%   Thin wrapper around generate_manuscript_fig('supp').
%
%   Edit DATA_ROOT below to point at your local copy of the 1DRF dataset.
%   See MANUSCRIPT_FIGURES.md §2 for the expected directory layout.

DATA_ROOT = '/Users/reiserm/Documents/ttl_1DRF';   % <-- edit for your setup

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));
% Manuscript settings (also the defaults since 2026-09-21): per-cell position test, no FWHM bars.
generate_manuscript_fig('supp', struct('data_root', DATA_ROOT, 'pool_test', 'percell', 'show_fwhm', false));
