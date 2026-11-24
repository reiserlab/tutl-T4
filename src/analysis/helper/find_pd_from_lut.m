function pd_info = find_pd_from_lut(max_v, lut_directions, lut_orientations, ...
    lut_patterns, lut_functions, plot_order, Tbl, pattern_offset)
% FIND_PD_FROM_LUT  Determine preferred direction and map to bar flash columns.
%
%   PD_INFO = FIND_PD_FROM_LUT(MAX_V, LUT_DIRECTIONS, LUT_ORIENTATIONS, ...
%       LUT_PATTERNS, LUT_FUNCTIONS, PLOT_ORDER, TBL, PATTERN_OFFSET)
%   computes the vector sum of directional bar sweep responses, finds the
%   closest LUT direction, and maps the preferred direction (PD) to the
%   corresponding bar flash orientation column. Also identifies the
%   orthogonal orientation and determines spatial ordering (ND to PD).
%
%   INPUTS:
%     max_v            - 16x1 array of depolarization amplitudes, ordered
%                        by subplot position (from compute_bar_sweep_responses)
%     lut_directions   - 16x1 array of LUT motion directions in degrees
%     lut_orientations - 16x1 array of LUT bar orientations in degrees
%     lut_patterns     - 16x1 array of experiment pattern numbers
%     lut_functions    - 16x1 array of experiment function numbers
%     plot_order       - 1x16 subplot-to-data-row mapping
%     Tbl              - Full LUT table (for orthogonal orientation lookup)
%     pattern_offset   - Offset from experiment pattern to full-field pattern
%                        index (default: 2, so experiment pattern 3 maps to
%                        full-field pattern 1)
%
%   OUTPUT:
%     pd_info - Structure with fields:
%       .resultant_angle   - Vector sum angle in radians
%       .pd_direction      - Closest LUT direction in degrees
%       .pd_orientation    - Bar orientation at PD in degrees
%       .pd_pattern        - Experiment pattern number at PD
%       .pd_function       - Experiment function number at PD
%       .pd_data_row       - Data row index for PD
%       .bar_flash_col     - Bar flash column index for PD orientation (1-8)
%       .ortho_flash_col   - Bar flash column for orthogonal orientation (1-8)
%       .ortho_orientation - Orthogonal bar orientation in degrees
%       .pos_order         - 1x11 spatial ordering array (ND side to PD side)
%       .ortho_pos_order   - 1x11 spatial ordering for orthogonal axis
%                            (right-hand rule: 90° CW from PD defines +)
%
%   VECTOR SUM METHOD:
%     The preferred direction is computed as the angle of the weighted
%     vector sum across all 16 directions, using depolarization amplitudes
%     as weights. The closest discrete LUT direction is then selected.
%
%   SPATIAL ORDERING:
%     Forward functions (odd: 3, 5, 7) sweep positions 1 to 11 in the
%     ND-to-PD direction. Reverse functions (even: 4, 6, 8) sweep in the
%     opposite direction, so positions are flipped to maintain ND-to-PD
%     ordering in the output.
%
%   See also VECTOR_SUM_POLAR, VERIFY_LUT_DIRECTIONS, COMPUTE_BAR_SWEEP_RESPONSES

    if nargin < 8, pattern_offset = 2; end

    % LUT directions ordered by subplot position
    lut_dirs_ordered_rad = deg2rad(lut_directions(plot_order));

    % Compute vector sum
    rho_for_vecsum   = max_v';                  % 1x16
    theta_for_vecsum = lut_dirs_ordered_rad';   % 1x16
    [~, resultant_angle] = vector_sum_polar(rho_for_vecsum, theta_for_vecsum);

    % Find closest LUT direction to vector sum angle
    angle_diffs = abs(lut_dirs_ordered_rad - resultant_angle);
    angle_diffs = min(angle_diffs, 2*pi - angle_diffs);  % handle wraparound
    [~, pd_subplot_idx] = min(angle_diffs);

    pd_data_row   = plot_order(pd_subplot_idx);
    pd_direction  = lut_directions(pd_data_row);
    pd_direction_vs = mod(rad2deg(resultant_angle), 360);   % continuous vector-sum angle, before snapping
    pd_orientation = lut_orientations(pd_data_row);
    pd_pattern    = lut_patterns(pd_data_row);
    pd_function   = lut_functions(pd_data_row);

    % Map PD to bar flash column
    pd_ff_pattern = pd_pattern - pattern_offset;
    bar_flash_col = pd_ff_pattern;

    % Orthogonal orientation: 4 columns away (wraps around 8 orientations)
    ortho_flash_col = mod(bar_flash_col - 1 + 4, 8) + 1;
    ortho_exp_pat   = ortho_flash_col + pattern_offset;

    % Locate one orthogonal LUT row (prefer function==3, fall back to any odd function).
    % We derive both ortho_orientation and ortho_dir_fwd from the SAME row so the two
    % stay consistent in the rare fallback case.
    ortho_fwd_mask = Tbl.pattern == ortho_exp_pat & Tbl.function == 3;
    if any(ortho_fwd_mask)
        ortho_lut_idx = find(ortho_fwd_mask, 1);
    else
        ortho_odd_mask = Tbl.pattern == ortho_exp_pat & mod(Tbl.function, 2) == 1;
        if ~any(ortho_odd_mask)
            error('find_pd_from_lut:NoOrthoLutRow', ...
                'No LUT row found for orthogonal pattern %d (any odd function).', ortho_exp_pat);
        end
        ortho_lut_idx = find(ortho_odd_mask, 1);
    end
    ortho_orientation = Tbl.orientation(ortho_lut_idx);
    ortho_dir_fwd     = Tbl.direction(ortho_lut_idx);

    % Determine PD spatial ordering (ND to PD, left to right)
    is_forward = mod(pd_function, 2) == 1;
    if is_forward
        pos_order = 1:11;
    else
        pos_order = 11:-1:1;
    end

    % Determine ortho spatial ordering via right-hand rule.
    % Right-hand rule: 90 deg clockwise from PD defines the positive ortho direction.
    rh_target = mod(pd_direction - 90, 360);
    ang_diff  = abs(ortho_dir_fwd - rh_target);
    ang_diff  = min(ang_diff, 360 - ang_diff);  % handle wraparound
    if ang_diff < 90
        ortho_pos_order = 1:11;
    else
        ortho_pos_order = 11:-1:1;
    end

    % Print summary
    fprintf('\n=== Preferred Direction (PD via vector sum) ===\n');
    fprintf('Vector sum angle: %.1f°\n', rad2deg(resultant_angle));
    fprintf('Closest LUT direction: %.1f°\n', pd_direction);
    fprintf('Bar orientation:       %.1f°\n', pd_orientation);
    fprintf('Pattern:               %d\n', pd_pattern);
    fprintf('Function:              %d\n', pd_function);
    fprintf('Data row:              %d\n', pd_data_row);
    fprintf('Full-field pattern index: %d\n', pd_ff_pattern);
    fprintf('Bar flash column:         %d\n', bar_flash_col);
    fprintf('Orthogonal bar flash column:  %d (orientation: %.1f°)\n', ...
        ortho_flash_col, ortho_orientation);
    if is_forward
        dir_label = 'forward';
    else
        dir_label = 'reverse';
    end
    fprintf('PD function %d is %s → position order: %s\n', ...
        pd_function, dir_label, mat2str(pos_order));
    fprintf('Ortho pos_order: %s (right-hand-rule target: %.1f°, ortho fwd: %.1f°)\n', ...
        mat2str(ortho_pos_order), rh_target, ortho_dir_fwd);

    % Pack into output struct
    pd_info.resultant_angle   = resultant_angle;
    pd_info.pd_direction      = pd_direction;
    pd_info.pd_direction_vs = pd_direction_vs;
    pd_info.pd_orientation    = pd_orientation;
    pd_info.pd_pattern        = pd_pattern;
    pd_info.pd_function       = pd_function;
    pd_info.pd_data_row       = pd_data_row;
    pd_info.bar_flash_col     = bar_flash_col;
    pd_info.ortho_flash_col   = ortho_flash_col;
    pd_info.ortho_orientation = ortho_orientation;
    pd_info.pos_order         = pos_order;
    pd_info.ortho_pos_order   = ortho_pos_order;

end
