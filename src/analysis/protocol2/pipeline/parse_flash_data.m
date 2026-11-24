function [data_comb, cmap_id, var_across_reps, var_within_reps, diff_mean, max_data, min_data, mean_traces] = parse_flash_data(f_data, v_data, on_off, slow_fast, px_size, PROJECT_ROOT, filter_idx)
% PARSE_FLASH_DATA  Extract and analyze flash stimulus responses from recording.
%
%   [DATA_COMB, CMAP_ID, VAR_ACROSS, VAR_WITHIN, DIFF_MEAN, MAX_DATA, MIN_DATA]
%       = PARSE_FLASH_DATA(F_DATA, V_DATA, ON_OFF, SLOW_FAST, PX_SIZE, PROJECT_ROOT)
%       = PARSE_FLASH_DATA(..., FILTER_IDX)
%   identifies flash stimulus timing, extracts responses, and computes
%   spatial response maps for receptive field analysis.
%
%   INPUTS:
%     f_data       - 1xN array of frame numbers at each time point
%     v_data       - 1xN array of voltage values at each time point (10kHz)
%     on_off       - 'on' or 'off' indicating flash polarity
%     slow_fast    - 'slow' for standard timing (currently only 'slow' used)
%     px_size      - Flash size in pixels (4 or 6)
%     PROJECT_ROOT - Base path for saving quality check figures
%     filter_idx   - (optional) If true (default), apply Autumn idx filtering
%                    for OFF cells. Set false for Summer 2025 data.
%
%   OUTPUTS:
%     data_comb      - NxN matrix of response values at each grid position
%                      N=14 for 4px flashes, N=10 for 6px flashes
%     cmap_id        - NxN matrix of response classification:
%                      1 = excitatory (positive), 2 = inhibitory (negative),
%                      3 = neutral (below threshold)
%     var_across_reps - NxN matrix of coefficient of variation across 3 reps
%     var_within_reps - NxN matrix of CV within each rep (temporal variance)
%     diff_mean       - NxN matrix of max-min response difference
%     max_data        - NxN matrix of 98th percentile of mean response
%     min_data        - NxN matrix of 2nd percentile of mean response
%     mean_traces     - (optional) NxNxT array: the 3-rep mean, median-referenced
%                      voltage trace per grid position (window as extracted below)
%
%   RESPONSE CLASSIFICATION:
%     Excitatory: |max| >= |min| AND diff_resp > 3
%     Inhibitory: |max| < |min| AND diff_resp > 2.8
%     Neutral: response below threshold
%
%   TIMING:
%     Extracts a 0.7 s window per flash (100 ms before the flash frame, 160 ms flash, ~440 ms after)
%     3 repetitions averaged for final response estimate
%
%   GRID STRUCTURE:
%     4px flashes: 14x14 grid = 196 positions, 50% overlap
%     6px flashes: 10x10 grid = 100 positions, 50% overlap
%     Rows correspond to vertical position (y), columns to horizontal (x)
%
%   See also PROCESS_FLASH_P2, GAUSSIAN_RF_ESTIMATE, LOAD_PROTOCOL2_DATA
% ________________________________________________________________________

    if nargin < 7
        filter_idx = true;  % backward compatible: existing callers get Autumn behavior
    end

    diff_f_data = diff(f_data);
    median_v = median(v_data);
    v2_data = v_data - median_v;

    if px_size == 4
        n_flashes = 196;
        n_rows_cols = 14;
    elseif px_size == 6
        n_flashes = 100;
        n_rows_cols = 10;
    end 

    % Find the timepoints when the flash stimuli end. 
    if on_off == "on"
        % 1 - find the start of the 4px flashes.
        idx1 = find(diff_f_data == 197 & f_data(2:end) == 197); 
        idx1([2,4,6]) = [];
        % 2 - find the start of the 6px flashes.
        idx2 = find(diff_f_data == 101 & f_data(2:end) == 101); 
        idx = [idx1, idx2];
        idx = sort(idx);

        % For finding the end of the 6 pixel flashes
        drop_at_end = -200; % ON flashes. the last 6 pixel flash is frame 200.

    elseif on_off == "off"
        % Both 4px and 6px flashes start with frame 1. 
        % 'idx' is an array of the timepoints when the 4px flashes start
        % and the 6px flashes start. 'idx' should be an array of size [1,6]
        idx = find(diff_f_data == 1 & f_data(2:end) == 1); % First flash.
        if filter_idx
            idx = idx([1,2,5,6, 9, 10]);  % Autumn: filter to 6 rep-boundary transitions
        end

        % For finding the end of the 6 pixel flashes
        drop_at_end = -100; % ON flashes. the last 6 pixel flash is frame 200.
    end 


    % This finds when the difference in frame number == drop_at_end. 
    % This should find 6 values. Values 1, 3 and 5 are wthin the 4 pixel
    % flashes and values 2,4 and 6 correspond to the timing of the very end
    % of the last flash.
    idx_6px = find(diff_f_data == drop_at_end); % where the flash stimuli end.
    idx_6px([1,3,5]) = []; 

    % TEST - 'idx' values = start of the different groups of flash stimuli. 
    % figure; 
    % plot(f_data);
    % hold on;
    % for iii = 1:numel(idx)
    %     plot([idx(iii), idx(iii)], [0 400], 'm');
    % end 
    
    % % 340 ms OFF - bkg - 160 ms FLASH. 500 ms = 0.5s. 
    % % 170 ms OFF - bkg - 80 ms FLASH - 0.25s

    % % July 2025 onwards - increased interval time from 340ms to 440ms. 
    if slow_fast == "slow"
        slow_flashes_dur = 7000; % 0.6s * sampling rate.
        % flash_dur_ms_slow = 976700;
    % elseif slow_fast == "fast"
    %     fast_flashes_dur = 2500;
    %     flash_dur_ms_fast = 976700/2;
    end 
    
    % Empty arrays to fill - 196 flashes/values. 
    data_comb = zeros(n_rows_cols, n_rows_cols);
    cmap_id = zeros(n_rows_cols, n_rows_cols);
    var_across_reps = zeros(n_rows_cols, n_rows_cols);
    var_within_reps = zeros(n_rows_cols, n_rows_cols);
    diff_mean = zeros(n_rows_cols, n_rows_cols);
    max_data = zeros(n_rows_cols, n_rows_cols);
    min_data = zeros(n_rows_cols, n_rows_cols);
    % i_num = zeros(14,14);
    
    % Run through the responses to each flash.
    for i = 1:n_flashes

        if slow_fast == "slow"

            % Array to collect the data for each flash over the three repetitions.
            data_flash = ones(3, slow_flashes_dur); 
            data_flash_raw = ones(3, slow_flashes_dur); 
            data_frame = ones(3, slow_flashes_dur); 
        

            for r = 1:3
        
                if r == 1 % rep 1
                    if px_size == 4
                        rng_rep1 = (idx(1):idx(2));
                        start_idx = idx(1);
                    else 
                        rng_rep1 = (idx(2):idx_6px(1));
                        start_idx = idx(2);
                    end 
                    start_flash_idxs = find(diff(f_data(rng_rep1))>0)+start_idx-1;
                elseif r == 2 % rep 2 
                    if px_size == 4
                        rng_rep2 = idx(3):idx(4);
                        start_idx = idx(3);
                    else
                        rng_rep2 = (idx(4):idx_6px(2));
                        start_idx = idx(4);
                    end 
                    start_flash_idxs = find(diff(f_data(rng_rep2))>0)+start_idx-1;
                elseif r == 3 % rep3 
                    if px_size == 4
                        rng_rep3 = idx(5):idx(6);
                        start_idx = idx(5);
                    else
                        rng_rep3 = (idx(6):idx_6px(3));
                        start_idx = idx(6);
                    end
                    start_flash_idxs = find(diff(f_data(rng_rep3))>0)+start_idx-1;
                end 
        
                % Extract data 1000 timepoints before the flash starts
                % til the end of the interval before the next flash. 
                d = f_data(start_flash_idxs(i)-1000:start_flash_idxs(i)+6000-1); % frame data during flash. 
                v = v2_data(start_flash_idxs(i)-1000:start_flash_idxs(i)+6000-1); % median-subtracted voltage data
                v_raw = v_data(start_flash_idxs(i)-1000:start_flash_idxs(i)+6000-1);

                data_frame(r, :) = d;
                data_flash(r, :) = v;
                data_flash_raw(r, :) = v_raw;
            end

        end 
        
        %% Check - find beginning of each individual flash

        % figure; plot(f_data); hold on;
        % plot([start_flash_idxs(i), start_flash_idxs(i)], [0 400], 'r') % start of flash 1
        % plot([start_flash_idxs(2), start_flash_idxs(2)], [0 400], 'r') % start of flash 2

        % Show the actual range extracted
        % figure; plot(f_data); hold on;
        % plot([start_flash_idxs(i)-1000, start_flash_idxs(i)-1000], [0 400], 'c')
        % plot([start_flash_idxs(i)+6000-1, start_flash_idxs(i)+6000-1], [0 400], 'c')


        %%
        % Check for the variance across reps:
        % var_X_rep = var(var(data_flash));
    
        % Standard deviation within each repetition (across time points)
        std_within_trial = std(data_flash_raw, 0, 2);  % std along rows (across time points)
        
        % Mean within each repetition (across time points)
        mean_within_trial = mean(data_flash_raw, 2);  % mean along rows (across time points)
        
        % Coefficient of Variation (CV) within each repetition
        cv_within_trial = std_within_trial ./ mean_within_trial;
        
        % Mean CV across all repetitions (mean within trial CV)
        mean_cv_within_trial = abs(mean(cv_within_trial));
    
        % Standard deviation across time points (across trials)
        std_across_trials = std(data_flash_raw, 0, 1);  % std along columns (across trials)
        
        % Mean across time points (across trials)
        mean_across_trials = mean(data_flash_raw, 1);  % mean along columns (across trials)
        
        % Coefficient of Variation (CV) across repetitions (across trials)
        cv_across_trials = std_across_trials ./ mean_across_trials;
    
        % max_cv_across_trials = prctile(cv_across_trials, 98);
        mean_cv_across_trials = abs(mean(cv_across_trials));
    
        mean_data_flash = mean(data_flash);
        n_vals = size(mean_data_flash, 2);
        if nargout >= 8 && i == 1, mean_traces = NaN(n_rows_cols, n_rows_cols, n_vals); end   % optional 8th output: mean trace per square
    
        % % % % % % % TODO - - -- update this based on the 1000
        % datapoints before the start of the flash now. 

        % Max and min of the mean flash response:
        if slow_fast == "slow"
            max_val_flash = prctile(mean_data_flash(500:end), 98); % Max 
            % min_val_flash = prctile(mean_data_flash(1250:end), 2); % Min in the second half - ignore if min is early. 
            min_val_flash = prctile(mean_data_flash(2500:end), 2); % Min in the second half - ignore if min is early. 
        elseif slow_fast == "fast"
            max_val_flash = prctile(mean_data_flash(200:end), 98); % Max 
            % min_val_flash = prctile(mean_data_flash(1250:end), 2); % Min in the second half - ignore if min is early. 
            min_val_flash = prctile(mean_data_flash(1250:end), 2); % Min in the second half - ignore if min is early. 
        end 
    
        diff_resp = max_val_flash - min_val_flash; 
        diff_med = median_v - min_val_flash;
    
        if abs(max_val_flash)>=abs(min_val_flash) % larger excitatory peak
            
            if diff_resp>3
                    val = max_val_flash;
                    cm = 1;
            else 
                val = mean(mean_data_flash(n_vals*0.75:end));
                cm = 3;
            end 
        elseif abs(max_val_flash)<abs(min_val_flash) % larger inhibitory peak 
            if diff_resp>2.8
                val = min_val_flash;
                cm = 2;
            else 
                val = mean(mean_data_flash(n_vals*0.75:end));
                cm = 3;
            end 
        end 
    
        flash_frame_num = max(d)-1;
        % fnum(i) = flash_frame_num;
        % disp(flash_frame_num)
    
        if on_off == "on" % from 196
            rows = n_rows_cols - mod(flash_frame_num - n_flashes, n_rows_cols);   % Rows decrease from 14 to 1
            cols = floor((flash_frame_num - n_flashes) / n_rows_cols) + 1; % Columns increase normally
        elseif on_off == "off" % 1- 196
            rows = n_rows_cols - mod(flash_frame_num, n_rows_cols);   % Rows decrease from 14 to 1
            cols = floor(flash_frame_num / n_rows_cols) + 1; % Columns increase normally
        end

        data_comb(rows, cols) = val;
        if nargout >= 8, mean_traces(rows, cols, :) = mean_data_flash; end
        cmap_id(rows, cols) = cm;
        var_across_reps(rows, cols) = mean_cv_across_trials;
        var_within_reps(rows, cols) = mean_cv_within_trial;
        diff_mean(rows, cols) = diff_resp;
        max_data(rows, cols) = max_val_flash;
        min_data(rows, cols) = min_val_flash;
        % i_num(rows, cols) = i;
    end 

    % Quality check figures:
    % save_fig = 0;
    % plot_quality_check_var_reps_prctile(var_across_reps, var_within_reps, diff_mean, max_data, min_data, save_fig, PROJECT_ROOT)
 
end 