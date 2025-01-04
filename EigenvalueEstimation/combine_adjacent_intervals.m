function combined_intervals = combine_adjacent_intervals(intervals)
    % Vectorized function to merge adjacent or overlapping intervals
    % Input: intervals - an array of infsup-type intervals
    % Output: combined_intervals - an array of merged intervals

    if isempty(intervals)
        combined_intervals = intervals;
        return;
    end

    % Extract the lower and upper bounds of the intervals
    lower_bounds = I_inf(intervals);
    upper_bounds = I_sup(intervals);

    % Initialize merged bounds
    merged_lower = lower_bounds(1);
    merged_upper = upper_bounds(1);

    % Store results in arrays
    combined_lower = [];
    combined_upper = [];

    for k = 2:length(intervals)
        if ~(merged_upper < lower_bounds(k))
            % Merge overlapping or adjacent intervals
            merged_upper = max(merged_upper, upper_bounds(k));
        else
            % Save the current merged interval
            combined_lower = [combined_lower, merged_lower];
            combined_upper = [combined_upper, merged_upper];
            % Start a new merged interval
            merged_lower = lower_bounds(k);
            merged_upper = upper_bounds(k);
        end
    end

    % Add the last merged interval
    combined_lower = [combined_lower, merged_lower];
    combined_upper = [combined_upper, merged_upper];

    % Create the output as an array of infsup intervals
    combined_intervals = I_infsup(combined_lower, combined_upper)';
end
