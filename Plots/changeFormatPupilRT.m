function [avg_adhd_struct, avg_nonadhd_struct] = changeFormatPupilRT(avg_adhd, avg_nonadhd, num_before, num_after)

    % Initialize structs with empty arrays
    avg_adhd_struct = struct('a', [], 'a_simple', [], 'b', [], 'b_simple', []);
    avg_nonadhd_struct = struct('a', [], 'a_simple', [], 'b', [], 'b_simple', []);

    % Preallocate lists to accumulate data per group
    diam_around_stim_adhd = struct('a', [], 'a_simple', [], 'b', [], 'b_simple', []);
    diam_around_stim_nonadhd = struct('a', [], 'a_simple', [], 'b', [], 'b_simple', []);

    %% Process ADHD Group
    for i = 1:height(avg_adhd)
        % Extract participant data for each condition
        a_data = avg_adhd.a(i,:);
        a_simple_data = avg_adhd.a_simple(i,:);
        b_data = avg_adhd.b(i,:);
        b_simple_data = avg_adhd.b_simple(i,:);

        % Accumulate before/after data for ADHD group as rows
        diam_around_stim_adhd.a.before{i,:} = a_data(1:num_before);
        diam_around_stim_adhd.a.after{i,:} = a_data(num_before+1:num_before+num_after);

        diam_around_stim_adhd.a_simple.before{i,:} = a_simple_data(1:num_before);
        diam_around_stim_adhd.a_simple.after{i,:} = a_simple_data(num_before+1:num_before+num_after);

        diam_around_stim_adhd.b.before{i,:} = b_data(1:num_before);
        diam_around_stim_adhd.b.after{i,:} = b_data(num_before+1:num_before+num_after);

        diam_around_stim_adhd.b_simple.before{i,:} = b_simple_data(1:num_before);
        diam_around_stim_adhd.b_simple.after{i,:} = b_simple_data(num_before+1:num_before+num_after);
    end

    %% Process non-ADHD Group
    for i = 1:height(avg_nonadhd)
        % Extract participant data for each condition
        a_data = avg_nonadhd.a(i,:);
        a_simple_data = avg_nonadhd.a_simple(i,:);
        b_data = avg_nonadhd.b(i,:);
        b_simple_data = avg_nonadhd.b_simple(i,:);

        % Accumulate before/after data for non-ADHD group as rows
        diam_around_stim_nonadhd.a.before{i,:} = a_data(1:num_before);
        diam_around_stim_nonadhd.a.after{i,:} = a_data(num_before+1:num_before+num_after);

        diam_around_stim_nonadhd.a_simple.before{i,:} = a_simple_data(1:num_before);
        diam_around_stim_nonadhd.a_simple.after{i,:} = a_simple_data(num_before+1:num_before+num_after);

        diam_around_stim_nonadhd.b.before{i,:} = b_data(1:num_before);
        diam_around_stim_nonadhd.b.after{i,:} = b_data(num_before+1:num_before+num_after);

        diam_around_stim_nonadhd.b_simple.before{i,:} = b_simple_data(1:num_before);
        diam_around_stim_nonadhd.b_simple.after{i,:} = b_simple_data(num_before+1:num_before+num_after);
    end

    %% Write the results to the final structure with rows for each participant
    avg_adhd_struct.a.before = cell2mat(diam_around_stim_adhd.a.before);
    avg_adhd_struct.a.after = cell2mat(diam_around_stim_adhd.a.after);
    avg_adhd_struct.a_simple.before = cell2mat(diam_around_stim_adhd.a_simple.before);
    avg_adhd_struct.a_simple.after = cell2mat(diam_around_stim_adhd.a_simple.after);
    avg_adhd_struct.b.before = cell2mat(diam_around_stim_adhd.b.before);
    avg_adhd_struct.b.after = cell2mat(diam_around_stim_adhd.b.after);
    avg_adhd_struct.b_simple.before = cell2mat(diam_around_stim_adhd.b_simple.before);
    avg_adhd_struct.b_simple.after = cell2mat(diam_around_stim_adhd.b_simple.after);

    avg_nonadhd_struct.a.before = cell2mat(diam_around_stim_nonadhd.a.before);
    avg_nonadhd_struct.a.after = cell2mat(diam_around_stim_nonadhd.a.after);
    avg_nonadhd_struct.a_simple.before = cell2mat(diam_around_stim_nonadhd.a_simple.before);
    avg_nonadhd_struct.a_simple.after = cell2mat(diam_around_stim_nonadhd.a_simple.after);
    avg_nonadhd_struct.b.before = cell2mat(diam_around_stim_nonadhd.b.before);
    avg_nonadhd_struct.b.after = cell2mat(diam_around_stim_nonadhd.b.after);
    avg_nonadhd_struct.b_simple.before = cell2mat(diam_around_stim_nonadhd.b_simple.before);
    avg_nonadhd_struct.b_simple.after = cell2mat(diam_around_stim_nonadhd.b_simple.after);
end
