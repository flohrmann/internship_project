function [p, observeddifference, effectsize] = permutationTestWelch(sample1, sample2, varargin)
%https://de.mathworks.com/matlabcentral/fileexchange/63276-permutation-test
% parsing input
% permutationTest(rand(1,100), rand(1,100)-.25, 10000, ...
%          'plotresult', 1,
permutations = 10000;
p = inputParser;
addRequired(p, 'sample1', @isnumeric);
addRequired(p, 'sample2', @isnumeric);
addRequired(p, 'permutations', @isnumeric);
addParamValue(p, 'sidedness', 'both', @(x) any(validatestring(x,{'both', 'smaller', 'larger'})));
addParamValue(p, 'exact' , 0, @isnumeric);
addParamValue(p, 'plotresult', 0, @isnumeric);
addParamValue(p, 'showprogress', 0, @isnumeric);
parse(p, sample1, sample2, permutations, varargin{:})
sample1 = p.Results.sample1;
sample2 = p.Results.sample2;
permutations = p.Results.permutations;
sidedness = p.Results.sidedness;
exact = p.Results.exact;
plotresult = p.Results.plotresult;
showprogress = p.Results.showprogress;
% enforcing row vectors
if iscolumn(sample1), sample1 = sample1'; end
if iscolumn(sample2), sample2 = sample2'; end
allobservations = [sample1, sample2];
observeddifference = nanmean(sample1) - nanmean(sample2);
pooledstd = sqrt(  ( (numel(sample1)-1)*std(sample1)^2 + (numel(sample2)-1)*std(sample2)^2 )  /  ( numel(allobservations)-2 )  );
effectsize = observeddifference / pooledstd;
w = warning('off', 'MATLAB:nchoosek:LargeCoefficient');
if ~exact && permutations > nchoosek(numel(allobservations), numel(sample1))
    warning(['the number of permutations (%d) is higher than the number of possible combinations (%d);\n' ...
             'consider running an exact test using the ''exact'' argument'], ...
             permutations, nchoosek(numel(allobservations), numel(sample1)));
end
warning(w);
if showprogress, w = waitbar(0, 'Preparing test...', 'Name', 'permutationTest'); end
if exact
    % getting all possible combinations
    allcombinations = nchoosek(1:numel(allobservations), numel(sample1));
    permutations = size(allcombinations, 1);
end
% running test
randomdifferences = zeros(1, permutations);
if showprogress, waitbar(0, w, sprintf('Permutation 1 of %d', permutations), 'Name', 'permutationTest'); end
for n = 1:permutations
    if showprogress && mod(n,showprogress) == 0, waitbar(n/permutations, w, sprintf('Permutation %d of %d', n, permutations)); end
    
    % selecting either next combination, or random permutation
    if exact, permutation = [allcombinations(n,:), setdiff(1:numel(allobservations), allcombinations(n,:))];
    else, permutation = randperm(length(allobservations)); end
    
    % dividing into two samples
    randomSample1 = allobservations(permutation(1:length(sample1)));
    randomSample2 = allobservations(permutation(length(sample1)+1:length(permutation)));
    
    % saving differences between the two samples
    randomdifferences(n) = nanmean(randomSample1) - nanmean(randomSample2);
end
if showprogress, delete(w); end
% getting probability of finding observed difference from random permutations
if strcmp(sidedness, 'both')
    p = (length(find(abs(randomdifferences) > abs(observeddifference)))+1) / (permutations+1);
elseif strcmp(sidedness, 'smaller')
    p = (length(find(randomdifferences < observeddifference))+1) / (permutations+1);
elseif strcmp(sidedness, 'larger')
    p = (length(find(randomdifferences > observeddifference))+1) / (permutations+1);
end
% plotting result
if plotresult
    figure;
    if verLessThan('matlab', '8.4')
        % MATLAB R2014a and earlier
        hist(randomdifferences, 20);
    else
        % MATLAB R2014b and later
        histogram(randomdifferences, 20);
    end
    hold on;
    xlabel('Random differences');
    ylabel('Count')
    od = plot(observeddifference, 0, '*r', 'DisplayName', sprintf('Observed difference.\nEffect size: %.2f,\np = %f', effectsize, p));
    legend(od);
end
end

% function p_value = permutationTestWelch(data1, data2)
%     % Performs a permutation test using the Welch t-statistic.
%     %   data1 - Data for group 1 (ADHD group scores)
%     %   data2 - Data for group 2 (nonADHD group scores)
%     %   num_permutations - Number of permutations (1000?)
%     %   p_value - Permutation-based p-value
%     num_permutations = 1000;
%     % Compute observed Welch t-statistic
%     [~, ~, ~, stats] = ttest2(data1, data2, 'Vartype', 'unequal');
%     observed_t = stats.tstat;
% 
%     % Combine data
%     combined_data = [data1; data2];
%     n1 = length(data1);
%     n2 = length(data2);
%     n = n1 + n2;
% 
%     % Perform permutations
%     permuted_t_stats = zeros(num_permutations, 1);
%     for i = 1:num_permutations
%         % Shuffle indices
%         permuted_indices = randperm(n);
%         perm_group1 = combined_data(permuted_indices(1:n1));
%         perm_group2 = combined_data(permuted_indices(n1+1:end));
% 
%         % Compute Welch t-statistic for permuted data
%         [~, ~, ~, perm_stats] = ttest2(perm_group1, perm_group2, 'Vartype', 'unequal');
%         permuted_t_stats(i) = perm_stats.tstat;
%     end
% 
%     % Compute permutation p-value
%     p_value = mean(abs(permuted_t_stats) >= abs(observed_t));
% 
%     % Display results
%     fprintf('Observed Welch t-statistic: %.4f\n', observed_t);
%     fprintf('Permutation-based p-value: %.4f\n', p_value);
% end
