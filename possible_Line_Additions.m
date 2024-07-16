function [correct_line_for_LinesMatrix] = possible_Line_Additions()

    
    % Define the mapping of connections to their respective values
    connection_values = [1, 3, 4, 16, 18, 17];
    
    % List of connections
    connections = {[1 2], [2 3], [2 4], [9 10], [10 11], [9 14]};
       
    % Number of connections
    n = length(connections);
    
    % Initialize cell array to hold all combinations
    all_combinations = {};
    
    % Loop over all possible sizes of subsets (from 1 to n)
    for subset_size = 1:n
        % Get all combinations of indices for the current subset size
        index_combinations = nchoosek(1:n, subset_size);
        
        % Loop over each combination of indices
        for comb_idx = 1:size(index_combinations, 1)
            % Get the current combination of indices
            current_combination = index_combinations(comb_idx, :);
            
            % Use the current combination of indices to get the actual connections
            subset_indices = current_combination;
            
            % Map each connection index to its corresponding value
            mapped_subset = connection_values(subset_indices);
            
            % Add the current subset to the list of all combinations
            all_combinations{end+1} = mapped_subset;
        end
    end
    
    % Determine the maximum length of the subsets
    max_length = max(cellfun(@length, all_combinations));
    
    % Initialize a matrix to hold all combinations with padding
    output_matrix = NaN(length(all_combinations), max_length);
    
    % Populate the matrix, padding shorter rows with NaNs
    for i = 1:length(all_combinations)
        output_matrix(i, 1:length(all_combinations{i})) = all_combinations{i};
    end
    
    % Replace NaN values with zeros
    output_matrix(isnan(output_matrix)) = 0;
    
    % Transpose the matrix if needed
    output_matrix = output_matrix';
    
    % Return it 
    correct_line_for_LinesMatrix = output_matrix;
    
    %[row1,column1] = size(correct_line_for_LinesMatrix)
%     for k = 1:1:column1
%       
%        for j =1:1:row1
%            val = correct_line_for_LinesMatrix(j,k);
%            if(val==0)
%                break;
%            else              
%                %new_row = [Lines_matrix(val,:)];
%                matrix = [matrix; val];
%            end
%             
%        end
%     end
    
    
end
