function[flag_Vout_of_bounts] = vol_limits(Bus_matrix_v)
    flag_Vout_of_bounts = 0; 

    % Check V limits from Bus_matrix
    for i = 1 : 1 : length(Bus_matrix_v)      
        % Vmax = helping1_Bus_matrix(i, 8)
        % Vmin = helping1_Bus_matrix(i, 9)
        if(Bus_matrix_v(i, 6) <= Bus_matrix_v(i, 8) && Bus_matrix_v(i, 6) >= Bus_matrix_v(i, 9)) 
            % V(i) is in limits    
            flag_Vout_of_bounts = flag_Vout_of_bounts + 1;           
        else         
            disp('V out of bounds')     
        end
    end 
end