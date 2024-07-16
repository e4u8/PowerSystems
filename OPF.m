% Project 2 - Omada 8
% Optimal Power Flow - OPF.m

% Initialize Data of Power System

S_B = 100;          % Power Base = 100MVA
A = 0;              % cost parameter
lamda_Grid = 72;    % 72 €/MWh
matrix_with_accepted_data = [];         % Format: P_G1 | P_G2 | P_grid |L_for_Lagrange| F(P_G) | F_optimal με αγο�?οπωλησία |
volt_matrix_with_accepted_data = [];


% Bus data for step v
% bus | Type of bus | P_L (MW) | Q_L (MVAr) | B_s (MVAr) | Vm | Va | Vmax | Vmin
% Initializations for Load busses: V_m = 1 and δ_angles = Va = 0.
% Initializations for Voltage-controlled busses: V_m = Vg -> Vspec from Generators_matrix and δ_angles = Va = 0.

Bus_matrix_v = [ 1    3     0      0      0    1.06    0   1.1236   0.9964;              % Slack Bus
                 2    2    43.4    25.4   0    1.06    0   1.1236   0.9964;              % Voltage-controlled bus
                 3    2    188.4   38     0    1.03    0   1.0918   0.9682;              % Voltage-controlled bus
                 4    1    95.6    -7.8   0    1       0   1.06     0.94;
                 5    1    15.2    3.2    0    1       0   1.06     0.94;
                 6    2    22.4    15     0    1.07    0   1.1342   1.0058;              % Voltage-controlled bus
                 7    1    0       0      0    1       0   1.06     0.94;
                 8    2    0       0      0    1.06    0   1.1236   0.9964;              % Voltage-controlled bus
                 9    1    59      33.2   19   1       0   1.06     0.94;
                 10   1    18      11.6   0    1       0   1.06     0.94;
                 11   1    7       3.6    0    1       0   1.06     0.94;
                 12   1    12.2    3.2    0    1       0   1.06     0.94;
                 13   1    27      11.6   0    1       0   1.06     0.94;
                 14   1    29.8    10     0    1       0   1.06     0.94
               ];

% Bus data for step v+1
Bus_matrix_v_plus_1 = [ 1    3     0      0      0    1.06    0   1.1236   0.9964;              % Slack Bus
                        2    2    43.4    25.4   0    1.06    0   1.1236   0.9964;              % Voltage-controlled bus
                        3    2    188.4   38     0    1.03    0   1.0918   0.9682;              % Voltage-controlled bus
                        4    1    95.6    -7.8   0    1       0   1.06     0.94;
                        5    1    15.2    3.2    0    1       0   1.06     0.94;
                        6    2    22.4    15     0    1.07    0   1.1342   1.0058;              % Voltage-controlled bus
                        7    1    0       0      0    1       0   1.06     0.94;
                        8    2    0       0      0    1.06    0   1.1236   0.9964;              % Voltage-controlled bus
                        9    1    59      33.2   19   1       0   1.06     0.94;
                        10   1    18      11.6   0    1       0   1.06     0.94;
                        11   1    7       3.6    0    1       0   1.06     0.94;
                        12   1    12.2    3.2    0    1       0   1.06     0.94;
                        13   1    27      11.6   0    1       0   1.06     0.94;
                        14   1    29.8    10     0    1       0   1.06     0.94
                       ]; 

% Data of Transmission Lines
% from bus | to bus | R(p.u.) | X(p.u.)| B(p.u.)| transf.limit(p.u.)| δ(degrees) | δmin(degrees) | δmax(degrees) |State
% δ_angles are initialized to 0 degrees
Lines_matrix = [1   2   0.01938    0.05917   0.0528   0  0  0  0  1;
                1   5   0.05403    0.22304   0.0492   0  0  0  0  1;  
                2   3   0.04699    0.19797   0.0438   0  0  0  0  1;
                2   4   0.05811    0.17632   0.034    0  0  0  0  1;
                2   5   0.05695    0.17388   0.0346   0  0  0  0  1;
                3   4   0.06701    0.17103   0.0128   0  0  0  0  1;
                4   5   0.01335    0.04211   0        0  0  0  0  1;
                4   7   0          0.20912   0        0  0  0  0  1;
                4   9   0          0.55618   0        0  0  0  0  1;
                5   6   0          0.25202   0        0  0  0  0  1;
                6   11  0.09498    0.1989    0        0  0  0  0  1;
                6   12  0.12291    0.25581   0        0  0  0  0  1;
                6   13  0.06615    0.13027   0        0  0  0  0  1;
                7   8   0          0.17615   0        0  0  0  0  1;
                7   9   0          0.11001   0        0  0  0  0  1;
                9   10  0.03181    0.0845    0        0  0  0  0  1;
                9   14  0.12711    0.27038   0        0  0  0  0  1;
                10  11  0.08205    0.19207   0        0  0  0  0  1;
                12  13  0.22092    0.19988   0        0  0  0  0  1;
                13  14  0.17093    0.34802   0        0  0  0  0  1];

                
% Construction cost of lines 
% From | To | €
Cost_of_lines = [ 1  2    A
                  2  3    0.7*A
                  2  4    0.5*A
                  9  10   0.2*A
                  9  14   0.35*A
                  10 11   0.15*A
                 ];

% Data of Generators
% bus | P_G(MW) | Q_G(MVAr) | P_G,max(MW) | P_G,min(MW) | Q_G_max(MVAr)| Q_G_min(MVAr) |V_G(p.u.) | State
% Q_G is initialized to 0
Generators_matrix = [1 430 0 650 0  200  -100 1.06 1;
                     2 155 0 250 0  100  -40  1.06 1;
                     3   0 0 100 0  60   -12  1.03 1;
                     6   0 0 100 0  30   -6   1.07 1;
                     8   0 0 100 0  40   -6   1.06 1];

% Operation cost of Generators
% n | c_2 | c_1 | c_0
Cost_matrix = [3 0.0091044 63.1526 2247.9968
               3 0.024254 66.4318 688.9537 ];


% Calculate Ybus from the function calculate_Ybus()
[Y_Final, Z, B] = calculate_Ybus(Lines_matrix,Bus_matrix_v);  


% Gauss Seidel for initial Data of Power System
[not_used_1, not_used_2, P_total_Loads, P_total_Losses] = GaussSeidel(Bus_matrix_v,Bus_matrix_v_plus_1,Lines_matrix,Generators_matrix, Cost_matrix, Y_Final, Z, B);

L_for_Lagrange = P_total_Loads + P_total_Losses;

% call Lagrange function
[P_G1_lagrance, P_G2_lagrance, P_grid] = Lagrange_method(L_for_Lagrange, Cost_matrix, Generators_matrix);

% Gauss Seidel again for P_G1_lagrance, P_G2_lagrance
helping_Generators_matrix = Generators_matrix;   % helping_Generators_matrix will have the correct Data for Gauss Seidel
helping_Generators_matrix(1, 2) = P_G1_lagrance;   
helping_Generators_matrix(2, 2) = P_G2_lagrance;   

[helping1_Bus_matrix, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(Bus_matrix_v, Bus_matrix_v_plus_1, Lines_matrix, helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);

flag_Vout_of_bounts = 0; 

% Check V limits from Bus_matrix
for i = 1 : 1 : length(helping1_Bus_matrix)      
    % Vmax = helping1_Bus_matrix(i, 8)
    % Vmin = helping1_Bus_matrix(i, 9)
    if(helping1_Bus_matrix(i, 6) <= helping1_Bus_matrix(i, 8) && helping1_Bus_matrix(i, 6) >= helping1_Bus_matrix(i, 9)) 
        % V(i) is in limits    
        flag_Vout_of_bounts = flag_Vout_of_bounts +1;           
    else         
        disp('V out of bounds')     
    end
end 

if (flag_Vout_of_bounts == length(helping1_Bus_matrix))
    newRow = [helping_Generators_matrix(1, 2) helping_Generators_matrix(2, 2) P_grid L_for_Lagrange 0 0 helping1_Bus_matrix(1,6) helping1_Bus_matrix(2,6)];
    matrix_with_accepted_data = [matrix_with_accepted_data; newRow];   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%  Change P_G2 and V1,V2
helping_Generators_matrix = Generators_matrix;          % initialize helping_Generators_matrix[] with Generators_matrix data
helping1_Bus_matrix = Bus_matrix_v;
helping1_Bus_matrix_next_step = Bus_matrix_v_plus_1;

P_G2min = helping_Generators_matrix(2, 5);
P_G2max = helping_Generators_matrix(2, 4);
for  increasingPG2 = P_G2min:2.5:P_G2max  
    for v = 1.04:0.01:1.08
        helping1_Bus_matrix(1:2, 6)= v;                 % Put V from iteration
        helping1_Bus_matrix_next_step(1:2, 6)= v;
        helping1_Bus_matrix(1:2, 8)= v + v*0.06;         % Define new limits according to new V from iteration
        helping1_Bus_matrix_next_step(1:2, 8)= helping1_Bus_matrix(1:2, 8);
        helping1_Bus_matrix(1:2, 9)= v - v*0.06;
        helping1_Bus_matrix_next_step(1:2, 9)= helping1_Bus_matrix(1:2, 9);
        %helping_Generators_matrix = Generators_matrix;
        helping_Generators_matrix(1:2, 8) = v;
        helping_Generators_matrix(2, 2) = increasingPG2;    % initialize helping_Generators_matrix with new P_G2
        % Gauss Seidel for initial Data of Power System
        [not_used_1, not_used_2, P_total_Loads, P_total_Losses] = GaussSeidel(helping1_Bus_matrix,helping1_Bus_matrix_next_step,Lines_matrix,helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);
        L_for_Lagrange = P_total_Loads + P_total_Losses;
        % call Lagrange function
        [P_G1_lagrance, P_G2_lagrance, P_grid] = Lagrange_method(L_for_Lagrange, Cost_matrix, Generators_matrix);
        % Gauss Seidel again for P_G1_lagrance, P_G2_lagrance
        helping_Generators_matrix(1, 2) = P_G1_lagrance;   
        helping_Generators_matrix(2, 2) = P_G2_lagrance;   
        [helping1_Bus_matrix, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(helping1_Bus_matrix, helping1_Bus_matrix_next_step, Lines_matrix, helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);
        flag_Vout_of_bounts = 0; 
        % Check V limits from Bus_matrix
        for i = 1 : 1 : length(helping1_Bus_matrix)      
            % Vmax = helping1_Bus_matrix(i, 8)
            % Vmin = helping1_Bus_matrix(i, 9)
            if(helping1_Bus_matrix(i, 6) <= helping1_Bus_matrix(i, 8) && helping1_Bus_matrix(i, 6) >= helping1_Bus_matrix(i, 9)) 
                % V(i) is in limits    
                flag_Vout_of_bounts = flag_Vout_of_bounts +1;           
            else         
                disp('V out of bounds')     
            end
        end 
        if (flag_Vout_of_bounts == length(helping1_Bus_matrix))
            newRow = [helping_Generators_matrix(1, 2) helping_Generators_matrix(2, 2) P_grid L_for_Lagrange 0 0 helping1_Bus_matrix(1,6) helping1_Bus_matrix(2,6)];
            matrix_with_accepted_data = [matrix_with_accepted_data; newRow]; 
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%     Change P_G1 and V1,V2
helping_Generators_matrix = Generators_matrix;          % initialize helping_Generators_matrix[] with Generators_matrix data
helping1_Bus_matrix = Bus_matrix_v;
helping1_Bus_matrix_next_step = Bus_matrix_v_plus_1;

P_G1min = helping_Generators_matrix(1, 5);
P_G1max = helping_Generators_matrix(1, 4);
for  increasingPG1 = P_G1min:2.5:P_G1max  
    for v = 1.04:0.01:1.08
        helping1_Bus_matrix(1:2, 6)= v;                 % Put V from iteration
        helping1_Bus_matrix_next_step(1:2, 6)= v;
        helping1_Bus_matrix(1:2, 8)= v + v*0.06;         % Define new limits according to new V from iteration
        helping1_Bus_matrix_next_step(1:2, 8)= helping1_Bus_matrix(1:2, 8);
        helping1_Bus_matrix(1:2, 9)= v - v*0.06;
        helping1_Bus_matrix_next_step(1:2, 9)= helping1_Bus_matrix(1:2, 9);
        %helping_Generators_matrix = Generators_matrix;
        helping_Generators_matrix(1:2, 8) = v;
        helping_Generators_matrix(1, 2) = increasingPG1;    % initialize helping_Generators_matrix with new P_G1
        % Gauss Seidel for initial Data of Power System
        [not_used_1, not_used_2, P_total_Loads, P_total_Losses] = GaussSeidel(helping1_Bus_matrix,helping1_Bus_matrix_next_step,Lines_matrix,helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);
        L_for_Lagrange = P_total_Loads + P_total_Losses;
        % call Lagrange function
        [P_G1_lagrance, P_G2_lagrance, P_grid] = Lagrange_method(L_for_Lagrange, Cost_matrix, Generators_matrix);
        % Gauss Seidel again for P_G1_lagrance, P_G2_lagrance
        helping_Generators_matrix(1, 2) = P_G1_lagrance;   
        helping_Generators_matrix(2, 2) = P_G2_lagrance;   
        [helping1_Bus_matrix, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(helping1_Bus_matrix, helping1_Bus_matrix_next_step, Lines_matrix, helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);
        flag_Vout_of_bounts = 0; 
        % Check V limits from Bus_matrix
        for i = 1 : 1 : length(helping1_Bus_matrix)      
            % Vmax = helping1_Bus_matrix(i, 8)
            % Vmin = helping1_Bus_matrix(i, 9)
            if(helping1_Bus_matrix(i, 6) <= helping1_Bus_matrix(i, 8) && helping1_Bus_matrix(i, 6) >= helping1_Bus_matrix(i, 9)) 
                % V(i) is in limits    
                flag_Vout_of_bounts = flag_Vout_of_bounts +1;           
            else         
                disp('V out of bounds')     
            end
        end 
        if (flag_Vout_of_bounts == length(helping1_Bus_matrix))
            newRow = [helping_Generators_matrix(1, 2) helping_Generators_matrix(2, 2) P_grid L_for_Lagrange 0 0 helping1_Bus_matrix(1,6) helping1_Bus_matrix(2,6)];
            matrix_with_accepted_data = [matrix_with_accepted_data; newRow]; 
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%     Change P_G2
helping_Generators_matrix = Generators_matrix;          % initialize helping_Generators_matrix[] with Generators_matrix data
P_G2min = helping_Generators_matrix(2, 5);
P_G2max = helping_Generators_matrix(2, 4);

for  increasingPG2 = P_G2min:2.5:P_G2max    

    helping_Generators_matrix(2, 2) = increasingPG2;    % initialize helping_Generators_matrix with new P_G2

    % Gauss Seidel for initial Data of Power System
    [not_used_1, not_used_2, P_total_Loads, P_total_Losses] = GaussSeidel(Bus_matrix_v,Bus_matrix_v_plus_1,Lines_matrix,helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);

    L_for_Lagrange = P_total_Loads + P_total_Losses;

    

    % call Lagrange function
    [P_G1_lagrance, P_G2_lagrance, P_grid] = Lagrange_method(L_for_Lagrange, Cost_matrix, Generators_matrix);

    % Gauss Seidel again for P_G1_lagrance, P_G2_lagrance
    helping_Generators_matrix(1, 2) = P_G1_lagrance;   
    helping_Generators_matrix(2, 2) = P_G2_lagrance;   

    % o helping_Generators_matrix exei tis sostes times? Afou ta P_G1_lagrance ypologizontai swsta
    [helping1_Bus_matrix, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(Bus_matrix_v, Bus_matrix_v_plus_1, Lines_matrix, helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);

    flag_Vout_of_bounts = 0; 

    % Check V limits from Bus_matrix
    for i = 1 : 1 : length(helping1_Bus_matrix)      
        % Vmax = helping1_Bus_matrix(i, 8)
        % Vmin = helping1_Bus_matrix(i, 9)
        if(helping1_Bus_matrix(i, 6) <= helping1_Bus_matrix(i, 8) && helping1_Bus_matrix(i, 6) >= helping1_Bus_matrix(i, 9)) 
            % V(i) is in limits    
            flag_Vout_of_bounts = flag_Vout_of_bounts +1;           
        else         
            disp('V out of bounds')     
        end
    end 

    if (flag_Vout_of_bounts == length(helping1_Bus_matrix))
        newRow = [helping_Generators_matrix(1, 2) helping_Generators_matrix(2, 2) P_grid L_for_Lagrange 0 0 helping1_Bus_matrix(1,6) helping1_Bus_matrix(2,6)];
        matrix_with_accepted_data = [matrix_with_accepted_data; newRow]; 
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%     Change P_G1
helping_Generators_matrix = Generators_matrix;          % initialize helping_Generators_matrix[] with Generators_matrix data
P_G1min = helping_Generators_matrix(1, 5);
P_G1max = helping_Generators_matrix(1, 4);

for  increasingPG1 = P_G1min:2.5:P_G1max    

    helping_Generators_matrix(1, 2) = increasingPG1;    % initialize helping_Generators_matrix with new P_G1

    % Gauss Seidel for initial Data of Power System
    [not_used_1, not_used_2, P_total_Loads, P_total_Losses] = GaussSeidel(Bus_matrix_v,Bus_matrix_v_plus_1,Lines_matrix,helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);

    L_for_Lagrange = P_total_Loads + P_total_Losses;

    % call Lagrange function
    [P_G1_lagrance, P_G2_lagrance, P_grid] = Lagrange_method(L_for_Lagrange, Cost_matrix, Generators_matrix);

    % Gauss Seidel again for P_G1_lagrance, P_G2_lagrance
    helping_Generators_matrix(1, 2) = P_G1_lagrance;   
    helping_Generators_matrix(2, 2) = P_G2_lagrance;   

    % o helping_Generators_matrix exei tis sostes times? Afou ta P_G1_lagrance ypologizontai swsta
    [helping1_Bus_matrix, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(Bus_matrix_v, Bus_matrix_v_plus_1, Lines_matrix, helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);

    flag_Vout_of_bounts = 0; 

    % Check V limits from Bus_matrix
    for i = 1 : 1 : length(helping1_Bus_matrix)      
        % Vmax = helping1_Bus_matrix(i, 8)
        % Vmin = helping1_Bus_matrix(i, 9)
        if(helping1_Bus_matrix(i, 6) <= helping1_Bus_matrix(i, 8) && helping1_Bus_matrix(i, 6) >= helping1_Bus_matrix(i, 9)) 
            % V(i) is in limits    
            flag_Vout_of_bounts = flag_Vout_of_bounts +1;           
        else         
            disp('V out of bounds')     
        end
    end 

    if (flag_Vout_of_bounts == length(helping1_Bus_matrix))
        newRow = [helping_Generators_matrix(1, 2) helping_Generators_matrix(2, 2) P_grid L_for_Lagrange 0 0 helping1_Bus_matrix(1,6) helping1_Bus_matrix(2,6)];
        matrix_with_accepted_data = [matrix_with_accepted_data; newRow]; 
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Change V1, V2
helping_Generators_matrix = Generators_matrix;
helping1_Bus_matrix = Bus_matrix_v;
helping1_Bus_matrix_next_step = Bus_matrix_v_plus_1;

 for v = 1.04:0.01:1.08
     
     helping1_Bus_matrix(1:2, 6)= v;                 % Put V from iteration
     helping1_Bus_matrix_next_step(1:2, 6)= v;
     helping1_Bus_matrix(1:2, 8)= v + v*0.06;         % Define new limits according to new V from iteration
     helping1_Bus_matrix_next_step(1:2, 8)= helping1_Bus_matrix(1:2, 8);
     helping1_Bus_matrix(1:2, 9)= v - v*0.06;
     helping1_Bus_matrix_next_step(1:2, 9)= helping1_Bus_matrix(1:2, 9);
     helping_Generators_matrix(1:2,8) = v;
     
     
     % Gauss Seidel for initial Data of Power System
     [not_used_1, not_used_2, P_total_Loads, P_total_Losses] = GaussSeidel(helping1_Bus_matrix,helping1_Bus_matrix_next_step,Lines_matrix,helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);
 
     L_for_Lagrange = P_total_Loads + P_total_Losses;     
 
     % call Lagrange function
     [P_G1_lagrance, P_G2_lagrance, P_grid] = Lagrange_method(L_for_Lagrange, Cost_matrix, Generators_matrix);
 
     % Gauss Seidel again for P_G1_lagrance, P_G2_lagrance
     helping_Generators_matrix(1, 2) = P_G1_lagrance;   
     helping_Generators_matrix(2, 2) = P_G2_lagrance;   
 
     % o helping_Generators_matrix exei tis sostes times? Afou ta P_G1_lagrance ypologizontai swsta
     [helping1_Bus_matrix, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(helping1_Bus_matrix, helping1_Bus_matrix_next_step, Lines_matrix, helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);
 
     flag_Vout_of_bounts = 0; 
 
     % Check V limits from Bus_matrix
     for i = 1 : 1 : length(helping1_Bus_matrix)      
         % Vmax = helping1_Bus_matrix(i, 8)
         % Vmin = helping1_Bus_matrix(i, 9)
         if(helping1_Bus_matrix(i, 6) <= helping1_Bus_matrix(i, 8) && helping1_Bus_matrix(i, 6) >= helping1_Bus_matrix(i, 9)) 
             % V(i) is in limits    
             flag_Vout_of_bounts = flag_Vout_of_bounts +1;           
         else         
             disp('V out of bounds')     
         end
     end 
 
     if (flag_Vout_of_bounts == length(helping1_Bus_matrix))
         newRow = [helping_Generators_matrix(1, 2) helping_Generators_matrix(2, 2) P_grid L_for_Lagrange 0 0 helping1_Bus_matrix(1,6) helping1_Bus_matrix(2,6)];
         matrix_with_accepted_data = [matrix_with_accepted_data; newRow]; 
     end
 end

%%%%%%%%%%%%%%   All the possible cases
% helping_Generators_matrix = Generators_matrix;          % initialize helping_Generators_matrix[] with Generators_matrix data
% helping1_Bus_matrix = Bus_matrix_v;
% helping1_Bus_matrix_next_step = Bus_matrix_v_plus_1;
% 
% P_G1min = helping_Generators_matrix(1, 5);
% P_G1max = helping_Generators_matrix(1, 4);
% P_G2min = helping_Generators_matrix(2, 5);
% P_G2max = helping_Generators_matrix(2, 4);
% 
% for  increasingPG1 = P_G1min:2.5:P_G1max 
%     helping_Generators_matrix(1, 2) = increasingPG1;        % initialize helping_Generators_matrix with new P_G1
%     for increasingPG2 = P_G2min:2.5:P_G2max
%         helping_Generators_matrix(2, 2) = increasingPG2;    % initialize helping_Generators_matrix with new P_G2
%         for v1 = 1.04:0.01:1.08
%             helping1_Bus_matrix(1, 6)= v1;                 % Put V from iteration
%             helping1_Bus_matrix_next_step(1, 6)= v1;
%             helping1_Bus_matrix(1, 8)= v1 + v1*0.06;         % Define new limits according to new V from iteration
%             helping1_Bus_matrix_next_step(1, 8)= helping1_Bus_matrix(1, 8);
%             helping1_Bus_matrix(1, 9)= v1 - v1*0.06;
%             helping1_Bus_matrix_next_step(1, 9)= helping1_Bus_matrix(1, 9);
%             helping_Generators_matrix(1, 8) = v1;
%             for v2 = 1.04:0.01:1.08
%                 helping1_Bus_matrix(2, 6)= v2;                 % Put V from iteration
%                 helping1_Bus_matrix_next_step(2, 6)= v2;
%                 helping1_Bus_matrix(2, 8)= v2 + v2*0.06;         % Define new limits according to new V from iteration
%                 helping1_Bus_matrix_next_step(2, 8)= helping1_Bus_matrix(2, 8);
%                 helping1_Bus_matrix(2, 9)= v2 - v2*0.06;
%                 helping1_Bus_matrix_next_step(2, 9)= helping1_Bus_matrix(2, 9);
%                 helping_Generators_matrix(2, 8) = v2;
%                 % Gauss Seidel for initial Data of Power System
%                 [not_used_1, not_used_2, P_total_Loads, P_total_Losses] = GaussSeidel(helping1_Bus_matrix,helping1_Bus_matrix_next_step,Lines_matrix,helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);
%                 L_for_Lagrange = P_total_Loads + P_total_Losses;
%                 % call Lagrange function
%                 [P_G1_lagrance, P_G2_lagrance, P_grid] = Lagrange_method(L_for_Lagrange, Cost_matrix, Generators_matrix);
%                 % Gauss Seidel again for P_G1_lagrance, P_G2_lagrance
%                 helping_Generators_matrix(1, 2) = P_G1_lagrance;   
%                 helping_Generators_matrix(2, 2) = P_G2_lagrance;   
%                 [helping1_Bus_matrix, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(helping1_Bus_matrix, helping1_Bus_matrix_next_step, Lines_matrix, helping_Generators_matrix, Cost_matrix, Y_Final, Z, B);
%                 flag_Vout_of_bounts = 0; 
%                 % Check V limits from Bus_matrix
%                 for i = 1 : 1 : length(helping1_Bus_matrix)      
%                     % Vmax = helping1_Bus_matrix(i, 8)
%                     % Vmin = helping1_Bus_matrix(i, 9)
%                     if(helping1_Bus_matrix(i, 6) <= helping1_Bus_matrix(i, 8) && helping1_Bus_matrix(i, 6) >= helping1_Bus_matrix(i, 9)) 
%                         % V(i) is in limits    
%                         flag_Vout_of_bounts = flag_Vout_of_bounts +1;           
%                     else         
%                         disp('V out of bounds')     
%                     end
%                 end 
%                 if (flag_Vout_of_bounts == length(helping1_Bus_matrix))
%                     newRow = [helping_Generators_matrix(1, 2) helping_Generators_matrix(2, 2) P_grid L_for_Lagrange 0 0 helping1_Bus_matrix(1,6) helping1_Bus_matrix(2,6)];
%                     matrix_with_accepted_data = [matrix_with_accepted_data; newRow]; 
%                 end
%             end
%         end
%     end    
% end


%% Calculate F(PG) = c2 PG^2 + c1 PG + c0               (where c2 = a, c1 = b, c0 = c)
a1 = Cost_matrix(1, 2);
b1 = Cost_matrix(1, 3);
c1 = Cost_matrix(1, 4);

a2 = Cost_matrix(2, 2);
b2 = Cost_matrix(2, 3);
c2 = Cost_matrix(2, 4);

for i = 1:1:length(matrix_with_accepted_data)
    P_G1 = matrix_with_accepted_data(i, 1);
    P_G2 = matrix_with_accepted_data(i, 2);

    % Withoud Pgrid
    % F_PG = F_PG1 + F_PG2
    F_PG1 = a1 * (P_G1^2) + b1 * P_G1 + c1;
    F_PG2 = a2 * (P_G2^2) + b2 * P_G2 + c2;
    F_PG = F_PG1 + F_PG2;

    matrix_with_accepted_data(i, 5) = F_PG;
end


%% Calculate F_PG_with_Pgrid 
new_P1 = (lamda_Grid - b1) / (2*a1); 
new_P2 = (lamda_Grid - b2) / (2*a2);

F_PG1 = a1 * (new_P1^2) + b1 * new_P1 + c1;         % F_PG1 is overwritten, but we don't need the previous value any more
F_PG2 = a2 * (new_P2^2) + b2 * new_P2 + c2;

for i = 1:1:length(matrix_with_accepted_data)
    P_Grid_estimated = matrix_with_accepted_data(i, 4) - (new_P1 + new_P2);         % P_grid = L - (P1 + P2)
    F_PG_with_Pgrid = F_PG1 + F_PG2 + P_Grid_estimated * lamda_Grid;

    matrix_with_accepted_data(i, 6) = F_PG_with_Pgrid;
end    

[min_cost, index1] = min(matrix_with_accepted_data(:, 5));
disp( 'Minimum cost without buying or selling from other Power Systems')
disp(min_cost)
disp('Optimal PG1 and PG2');
disp(matrix_with_accepted_data(index1,1:2))
disp('Optimal VG1 and VG2');
disp(matrix_with_accepted_data(index1,7:8))

[min_cost_grid, index2] = min(matrix_with_accepted_data(:, 6));
disp( 'Minimum cost with buying or selling from other Power Systems with 72 €/MWh')
disp(min_cost_grid)
disp('Optimal PG1 and PG2');
disp(matrix_with_accepted_data(index2,1:2))
disp('Optimal VG1 and VG2');
disp(matrix_with_accepted_data(index2,7:8))


% matrix_with_accepted_data
% Format: P_G1 | P_G2 | P_grid |L_for_Lagrange| F(P_G) | F_optimal με αγο�?οπωλησία |




%    % Cost Matrix format
   %  n   |  C2  |   C1   |  C0    
   %      |   a  |   b    |   c 












