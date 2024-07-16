% Project 2 - Omada 8
% Planning - planning.m
clear all;
clc;
% Initialize Data of Power System

S_B = 100;          % Power Base = 100MVA
A = 1000;              % cost parameter
years = 15;
% layout :
% Combination that matches with row(s) ( aka line(s) ) from LinesMatrix | Year of construction |  F(Pg) for the year 5 | F(Pg) for the year 10 | F(Pg) for the year 15 | F(lines) | Cost = F(Pg) + F(lines)  
results_matrix = [];


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
           
n_bus = length(Bus_matrix_v);


% Calculate Ybus from the function planning_Ybus()
[Y_Final, Z_initial, B_initial] = planning_Ybus(Lines_matrix,n_bus); 

all_possible_line_comb = possible_Line_Additions();         % this is an array with all the possible line combinations
[rows_comb,col_comb]= size(all_possible_line_comb);

helping_matrix = Bus_matrix_v;
% Loads_for_every_year = zeros(15,3);   % First row for year 0[  P_L (MW) | Q_L (MVAr) | B_s (MVAr)]

for i = 1 : 1 : years
    helping_matrix(:, 3:5) = 1.015 * helping_matrix(:, 3:5);    

    if(i==5)
        Bus_matrix_v_for_Load_year_5 = helping_matrix;          % Load at year 5 -> Bus_matrix_v_for_Load_year_5
        Bus_matrix_v_plus_1_for_Load_year_5 = helping_matrix;
    end

    if(i==10)
        
        Bus_matrix_v_for_Load_year_10 = helping_matrix;          % Load at year 10 -> Bus_matrix_v_for_Load_year_10
        Bus_matrix_v_plus_1_for_Load_year_10 = helping_matrix;
    end

    if(i==15)
        Bus_matrix_v_for_Load_year_15 = helping_matrix;         % Load at year 15 -> Bus_matrix_v_for_Load_year_15
        Bus_matrix_v_plus_1_for_Load_year_15 = helping_matrix;
    end    
end

helping_Lines_matrix = Lines_matrix;
helping_Generators_matrix = Generators_matrix;

% Calculate F(PG) = c2 PG^2 + c1 PG + c0               (where c2 = a, c1 = b, c0 = c)
a1 = Cost_matrix(1, 2);
b1 = Cost_matrix(1, 3);
c1 = Cost_matrix(1, 4);

a2 = Cost_matrix(2, 2);
b2 = Cost_matrix(2, 3);
c2 = Cost_matrix(2, 4);

cases_of_year_5 = zeros(rows_comb,1);
cases_of_year_10 = zeros(rows_comb,1);
cases_of_year_15 = zeros(rows_comb,1);

% Start checking all cases of constructing Lines in years 5, 10, 15
for k = 1 : 1 : length(all_possible_line_comb)  % for year 5
    helping_Lines_matrix = Lines_matrix;

    initializeMatr = zeros(1,1);          % [0]
    rows_to_add_year_5 = initializeMatr;

    % Decode which line or lines of Lines matrix should be added
    [values] = possible_Line_Additions();
    cases_of_year_5 = values(:,k);      

    for p = 1:1:rows_comb
        if(cases_of_year_5(p) == 0)
            break;
        else
            % keep all the lines to be added (coded as a number that declares the row of LinesMatrix) in an array 
            % array has only one row, and columns ara from 0 - 6   
            rows_to_add_year_5 =  [ rows_to_add_year_5 cases_of_year_5(p)];         % format : [ 1 2] or [ 1 2 16 17 18]    
        end
    end

    % Handle case of empty implied index       
    if (~isequal(rows_to_add_year_5, zeros(1,1)))
        % Add to helping_Lines_matrix the row we want (that means the line we want)
        for iter1 = 2:1:length(rows_to_add_year_5)
            helping_Lines_matrix = [helping_Lines_matrix;   Lines_matrix(rows_to_add_year_5(1, iter1), :)   ]; % increase helping_Lines_matrix with that line
        end 
    end     

    % calculate the new Ybus    
    [Y_bus_for_year_x, Z, B] = planning_Ybus(helping_Lines_matrix,n_bus);

    % Gauss Seidel -> for 5th year will have the correct total load
    [Bus_matrix_v_plus_1, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(Bus_matrix_v_for_Load_year_5, Bus_matrix_v_plus_1_for_Load_year_5, helping_Lines_matrix, Generators_matrix, Cost_matrix, Y_bus_for_year_x, Z, B);
    
    disp('Loads')
    disp(P_total_Loads)

    disp('Losses')
    disp(P_total_Losses)

    V_check_flag_year5 = 0;
    P_check_flag_year5 = 0;

    % Check for this case if limits at V and P are ok fro year 5
    if(P_total_Losses < 0.08*P_total_Loads)                     % Check P losses limits
        P_check_flag_year5 = 1;
        disp('Losses Okay')
    else
        P_check_flag_year5 = 0;
        disp('Losses not okay')
    end     

    % V limits check for year 5
    [flag_Vout_of_bounts] = vol_limits(Bus_matrix_v_plus_1);     % check V limits for all busses for year 5

    if (flag_Vout_of_bounts == length(Bus_matrix_v_plus_1))
        V_check_flag_year5 = 1;
        disp('Voltages okay');  
    else
        V_check_flag_year5 = 0; 
    end     
    
    V_check_flag = 0;
    P_check_flag = 0;
    helping_matrix = Bus_matrix_v;

    % Check for this case if limits at V and P are ok fro years 0 to 4
    for iter1 = 1:1:4           % [1, 2, 3, 4, 5]           these are 5 cases but 1 is for year 0

        helping_matrix(:, 3:5) = (1.015)^(iter1-1) * Bus_matrix_v(:, 3:5);    % calculate             

        % For years 0-4 no lines have been added, so only the Load has changed

        % Gauss-Seidel for each year, in order to check the results
        [Bus_matrix_v_plus_1, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(helping_matrix, helping_matrix, Lines_matrix, Generators_matrix, Cost_matrix, Y_Final, Z_initial, B_initial);
        
          disp('Loads')
          disp(P_total_Loads)

          disp('Losses')
          disp(P_total_Losses)
        % P limits check for year (1ter1 - 1)
        if(P_total_Losses < 0.12*P_total_Loads)                     % Check P losses limits
            P_check_flag = 1;
            disp('Losses Okay')
        else
            P_check_flag = 0;
            disp('Losses not okay')
        end     

        % V limits check for year (1ter1 - 1)
        [flag_Vout_of_bounts] = vol_limits(Bus_matrix_v_plus_1);     % check V limits for all busses for year 5

        if (flag_Vout_of_bounts == length(Bus_matrix_v_plus_1))
            V_check_flag = 1;
            disp('Voltages okay');  
        else
            V_check_flag = 0; 
        end             

    end

    % if limits are ok continue to the next for-loops, else if not ok go to next k
    if (V_check_flag_year5 == 1 && P_check_flag_year5 == 1 && V_check_flag == 1 && P_check_flag == 1 )          % both okay for years 0 - 5
 
        P_G1_year5 = helping_Generators_matrix(1, 2);
        P_G2_year5 = helping_Generators_matrix(2, 2);

        F_PG1_year5 = a1 * (P_G1_year5^2) + b1 * P_G1_year5 + c1;
        F_PG2_year5 = a2 * (P_G2_year5^2) + b2 * P_G2_year5 + c2;

        F_PG_year5 = F_PG1_year5 + F_PG2_year5;    
        
        % Calculate Flines
        F_cost_for_lines_only = 0;              

        for iter1 = 1:1:length(rows_to_add_year_5)
            if(rows_to_add_year_5(1,iter1) == 1)                        
                F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(1, 3);
            elseif(rows_to_add_year_5(1,iter1) == 3)
                F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(2, 3);
            elseif(rows_to_add_year_5(1,iter1) == 4)
                F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(3, 3);
            elseif(rows_to_add_year_5(1,iter1) == 16)
                F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(4, 3);
            elseif(rows_to_add_year_5(1,iter1) == 17)
                F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(5, 3);
            elseif(rows_to_add_year_5(1,iter1) == 18)
                F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(6, 3);
            end

        end

        % Calculate Cost
        F_TOTAL = F_PG_year5 + F_cost_for_lines_only;

        
        % Store to matrix: layout of matrix ->
        % Combination that matches with row(s) ( aka line(s) ) from LinesMatrix | Year of construction |  F(Pg) for the year 5 | F(Pg) for the year 10 | F(Pg) for the year 15 | F(lines) | Cost = F(Pg) + F(lines)
        all_cases = cases_of_year_5(:)';
        
        one_row_added = {all_cases, 5, F_PG_year5, 0, 0, F_cost_for_lines_only, F_TOTAL};               % Create the row with the accepted data 
        results_matrix = [ results_matrix; one_row_added];                    % Add row to results_matrix

    
        % nested check for year 5 and year 10 (year 5 wight have 0 lines, so the case of line additions only at year 10 is included)    
        for l = 1 : 1 : length(all_possible_line_comb)  % for year 10
            helping_Lines_matrix = Lines_matrix;
                        
            initializeMatr = zeros(1,1);          % [0]
            rows_to_add_year_10 = initializeMatr;

            % Decode which line or lines of Lines matrix should be added
            [values] = possible_Line_Additions();
            cases_of_year_10 = values(:,l);             % here iteeration value changes from k to l

            for p = 1:1:rows_comb
                if(cases_of_year_10(p) == 0)
                    break;
                else
                    % keep all the lines to be added (coded as a number that declares the row of LinesMatrix) in an array 
                    % array has only one row, and columns ara from 0 - 6   
                    rows_to_add_year_10 =  [ rows_to_add_year_10 cases_of_year_10(p)];         % format : [ 1 2] or [ 1 2 16 17 18]    
                end
            end

            % Handle case of empty implied index    
            
            % FOR 5TH YEAR
            %if (rows_to_add_year_5 ~= zeros(1,1) || rows_to_add_year_5 ~= zeros(1,2) || rows_to_add_year_5 ~= zeros(1,3) || rows_to_add_year_5 ~= zeros(1,4) || rows_to_add_year_5 ~= zeros(1,5) || rows_to_add_year_5 ~= zeros(1,6))
            if (~isequal(rows_to_add_year_5, zeros(1,1)))
                % Add to helping_Lines_matrix the row we want (that means the line we want)
                for iter1 = 2:1:length(rows_to_add_year_5)
                    helping_Lines_matrix = [helping_Lines_matrix;   Lines_matrix(rows_to_add_year_5(1, iter1), :)   ]; % increase helping_Lines_matrix with that line
                end 
            end
    
            helping_Lines_matrix_for_year5 = helping_Lines_matrix;
            % NOW helping_Lines_matrix HAS ALREADY THE LINES FROM YEAR 5 
            
            % FOR 10TH YEAR
            
            if (~isequal(rows_to_add_year_10, zeros(1,1)))
                % Add to helping_Lines_matrix the row we want (that means the line we want) 
                for m = 2:1:length(rows_to_add_year_10)
                    helping_Lines_matrix = [helping_Lines_matrix;   Lines_matrix(rows_to_add_year_10(1, m), :)   ]; % increase helping_Lines_matrix with that line
                end   
            end   

            % calculate the new Ybus  WITH LINES FROM YEARS 5 AND 10 
            [Y_bus_for_year_x, Z, B] = planning_Ybus(helping_Lines_matrix,n_bus);

            % Gauss Seidel -> for 10th year will have the correct total load
            [Bus_matrix_v_plus_1, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(Bus_matrix_v_for_Load_year_10, Bus_matrix_v_plus_1_for_Load_year_10, helping_Lines_matrix, Generators_matrix, Cost_matrix, Y_bus_for_year_x, Z, B);
            disp('Loads')
            disp(P_total_Loads)
  
            disp('Losses')
            disp(P_total_Losses)

            V_check_flag_year10 = 0;
            P_check_flag_year10 = 0;

            % P limits check 
            if(P_total_Losses < 0.08*P_total_Loads)                     % Check P losses limits year 10
                P_check_flag_year10 = 1;
                disp('Losses Okay')
            else
                P_check_flag_year10 = 0;
                disp('Losses not okay')
            end     

            % V limits check 
            [flag_Vout_of_bounts] = vol_limits(Bus_matrix_v_plus_1);     % check V limits for all busses for year 10

            if (flag_Vout_of_bounts == length(Bus_matrix_v_plus_1))
                V_check_flag_year10 = 1;
                disp('Voltages okay');  
            else
                V_check_flag_year10 = 0; 
            end     


            % Check limits for years 6-9
            V_check_flag = 0;
            P_check_flag = 0;
            helping_matrix = Bus_matrix_v;
        
            % Check for this case if limits at V and P are ok fro years 6 to 9
            for iter1 = 1:1:10          
        
                helping_matrix(:, 3:5) = (1.015)^(iter1-1) * Bus_matrix_v(:, 3:5);    % calculate             
        
                % For years 6 to 9 lines HAVE been added, so the Load has changed and, Ybus
                % HOWEVER LinesMatrix is the one with the lines only at the 5th year
                if(iter1 >= 6)
                    % Gauss Seidel -> for 10th year will have the correct total load
                [Bus_matrix_v_plus_1, helping_Generators_matrix, P_total_Loads, P_total_Losses] = GaussSeidel(helping_matrix, helping_matrix, helping_Lines_matrix_for_year5, Generators_matrix, Cost_matrix, Y_bus_for_year_x, Z, B);
                disp('Loads')
                disp(P_total_Loads)
    
                disp('Losses')
                disp(P_total_Losses)
                    
                    if(P_total_Losses < 0.08*P_total_Loads)                   % P limits check for year (1ter1 - 1)
                        P_check_flag = 1;
                        disp('Losses Okay')
                    else
                        P_check_flag = 0;
                        disp('Losses not okay')
                    end     
            
                    % V limits check for year (1ter1 - 1)
                    [flag_Vout_of_bounts] = vol_limits(Bus_matrix_v_plus_1);     % check V limits for all busses for year (1ter1 - 1)
            
                    if (flag_Vout_of_bounts == length(Bus_matrix_v_plus_1))
                        V_check_flag = 1;
                        disp('Voltages okay');  
                    else
                        V_check_flag = 0; 
                    end             
                end    
            end



            % if limits are ok continue to the next for-loop, else if not ok go to next k
            if (V_check_flag_year10 == 1 && P_check_flag_year10 == 1 && V_check_flag == 1 && P_check_flag == 1)          % both okay
                % calculate F(PG) 
                P_G1_year10 = helping_Generators_matrix(1, 2);
                P_G2_year10 = helping_Generators_matrix(2, 2);

                F_PG1_year10 = a1 * (P_G1_year10^2) + b1 * P_G1_year10 + c1;
                F_PG2_year10 = a2 * (P_G2_year10^2) + b2 * P_G2_year10 + c2;

                F_PG_year10 = F_PG1_year10 + F_PG2_year10 + F_PG_year5;         % add the cost from year 5  
                %F_PG_year10 = F_PG1_year10 + F_PG2_year10;                

                all_cases = [cases_of_year_5(:)', cases_of_year_10(:)'];

                % Calculate Flines
                F_cost_for_lines_only = 0;              

                for iter1 = 1:1:length(rows_to_add_year_5)
                    if(rows_to_add_year_5(1,iter1) == 1)                        
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(1, 3);
                    elseif(rows_to_add_year_5(1,iter1) == 3)
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(2, 3);
                    elseif(rows_to_add_year_5(1,iter1) == 4)
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(3, 3);
                    elseif(rows_to_add_year_5(1,iter1) == 16)
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(4, 3);
                    elseif(rows_to_add_year_5(1,iter1) == 17)
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(5, 3);
                    elseif(rows_to_add_year_5(1,iter1) == 18)
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(6, 3);
                    end

                end

                for iter1 = 1:1:length(rows_to_add_year_10)
                    if(rows_to_add_year_10(1,iter1) == 1)                        
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(1, 3);
                    elseif(rows_to_add_year_10(1,iter1) == 3)
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(2, 3);
                    elseif(rows_to_add_year_10(1,iter1) == 4)
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(3, 3);
                    elseif(rows_to_add_year_10(1,iter1) == 16)
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(4, 3);
                    elseif(rows_to_add_year_10(1,iter1) == 17)
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(5, 3);
                    elseif(rows_to_add_year_10(1,iter1) == 18)
                        F_cost_for_lines_only = F_cost_for_lines_only + Cost_of_lines(6, 3);
                    end

                end

                % Calculate Cost
                F_TOTAL = F_PG_year5 + F_cost_for_lines_only;
        

                % Store to matrix: layout of matrix ->
                % Combination that matches with row(s) ( aka line(s) ) from LinesMatrix | Year of construction |  F(Pg) for the year 5 | F(Pg) for the year 10 | F(Pg) for the year 15 | F(lines) | Cost = F(Pg) + F(lines)
                % here the first column has something like this : [3 0 0 0 0 0 1 4 0 0 0 0] And we know that the first 5 are refering to year 5 and the rest to year 10
                
                one_row_added = {all_cases, 10, F_PG_year5, F_PG_year10, 0, F_cost_for_lines_only, F_TOTAL};               % Create the row with the accepted data 
                results_matrix = [ results_matrix; one_row_added];               % Add row to results_matrix

    
                % here "cases_of_year_5" are the combinations of lines for year 5         

                % for m = 1 : 1 : length(all_possible_line_comb)  % for year 15
                %     % here k are the combinations of lines for year 5
                %     % here l are the combinations of lines for year 10

                %     calculate new Ybus - from the previous one which is for additions in year 5
                %     store Ybus to final_matrix
                %     store that combination of the forloop in final_matrix

                %     Gauss Seidel 

                %     P limits check 
                %     V limits check          
                % end
            end
        end 

    end   
end       




numerical_column = cell2mat(results_matrix(:, 7));

[min_value, min_index] = min(numerical_column);

% Display the most efficient line construction plan for the next 10 years
disp('Most efficient line construction plan for the next 10 years')

% Extract the corresponding row from results_matrix
min_row = results_matrix(min_index, :);


disp(min_row)


% Combination that matches with row(s) ( aka line(s) ) from LinesMatrix | Year of construction |  F(Pg) for the year 5 | F(Pg) for the year 10 | F(Pg) for the year 15 | F(lines) | Cost = F(Pg) + F(lines)