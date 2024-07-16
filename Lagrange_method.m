function [P_G1_lagrance, P_G2_lagrance, P_grid] = Lagrange_method(L, Cost_matrix, Generators_matrix)
% Arguments : L (total Load which is L_for_Lagrange = P_total_Loads + P_total_Losses)
%             Cost_matrix -> has c1, c2, c3 
% Outputs : P_G1_lagrance, P_G2_lagrance


    cost_sum1 = 0;
    cost_sum2 = 0;

   % Cost Matrix format
   %  n   |  C2  |   C1   |  C0    
   %      |   a  |   b    |   c 

    for i = 1:1:2
        cost_sum1 = cost_sum1 + Cost_matrix(i,3)/(2*Cost_matrix(i,2)); 
        cost_sum2 = cost_sum2 + 1/(2*Cost_matrix(i,2)); 
    end
    
    % calculate λ
    lamda = (L + cost_sum1 )/ cost_sum2;
   
    P_G1_lagrance = (lamda - Cost_matrix(1,3))/ (2 * Cost_matrix(1,2)); 
    P_G2_lagrance = (lamda - Cost_matrix(2,3))/ (2 * Cost_matrix(2,2));

    flag_P_G1_changed = 0; 
    flag_P_G2_changed = 0; 

   
    % bus | P_G(MW) | Q_G(MVAr) | P_G,max(MW) | P_G,min(MW) | Q_G_max(MVAr)| Q_G_min(MVAr) |V_G(p.u.) | State
    P_G1_min = Generators_matrix(1,5);
    P_G1_max =  Generators_matrix(1,4);
    P_G2_min = Generators_matrix(2,5);
    P_G2_max = Generators_matrix(2,4);


    if (P_G1_lagrance < P_G1_min && P_G2_lagrance > P_G2_min && P_G2_lagrance < P_G2_max)      % Only one out of bounds
        P_G1_lagrance = P_G1_min;
        P_G2_lagrance = L - P_G1_lagrance;
        flag_P_G2_changed = 1; 
    elseif (P_G1_lagrance > P_G1_max && P_G2_lagrance > P_G2_min && P_G2_lagrance < P_G2_max)
        P_G1_lagrance = P_G1_max;
        P_G2_lagrance = L - P_G1_lagrance;
        flag_P_G2_changed = 1; 
    elseif (P_G2_lagrance < P_G2_min && P_G1_lagrance > P_G1_min && P_G1_lagrance < P_G1_max)
        P_G2_lagrance = P_G2_min;
        P_G1_lagrance = L - P_G2_lagrance; 
        flag_P_G1_changed = 1; 
    elseif (P_G2_lagrance > P_G2_max && P_G1_lagrance > P_G1_min && P_G1_lagrance < P_G1_max)                                        
        P_G2_lagrance = P_G2_max;
        P_G1_lagrance = L - P_G2_lagrance; 
        flag_P_G1_changed = 1; 
    elseif (P_G1_lagrance < P_G1_min &&  P_G2_lagrance < P_G2_min)              % both out of bounds
        P_G1_lagrance = P_G1_min;
        P_G2_lagrance = P_G2_min;
    elseif (P_G1_lagrance > P_G1_max &&  P_G2_lagrance > P_G2_max)
        P_G1_lagrance = P_G1_max;
        P_G2_lagrance = P_G2_max;
    elseif (P_G1_lagrance > P_G1_max && P_G2_lagrance < P_G2_min)
        P_G1_lagrance = P_G1_max;
        P_G2_lagrance = P_G2_min;
    elseif (P_G1_lagrance < P_G1_min && P_G2_lagrance > P_G2_max)
        P_G1_lagrance = P_G1_min;
        P_G2_lagrance = P_G2_max;
    else                                                                        % both ok
        disp('Everything is fine')
    end                     

    if(flag_P_G2_changed == 1)
        if (P_G2_lagrance < P_G2_min)
            P_G2_lagrance = P_G2_min;
        elseif (P_G2_lagrance > P_G2_max)
            P_G2_lagrance = P_G2_max;
        end
    elseif (flag_P_G1_changed == 1)
        if (P_G1_lagrance < P_G1_min)
            P_G1_lagrance = P_G1_min;
        elseif (P_G1_lagrance > P_G1_max)
            P_G1_lagrance = P_G1_max;
        end    
    end

    P_grid = L - (P_G1_lagrance + P_G2_lagrance);

    % P_grid > 0,  we need to buy 
    % P_grid < 0,  we have more MW available, so we sell 
    % P_grid = 0,  

    