% Project 2 - Omada 8
% 4.** are refered to formulas from the book

clear;
clc

% Initialize Data of Power System

S_B = 100;     % Power Base = 100MVA

% Bus data for step v
% bus | Type of bus | P_L (MW) | Q_L (MVAr) | B_s (MVAr) | Vm | Va | Vmax | Vmin
% Initializations for Load busses: V_m = 1 and δ_angles = Va = 0.
% Initializations for Voltage-controlled busses: V_m = Vg -> Vspec from Generators_matrix and δ_angles = Va = 0.

Bus_matrix_v = [1   3     0        0       0    1.06    0    1.1236   0.9964;       % Slack Bus
                2   2     21.7     12.7    0    1.045   0    1.1077   0.9823        % Voltage-controlled bus
                3   2     94.2     19      0    1.01    0    1.0706   0.9494;       % Voltage-controlled bus
                4   1     47.8     -3.9    0    1       0    1.06     0.94;
                5   1     7.6      1.6     0    1       0    1.06     0.94;
                6   2     11.2     7.5     0    1.07    0    1.1342   1.0058;       % Voltage-controlled bus
                7   1     0        0       0    1       0    1.06     0.94;
                8   2     0        0       0    1.09    0    1.1554   1.0246;       % Voltage-controlled bus
                9   1     29.5     16.6    19   1       0    1.06     0.94;
                10  1     9        5.8     0    1       0    1.06     0.94;
                11  1     3.5      1.8     0    1       0    1.06     0.94;
                12  1     6.1      1.6     0    1       0    1.06     0.94;
                13  1     13.5     5.8     0    1       0    1.06     0.94;
                14  1     14.9     5       0    1       0    1.06     0.94];

% Bus data for step v+1
Bus_matrix_v_plus_1 = [1   3     0        0       0    1.06    0    1.1236   0.9964;       % Slack Bus
                       2   2     21.7     12.7    0    1.045   0    1.1077   0.9823        % Voltage-controlled bus
                       3   2     94.2     19      0    1.01    0    1.0706   0.9494;       % Voltage-controlled bus
                       4   1     47.8     -3.9    0    1       0    1.06     0.94;
                       5   1     7.6      1.6     0    1       0    1.06     0.94;
                       6   2     11.2     7.5     0    1.07    0    1.1342   1.0058;       % Voltage-controlled bus
                       7   1     0        0       0    1       0    1.06     0.94;
                       8   2     0        0       0    1.09    0    1.1554   1.0246;       % Voltage-controlled bus
                       9   1     29.5     16.6    19   1       0    1.06     0.94;
                       10  1     9        5.8     0    1       0    1.06     0.94;
                       11  1     3.5      1.8     0    1       0    1.06     0.94;
                       12  1     6.1      1.6     0    1       0    1.06     0.94;
                       13  1     13.5     5.8     0    1       0    1.06     0.94;
                       14  1     14.9     5       0    1       0    1.06     0.94];


% Data of Transmission Lines
% from bus | to bus | R(p.u.) | X(p.u.)| B(p.u.)| transf.limit(p.u.)| δ(degrees) | δmin(degrees) | δmax(degrees) |State
% δ_angles are initialized to 0 degrees
Lines_matrix = [1 2 0.01938 0.05917 0.0528 0 0 0 0 1;
                1 5 0.05403 0.22304 0.0492 0 0 0 0 1;  % this is in comments for question 4.
                2 3 0.04699 0.19797 0.0438 0 0 0 0 1;
                2 4 0.05811 0.17632 0.034 0 0 0 0 1;
                2 5 0.05695 0.17388 0.0346 0 0 0 0 1;
                3 4 0.06701 0.17103 0.0128 0 0 0 0 1;
                4 5 0.01335 0.04211 0 0 0 0 0 1;
                4 7 0 0.20912 0 0 0 0 0 1;
                4 9 0 0.55618 0 0 0 0 0 1;
                5 6 0 0.25202 0 0 0 0 0 1;
                6 11 0.09498 0.1989 0 0 0 0 0 1;
                6 12 0.12291 0.25581 0 0 0 0 0 1;
                6 13 0.06615 0.13027 0 0 0 0 0 1;
                7 8 0 0.17615 0 0 0 0 0 1;
                7 9 0 0.11001 0 0 0 0 0 1;
                9 10 0.03181 0.0845 0 0 0 0 0 1;
                9 14 0.12711 0.27038 0 0 0 0 0 1;
                10 11 0.08205 0.19207 0 0 0 0 0 1;
                12 13 0.22092 0.19988 0 0 0 0 0 1;
                13 14 0.17093 0.34802 0 0 0 0 0 1];

% Data of Generators
% bus | P_G(MW) | Q_G(MVAr) | P_G,max(MW) | P_G,min(MW) | Q_G_max(MVAr)| Q_G_min(MVAr) |V_G(p.u.) | State
% Q_G isi initialized to 0
Generators_matrix = [1 232.4 0 332.4 0 10 0 1.06 1;
     2 40 0 140 0 50 -40 1.045 1;
     3 0 0 100 0 40 0 1.01 1;
     6 0 0 100 0 24 -6 1.07 1;
     8 0 0 100 0 24 -6 1.09 1];

% Operation cost of Generators
% n | c_2 | c_1 | c_0
Cost_matrix = [3 0.1 20 0
               3 0.25 20 0 ];


% Calculate Ybus
% all these are in pu
Z = Lines_matrix(:,3)+1i*Lines_matrix(:,4);  % Z = R + 1i*X
B = 1i*Lines_matrix(:,5);          % this is Ye
Y = 1./Z;      % this is not Ye
n_bus = length(Bus_matrix_v);
L1 = Lines_matrix(:,1);
L2 = Lines_matrix(:,2);

Ybus = zeros(n_bus,n_bus);

for n1 = 1:n_bus
     for n2 = 1:length(Lines_matrix)
          if(L1(n2) == n1)
               k = L2(n2);
               Ybus(n1,k) = -Y(n2);
          elseif(L2(n2) == n1)
               v1 = L1(n2);
               Ybus(n1,v1) = -Y(n2);
          end
          
          if((L1(n2) == n1) || (L2(n2) == n1))
               Ybus(n1,n1) = sum(Ybus(n1,n1)) + Y(n2) + (B(n2)/2);
          end
     end
end
Y_final = Ybus;

 gen_vol = Bus_matrix_v(:,6);

% Gauss-Seidel method
v_max = 100;                    % maximum number of iterations
e_threshold = 0.0001;            % threshold of accuracy between two consecutive values
a = 1.5;                        % accelerator factor
iteration = 0;
abs_of_DeltaV = 0;              % initialize ΔV
Delta_V_max = 0;                % initialize ΔVmax
Qi_v_plus_1 = 0;                % initialize to 0
[widthGenerators_matrix, lnght] = size(Generators_matrix);

PowerFlow_lines_i_to_j = zeros(length(Lines_matrix),2);  % From | to |  Sij
PowerFlow_lines_j_to_i = zeros(length(Lines_matrix),2);  % From | to |  Sji


for v = 0:1:v_max
     Delta_V_max = 0;                % initialize DVmax, at each iteration, reset DVmax

     for i = 2:1:length(Bus_matrix_v)


          % Calculate 2 sums to use later
          sum1 = 0;                          
         for k = 1:1:i-1
             [Vreal_v_plus_1,Vimag_v_plus_1] = pol2cart(Bus_matrix_v_plus_1(k,7),Bus_matrix_v_plus_1(k,6));
              sum1 = sum1 + Y_final(i,k) * ( Vreal_v_plus_1  + 1i *Vimag_v_plus_1 );
         end

         sum2 = 0;
         for k = i+1:1:length(Bus_matrix_v)  
             [Vreal_v,Vimag_v] = pol2cart(Bus_matrix_v(k,7),Bus_matrix_v(k,6));
              sum2 = sum2 + Y_final(i,k) * ( Vreal_v  + 1i * Vimag_v );
         end

          if(Bus_matrix_v(i,2) == 2)         % 2: equals to Voltage-controlled bus type
               
               for j = 2:1:(widthGenerators_matrix)
                    Bus_matrix_v_plus_1(i,6) = Bus_matrix_v(i,6);
                    generator_bus = Generators_matrix(j,1);
                    if(i==generator_bus)               % check if it's a Voltage-controlled bus type and has a Generator
                         % In our system, every Voltage-controlled bus is also a Gnerator bus, so there is no else case for that. 
                         % 4.43 : calculate Qi_v_plus_1
                         
                         % Calculate the 2 sums for Qi_v_plus_1
                         sum_v_plus_1 = 0;   % first sum of Qi_v_plus_1 formula
                         for k = 1:1:i-1
                              yij_real = real(Y_final(i,k));
                              yij_imag = imag(Y_final(i,k));
                              [yij_angle, yij_meter] = cart2pol(yij_real,yij_imag);
                              
                              sum_v_plus_1 = sum_v_plus_1 + abs(Bus_matrix_v_plus_1(k,6)) * yij_meter * sin(Bus_matrix_v_plus_1(k,7) - Bus_matrix_v(i,7) + yij_angle );
                         end
                         
                         sum_v = 0;   % second sum of Qi_v_plus_1 formula
                         for k = i:1:length(Bus_matrix_v)
                              yij_real = real(Y_final(i,k));
                              yij_imag = imag(Y_final(i,k));
                              [yij_angle, yij_meter] = cart2pol(yij_real,yij_imag);
                              
                              sum_v = sum_v + abs(Bus_matrix_v(k,6)) * yij_meter * sin(Bus_matrix_v(k,7) - Bus_matrix_v(i,7) + yij_angle );
                         end
                         
                         %Bs_pu = Bus_matrix_v(i,5) / S_B;
                         
                         % we changed this!!!!!!!!!!!!!
                         %Qi_v_plus_1_pu = - abs(gen_vol(i)) * (sum_v_plus_1 + sum_v) + Bs_pu;        % calculated with pu values, so the result is in pu
                         Qi_v_plus_1_pu = - abs(gen_vol(i)) * (sum_v_plus_1 + sum_v);        % calculated with pu values, so the result is in pu

                         % % turn Qi_v_plus_1, Pi into p.u.
                         % Qi_v_plus_1_pu = Qi_v_plus_1 / S_B;
                         Pi_pu = (Generators_matrix(j,2) / S_B) - ( Bus_matrix_v_plus_1(i,3)/ S_B);
                         
                                                 
                         Qi_max = Generators_matrix(j,6) - Bus_matrix_v_plus_1(i,4);
                         Qi_max_pu = Qi_max / S_B;

                         Qi_min = Generators_matrix(j,7) - Bus_matrix_v_plus_1(i,4);
                         Qi_min_pu = Qi_min / S_B;
                         
                         if(Qi_v_plus_1_pu < Qi_max_pu && Qi_v_plus_1_pu > Qi_min_pu)     % if it obeys the limits
                              % 4.45 + 4.42
                              %Bus_matrix_v(i,6) = gen_vol(i);
                              [Vreal_v,Vimag_v] = pol2cart(Bus_matrix_v(i,7),gen_vol(i));
                              %[Vreal_v,Vimag_v] = pol2cart(Bus_matrix_v(i,7),Generators_matrix(j,8));
                              
                              Vi_v = Vreal_v + 1i * Vimag_v;     % Vi,(v) = |Vi,(v)|δi,(v), but cartesian coordinates
                              
                              Vi_v_plus_1 = (1/Y_final(i,i)) * (((Pi_pu - 1i*Qi_v_plus_1_pu) / conj(Vi_v)) - sum1 - sum2);    % Vi,(v+1) = |Vi,(v+1)|δi,(v+1)
                              %[Va_v_plus_1,Vm_v_plus_1] = cart2pol(real(Vi_v_plus_1),imag(Vi_v_plus_1));
                              %[Vreal_v_plus_1,Vimag_v_plus_1] = pol2cart(Va_v_plus_1,Bus_matrix_v(i,6));
                              %Vi_v_plus_1 = Vreal_v_plus_1 + 1i * Vimag_v_plus_1;
                              
                              Vi_acc_v_plus_1 = Vi_v + a * (Vi_v_plus_1 - Vi_v);      % that is Vi,acc(v+1) = Vi,acc(v) + α(Vi(v+1) - Vi(v)) 
                              
                              real_part = real(Vi_acc_v_plus_1);           % Split complex Vi_acc_v_plus_1 to cartesian coordinates
                              imag_part = imag(Vi_acc_v_plus_1);
                             
%                               helping_var = zeros(i,1);          % this won't be used
%                               disp('This the Voltage that would have overwritten Vspec')
%                               disp(helping_var)
%                               disp('for bus ')
%                               disp(i)
%                               disp('times of iteration ')
%                               disp(v)
                              [Bus_matrix_v_plus_1(i,7), helping_var] = cart2pol(real_part,imag_part);      % Store the new accelerated values
                              % Vspec stays the same 
                                                            
                         elseif(Qi_v_plus_1_pu > Qi_max_pu)                                      % if it exceed Gi,max
                              Qi_v_plus_1 = Generators_matrix(j,6) - Bus_matrix_v_plus_1(i,4);            % Q_G - Q_L = Qi
                              Qi_v_plus_1_pu = Qi_v_plus_1 / S_B;                      % pu
                              
                              % 4.42
                              [Vreal_v,Vimag_v] = pol2cart(Bus_matrix_v(i,7),Bus_matrix_v(i,6));
                              
                              Vi_v = Vreal_v + 1i * Vimag_v;     % Vi,(v) = |Vi,(v)|δi,(v), but cartesian coordinates

                              Vi_v_plus_1 = (1/Y_final(i,i)) * (((Pi_pu - 1i*Qi_v_plus_1_pu) / conj(Vi_v)) - sum1 - sum2);    % Vi,(v+1) = |Vi,(v+1)|δi,(v+1)
                              
                              Vi_acc_v_plus_1 = Vi_v + a * (Vi_v_plus_1 - Vi_v);      % that is Vi,acc(v+1) = Vi,acc(v) + α(Vi(v+1) - Vi(v)) 
                              
                              real_part = real(Vi_acc_v_plus_1);           % Split complex Vi_acc_v_plus_1 to cartesian coordinates
                              imag_part = imag(Vi_acc_v_plus_1);

                              [Bus_matrix_v_plus_1(i,7),Bus_matrix_v_plus_1(i,6)] = cart2pol(real_part,imag_part);      % Store the new accelerated values

                              
                         else                                                                            % if it is below Gi,min
                              Qi_v_plus_1 = Generators_matrix(j,7) - Bus_matrix_v_plus_1(i,4);            % Q_G - Q_L = Qi
                              Qi_v_plus_1_pu = Qi_v_plus_1 / S_B;                      % pu
                              
                              % 4.42
                              [Vreal_v,Vimag_v] = pol2cart(Bus_matrix_v(i,7),Bus_matrix_v(i,6));
                              
                              Vi_v = Vreal_v + 1i * Vimag_v;     % Vi,(v) = |Vi,(v)|δi,(v), but cartesian coordinates

                              Vi_v_plus_1 = (1/Y_final(i,i)) * (((Pi_pu - 1i*Qi_v_plus_1_pu) / conj(Vi_v)) - sum1 - sum2);    % Vi,(v+1) = |Vi,(v+1)|δi,(v+1)
                              
                              Vi_acc_v_plus_1 = Vi_v + a * (Vi_v_plus_1 - Vi_v);      % that is Vi,acc(v+1) = Vi,acc(v) + α(Vi(v+1) - Vi(v)) 
                              
                              real_part = real(Vi_acc_v_plus_1);           % Split complex Vi_acc_v_plus_1 to cartesian coordinates
                              imag_part = imag(Vi_acc_v_plus_1);

                              [Bus_matrix_v_plus_1(i,7),Bus_matrix_v_plus_1(i,6)] = cart2pol(real_part,imag_part);      % Store the new accelerated values
                             
                         end

                         Qi_v_plus_1 = Qi_v_plus_1_pu * S_B;               % in Generators_matrix the values are not in pu

                         % turn into not pu, Bs is not in pu so we use it as it is
                         % Q_Gi = Q_Li + Qi - Bs
                         Generators_matrix(j,3) =  Qi_v_plus_1 + Bus_matrix_v_plus_1(i,4) - Bus_matrix_v(i,5);            % Store Q_G in Generators_matrix

                    end
               end
          else                               % 1: Load bus type
               Bus_matrix_v_plus_1(i,6) = Bus_matrix_v(i,6);
               Pi = - Bus_matrix_v_plus_1(i,3);        % Pi = P_G - P_L = 0 - P_L
               Qi = - Bus_matrix_v_plus_1(i,4);        % Qi = Q_G - Q_L = 0 - Q_L
               
               % pu
               Pi_pu = Pi /S_B;
               Qi_pu = Qi /S_B;
             
               % 4.42
               [Vreal_v,Vimag_v] = pol2cart(Bus_matrix_v(i,7),Bus_matrix_v(i,6));

               Vi_v = Vreal_v + 1i * Vimag_v;     % Vi,(v) = |Vi,(v)|δi,(v), but cartesian coordinates

               Vi_v_plus_1 = (1/Y_final(i,i)) * (((Pi_pu - 1i*Qi_pu) / conj(Vi_v)) - sum1 - sum2);

               Vi_acc_v_plus_1 = Vi_v + a * (Vi_v_plus_1 - Vi_v);      % that is Vi,acc(v+1) = Vi,acc(v) + α(Vi(v+1) - Vi(v)) 
                              
               real_part = real(Vi_acc_v_plus_1);           % Split complex Vi_acc_v_plus_1 to cartesian coordinates
               imag_part = imag(Vi_acc_v_plus_1);

               [Bus_matrix_v_plus_1(i,7),Bus_matrix_v_plus_1(i,6)] = cart2pol(real_part,imag_part);      % Store the new accelerated values

               
          end
          
          % Calculate |DVi,v+1|
          % disp('current \Delta V for i = ')
          % disp(i)
          abs_of_DeltaV = abs(Bus_matrix_v_plus_1(i,6) -  Bus_matrix_v(i,6));
          
          % DVmax
          if(abs_of_DeltaV >= Delta_V_max)
               Delta_V_max = abs_of_DeltaV;
          end
          
          Bus_matrix_v = Bus_matrix_v_plus_1;     % ready for the next step i, store new values
     end      % exit for loop i

     
     if(Delta_V_max <= e_threshold)          % if algorithms converges
          % % Q_G1 at Slack Bus
          % sum_G1 = 0;   % Calculate Sum of Q_G1 formula
          % for j = 1:1:length(Y_final)
          %      yij_real = real(Y_final(1,j));
          %      yij_imag = imag(Y_final(1,j));
          %      [yij_angle, yij_meter] = cart2pol(yij_real,yij_imag);
               
          %      sum_G1 = sum_G1 + ( Bus_matrix_v_plus_1(1, 6) *  Bus_matrix_v_plus_1(j, 6) * yij_meter * sin( Bus_matrix_v_plus_1(j, 7) - Bus_matrix_v_plus_1(1, 7) + yij_angle ));
          % end
          
          % Q_G1_pu = -sum_G1;                   %  Q_1 = Q_G1 - Q_L1 => (Q_L1 =0) => Q_1 = Q_G1  

          % Q_G1 = S_B * Q_G1_pu;                % turn Q_G1_pu to MVAr
          % Generators_matrix(1,3) = Q_G1;       % Store correct value to matrix


          % Calculate S_G1 at Slack Bus for Erotima 7
          % V_slackbus = Vm (meter), Va (angle)
          [V_slackbus_pu_real, V_slackbus_pu_imag] = pol2cart(Bus_matrix_v_plus_1(1,7), Bus_matrix_v_plus_1(1,6));     % V_slackbus in Cartesian  
         
          sum_S_slackbus_pu = 0;   % Calculate Sum of S_G_SlackBus formula
          for j = 1:1:length(Y_final)
               % yij_real = real(Y_final(1,j));
               % yij_imag = imag(Y_final(1,j));
               % [yij_angle, yij_meter] = cart2pol(yij_real,yij_imag);

               % All these are in pu
               [V_j_real_pu, V_j_imag_pu] = pol2cart(Bus_matrix_v_plus_1(j,7), Bus_matrix_v_plus_1(j,6));     % V_slackbus in Cartesian               
               sum_S_slackbus_pu = sum_S_slackbus_pu + (Y_final(1,j) * (V_j_real_pu+ 1i*V_j_imag_pu));
          end
          % We want conj(V_slackbus) that's why we have -1i*V_slackbus_imag
          S_G_SlackBus_pu = (V_slackbus_pu_real-1i*V_slackbus_pu_imag)* sum_S_slackbus_pu;
          
          S_G_SlackBus = S_G_SlackBus_pu * S_B;

          P_G_SlackBus = real(S_G_SlackBus);
          Q_G_SlackBus = imag(S_G_SlackBus);

          Generators_matrix(1,3) = Q_G_SlackBus;       % Store correct value to matrix

          % Compare Q_G_SlackBus with data from Generators_matrix(1,2) 
          if (P_G_SlackBus == Generators_matrix(1,2))
               disp('Data of P_G_SlackBus are the same')
          else
               disp('Data of P_G_SlackBus are NOT the` same')
          end

          % 4.30 
          from_bus = Lines_matrix(:, 1);     % Returns a vector with all the busses in the first row
          to_bus =  Lines_matrix(:, 2);      % Returns a vector with all the busses in the second row
          S_ij_pu = zeros(length(Lines_matrix),1);     % S_ij in pu
          S_ji_pu = zeros(length(Lines_matrix),1);     % S_ji in pu
          P_ij = zeros(length(Lines_matrix),1);        % P_ij
          Q_ij = zeros(length(Lines_matrix),1);        % Q_ij
          P_ji = zeros(length(Lines_matrix),1);        % P_ji
          Q_ji = zeros(length(Lines_matrix),1);        % Q_ji
          S_ij = zeros(length(Lines_matrix),1);
          S_ji = zeros(length(Lines_matrix),1);
          S_losses_total = 0; 
          
          for k = 1:1:length(Lines_matrix)            
               i = from_bus(k);                   % i has correct number of from_bus
               j = to_bus(k);                     % j has correct number of to_bus

               [Vi_real_pu,Vi_imag_pu] = pol2cart(Bus_matrix_v_plus_1(i,7),Bus_matrix_v_plus_1(i,6));
               Vi_phasor_pu = Vi_real_pu + 1i * Vi_imag_pu;            % Vi phasor pu
               
               [Vj_real_pu,Vj_imag_pu] = pol2cart(Bus_matrix_v_plus_1(j,7),Bus_matrix_v_plus_1(j,6));
               Vj_phasor_pu = Vj_real_pu + 1i * Vj_imag_pu;            % Vj phasor pu 

               % Erotima 4
               % conj() wants cartesian | Z(k) is Ze from p-equivalent of TL | B(k) is Ye from p-equivalent of TL
               S_ij_pu(k) = Vi_phasor_pu * (conj(Vi_phasor_pu) - conj(Vj_phasor_pu)) * ( 1 / conj(Z(k)) ) + ( (Vi_real_pu^2) *  (conj(B(k)) ));          % formula Sij for pu
               S_ji_pu(k) = Vj_phasor_pu * (conj(Vj_phasor_pu) - conj(Vi_phasor_pu)) * ( 1 / conj(Z(k)) ) + ( (Vj_real_pu^2) *  (conj(B(k)) ));          % formula Sji for pu

               S_ij(k) = S_ij_pu(k) * S_B;
               S_ji(k) = S_ji_pu(k) * S_B;

               P_ij(k) = real(S_ij(k));     % turn Sij to cartesian coordinates
               Q_ij(k) = imag(S_ij(k));
               
               PowerFlow_lines_i_to_j(k,1) = P_ij(k);       % Pij (MW)
               PowerFlow_lines_i_to_j(k,2) = Q_ij(k);       % Qij (MVAr)

               P_ji(k) = real(S_ji(k));     % turn Sij to cartesian coordinates
               Q_ji(k) = imag(S_ji(k));
               
               PowerFlow_lines_j_to_i(k,1) = P_ji(k);       % Pji (MW)
               PowerFlow_lines_j_to_i(k,2) = Q_ji(k);       % Qji (MVAr)

               % Erotima 5 -> katw kateytheian

               % Erotima 6
               S_losses_total = S_losses_total + S_ij(k) + S_ji(k);        % this iterates k times  
          end
          
          disp('The algorithm converges.');         
          
          % Turn all angles from rad to degrees
          angles_deg = rad2deg(Bus_matrix_v_plus_1(:,7));        
          
          % Put delta in Lines_matrix -> All angles are in degrees
          for k = 1:1:length(Lines_matrix)     
               i = from_bus(k);                   % i has correct number of from_bus
               j = to_bus(k);                     % j has correct number of to_bus

               % delta_i-delta_j (k)
               Lines_matrix(k,7) = Bus_matrix_v_plus_1(i,7) - Bus_matrix_v_plus_1(j,7);      % here the angles are in degrees
             
               delta_i_max = Bus_matrix_v_plus_1(i,7) + Bus_matrix_v_plus_1(i,7) * 0.06;
               delta_j_max = Bus_matrix_v_plus_1(j,7) + Bus_matrix_v_plus_1(j,7) * 0.06;
               
               Lines_matrix(k,9) = delta_i_max - delta_j_max;
               
               delta_i_min = Bus_matrix_v_plus_1(i,7) - Bus_matrix_v_plus_1(i,7) * 0.06;
               delta_j_min = Bus_matrix_v_plus_1(j,7) - Bus_matrix_v_plus_1(j,7) * 0.06;
               
               Lines_matrix(k,8)  = delta_i_min - delta_j_min;
          end

          % gia to erotima 6
          P_losses_ij = real(S_losses_total);    % P line losses
          Q_losses_ij = imag(S_losses_total);    % Q line losses

          % Print Results
          
          % Erotima 3
          disp('|Vi|in pu and delta_ i in degrees for each node ')  
          disp([Bus_matrix_v_plus_1(:,6)   angles_deg])
                  
          % Erotima 4
          disp('Power Flow in Transmission Lines [i -> j]')     
          disp('     from     | to | Pij (MW) | Qij (MVAr)')    
          disp([Lines_matrix(:,1) Lines_matrix(:,2) PowerFlow_lines_i_to_j(:,1) PowerFlow_lines_i_to_j(:,2)])

          disp('Power Flow in Transmission Lines [j -> i]')     
          disp('     from     | to | Pji (MW) | Qji (MVAr)')    
          disp([Lines_matrix(:,2) Lines_matrix(:,1) PowerFlow_lines_j_to_i(:,1) PowerFlow_lines_j_to_i(:,2)])

          % Erotima 5
          disp( 'Losses in each Transmission Line' )
          disp('     from     | to | Pji (MW) | Qji (MVAr)') 
          disp([Lines_matrix(:,1) Lines_matrix(:,2) (PowerFlow_lines_i_to_j(:,1) + PowerFlow_lines_j_to_i(:,1)) (PowerFlow_lines_i_to_j(:,2) + PowerFlow_lines_j_to_i(:,2))])
          
          % Erotima 6 
          disp('Total Losses in Transmission Lines ')    
          disp('Plosses (MW) in lines:')          
          disp(P_losses_ij)

          disp('Qlosses (MVar) in lines:')
          disp(Q_losses_ij)
     
          % % Erotima 6 - 2os tropos
          % disp('Prepei na einai ta idia me prin')
          % disp('Power Losses in System ')    
          % disp(' P_losses (MW) :')
          % P_losses = sum(Generators_matrix(:,2)) - sum(Bus_matrix_v_plus_1(:,3));           % sum P_G - sum P_L in (MW)
          % disp(P_losses)

          % disp('Q_losses (MVAr) :')          
          % Q_losses = sum(Generators_matrix(:,3)) - sum(Bus_matrix_v_plus_1(:,4));
          % disp(Q_losses)

          % Erotima 7           
          disp('For Slack Bus: S_1 = S_G1 - S_L1  (MVA)') 

          S_G1 = Generators_matrix(1,2) + 1i * Generators_matrix(1,3);     % complex

          S_L1 = Bus_matrix_v_plus_1(1, 3) + 1i * Bus_matrix_v_plus_1(1, 4); % complex
          %S_1 = S_G_SlackBus - S_L1;
          S_1 = S_G1 - S_L1;
          disp(S_1);
          return;
     else
         %disp(Delta_V_max)
         %disp('Still in the loop')
         iteration = iteration + 1;
         %disp(iteration);
     end
end      % exit for loop v

% since this code is out of the loop, v = v_max
disp('Maximum number of iterations exceeded.');
return;