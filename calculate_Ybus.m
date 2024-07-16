function [Y_busFinal, Z, B] = calculate_Ybus(Lines_matrix,Bus_matrix_v)
% Arguments : matrix of lines and matrix of busses
% Outputs : The Ybus of the given system

% Calculate Ybus
Z = Lines_matrix(:,3)+1i*Lines_matrix(:,4);
B = 1i*Lines_matrix(:,5);          % this is Ye
Y = 1./Z;                          % this is not Ye
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
               v = L1(n2);
               Ybus(n1,v) = -Y(n2);
          end
          
          if((L1(n2) == n1) || (L2(n2) == n1))
               Ybus(n1,n1) = sum(Ybus(n1,n1)) + Y(n2) + (B(n2)/2);
          end
     end
end
Y_busFinal = Ybus;
end