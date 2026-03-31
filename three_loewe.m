% Inputs for three drugs
A = 4;        % Concentration of drug A
B = 0;        % Concentration of drug B
C = 0;        % Concentration of drug C
EC50_A = 2;   % Half-maximal effective concentration for drug A
EC50_B = 2;   % Half-maximal effective concentration for drug B
EC50_C = 2;  % Half-maximal effective concentration for drug C
h_A = 1;      % Hill slope for drug A
h_B = 1;      % Hill slope for drug B
h_C = 1;      % Hill slope for drug C

% Calculate Loewe additivity terms
Loewe_A = (A / EC50_A)^h_A; 
Loewe_B = (B / EC50_B)^h_B;
Loewe_C = (C / EC50_C)^h_C;

% Combine according to Loewe additivity
Combined_Loewe = Loewe_A + Loewe_B + Loewe_C;

% Display the result
disp(['Combined Loewe term for three drugs: ', num2str(Combined_Loewe)]);

Combo_block = (Loewe_A + Loewe_B + Loewe_C) / (Loewe_A + Loewe_B + Loewe_C + 1)