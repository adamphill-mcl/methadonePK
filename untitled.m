% Inputs
A = 2;        % Concentration of drug A
B = 0;        % Concentration of drug B
EC50_A = 1;   % Half-maximal effective concentration for drug A
EC50_B = 1;   % Half-maximal effective concentration for drug B
h_A = 2;      % Hill slope for drug A
h_B = 0.5;      % Hill slope for drug B

% Calculate Loewe additivity term
% Each term is normalized to the respective EC50 and scaled by the Hill coefficient
Loewe_A = (A / EC50_A)^h_A; 
Loewe_B = (B / EC50_B)^h_B;

% Combine according to Loewe additivity
Combined_Loewe = Loewe_A + Loewe_B;

% The Loewe additivity relationship equals 1 when the combined effect is achieved
% In a dose-response context, we’d often normalize this by finding the point where Combined_Loewe = 1
% Here, we display the raw combined Loewe term
disp(['Combined Loewe term: ', num2str(Combined_Loewe)]);

combo_block = (Loewe_A + Loewe_B)/(Loewe_A + Loewe_B + 1)