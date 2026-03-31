%make table of IC50s etc. These are taken from Cliff's syncropatch data.

parameters =    ["IKr_methR"; 
                "IKr_methS";
                "IKr_metab"; 
                "INaL_methR";
                "INaL_methS";
                "INaL_metab"; 
                "ICaL_methR"; 
                "ICaL_methS";
                "ICaL_metab"];

        IC50s = [1.89e-6;
                2.09e-6;
                3.29e-5;
                4.5e-6;
                9.34e-6;
                4.13e-5;
                3.08e-5;
                2.07e-5;
                1.19e-4];

            h = [1; %0.7631;
                1; %0.8279;
                1; %0.7941;
                1; %0.988;
                1; %1.411;
                1; %1.562;
                1; %0.754;
                1; %0.6664;
                1;]; %0.9776];
DrugPars2 = table(IC50s,h,'RowNames',parameters);
save('DrugPars.mat', 'DrugPars2');
    


