%make table of IC50s etc

parameters =    ["IKr_methR"; 
                "IKr_methS";
                "IKr_metab"; 
                "INaL_methR";
                "INaL_methS";
                "INaL_metab"; 
                "ICaL_methR"; 
                "ICaL_methS";
                "ICaL_metab"];

        IC50s = [7e-6;
                7e-6;
                7e-6;
                31.8e-6;
                31.8e-6;
                31.8e-6;
                37.4e-6;
                37.4e-6;
                37.4e-6];

            h = [0.82;
                0.82;
                0.82
                1.37;
                1.37;
                1.37;
                1.67;
                1.67;
                1.67];
DrugPars = table(IC50s,h,RowNames=parameters);
save('DrugPars.mat', 'DrugPars');
    


