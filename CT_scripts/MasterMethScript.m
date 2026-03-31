for formulation = ["R"]
    for CypScore = [0, 0.5, 1, 1.5, 2]
        for RF = [0.25, 0.5, 0.75, 1]
            for BW = [50, 60, 70, 80, 90]
                ModelMain_Loewe_CT2(formulation, CypScore, RF, BW)
            end
        end
    end
end
