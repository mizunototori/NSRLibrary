%import matlab.unittest.TestSuite
classdef CalcrecovTest < matlab.unittest.TestCase

	properties
		param
        result_fname=['./results/test_RunMU_NSR_Wl1.mat'];		
        %result_fname=['./results/test_RunKSVD.mat'];		

        TrueDic
        TrueCoef
        Dictionary
        CoefMatrix
        Aatoms  % True atom
        Eatoms  % Estimated atom      
	end

	methods (TestMethodSetup)
		function importLearnedDic(testCase)
			importeddata =  import_history(testCase.result_fname);
			testCase.TrueDic    = importeddata.TrueDic;
			testCase.TrueCoef   = importeddata.TrueCoef;
			testCase.Dictionary = importeddata.Dictionary;
			testCase.CoefMatrix = importeddata.CoefMatrix;
		end
	end


	methods (TestMethodTeardown)
		function createFigure(testCase)
            N = 10;
            K = 4;            
            for k = 0:K
                for i = 0:N-1
                    figure(k+1)
                    num = k*N + i+1;
                    subplot(N,1,i+1)
                    plot(testCase.Aatoms(:,num)); hold on; plot(testCase.Eatoms(:,num))
                end
            end
		end
	end

	methods (Test)
		function test_compare_atoms(testCase)
			[ratio,totalDistances, testCase.Aatoms, testCase.Eatoms] = syntheticexpt.calcrecovratio(testCase.TrueDic, testCase.Dictionary)            
        end
	end
end