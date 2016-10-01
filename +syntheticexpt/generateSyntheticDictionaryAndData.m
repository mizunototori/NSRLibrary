function [Dictionary, data, coefs] = gererateSyntheticDictionaryAndData(N, L, dim, K, SNRdB)
%% Returns dictionary, a data matrix and a coefficient matrix

%      		Data 		= 		Dictionary 		x 	Coefficient matrix
%      _____ ... _____      _________________       _____ ... ______ 
%     |               |    |				 |     |				| 
%     |               |    |				 |     |				| 
%     |               | =  |				 |  x  |				| 
%     |               |    |				 |     |				| 
%     |               |    |				 |     |				| 
%      _____ ... _____      _________________      |				|
%												   |				|
%												    _____ ... ______ 
%			
%		  (dim x N)			    (dim x K)				 (K x N)
%							 Over complete			 Coefficients has L nonzero elements
%
%		Paramaters:
%			size of matrices : 
%				N : scala
%					number of signals to generate 
%				dim : scala
%					dimension of each data 
%				K : scala
%					number of dictionary elements
%			sparsity : 
%				L : scala
%				number of elements in each linear combination. (number of non-zero elements)
%			Noise level :	
%				SNRdB : scala
%					level of noise be added. The range of level is [0, 80)	

	% Generate a dictionary
	Dictionary = rand(dim,K); 
	Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));

	% Generate a data and coefficients matrix
	[data,coefs] = CreateDataFromDictionary(Dictionary, N, L);

	% Add noise
	if (SNRdB==0) | (SNRdB == 80) 
	    return
	else
		% Generate noise matrix
	    noise = rand(size(data));

	    % Generate actual noise
	    actualNoise = calcNoiseFromSNR(SNRdB,data, noise);

	    % Add SNRdB noise 
	    SNR = calcSNR(data, data+actualNoise);
	    data =  data + actualNoise*SNR/SNRdB;   
	end
end

function [data, coefficients] = CreateDataFromDictionary(dictionary, N, L)
% Returns Data  and  Coefficients
%		Parameters : 
%			dictionary : matrix
%				a given dictionary 
%			N numElements : scala
%				size of column of coefficient
%			L numCoef : scala
%				number of nonzero elements of coefficient

	%	resolution = 0.0001; % TODO: Remove this if this doesn't be used
	
	%	number of dictionary elements
	K = size(dictionary,2);

	% Initialize by zero matrix ()
	coefficients = zeros(K,N);

	% Generate L x N matrix 
	max_coef = 1;
	coefs = rand(L,N)*max_coef; 

	% Insert non-zero elements to coefficient matrix
	coefficients(1:L,:) = coefs;	


	% Arrange non-zero elements randomly
	for i=1:N
	    coefficients(:,i) = coefficients(randperm(size(coefficients,1)),i);
	end

	% Generate data 
	data = dictionary*coefficients;
end

function  actualNoise = calcNoiseFromSNR(TargerSNR, signal, randomNoise)
%% Returns actual noise
%		Parameters : 
%			TargerSNR : scalar
%			signal : matrix (or vector)
%				given data
%			randomNoise : matrix (or vector)
%				given noise to add 

	% matrix to vector 
	signalRow = signal(:);
	randomNoiseRow = randomNoise(:);

	% Calculate actual noise ratio
	signal_power = sum(signalRow.^2); % scalar
	actualNoise = signal_power/(10^(TargerSNR/10)); % scalar
	noise = sum(randomNoiseRow.^2); % scaler
	ratio = actualNoise./noise; % scalar

	% Generate actual noise matrix
	actualNoise = randomNoiseRow.*repmat(sqrt(ratio),size(randomNoiseRow,1),1);
	actualNoise = reshape(actualNoise,size(randomNoise));
end


function SNR = calcSNR(origSignal, noisySignal)
%% Returns Signal Noise Ratio (SNR)
%		Parameters : 
%			origSignal : matrix (or vector)
%				original signal
%			noisySignal	: matrix (or vector)
%				signal added noise

	errorSignal = origSignal-noisySignal;
	signal_2 = sum(origSignal.^2);
	noise_2 = sum(errorSignal.^2);
	SNRValues = 10*log10(signal_2./noise_2);
	SNR = mean(SNRValues);
end
