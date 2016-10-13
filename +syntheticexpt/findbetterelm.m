function [betterDictionaryElement,CoefMatrix,NewVectorAdded] = I_findBetterDictionaryElement(Data,Dictionary,j,CoefMatrix,numCoefUsed)
if (length(who('numCoefUsed'))==0)
    numCoefUsed = 1;
end
relevantDataIndices = find(CoefMatrix(j,:)); % the data indices that uses the j'th dictionary element.
if (length(relevantDataIndices)<1) %(length(relevantDataIndices)==0)
    ErrorMat = Data-Dictionary*CoefMatrix;
    ErrorNormVec = sum(ErrorMat.^2);
    [d,i] = max(ErrorNormVec);
    betterDictionaryElement = Data(:,i);%ErrorMat(:,i); %
    betterDictionaryElement = betterDictionaryElement./sqrt(betterDictionaryElement'*betterDictionaryElement);
    betterDictionaryElement = betterDictionaryElement.*sign(betterDictionaryElement(1));
    CoefMatrix(j,:) = 0;
    NewVectorAdded = 1;
    return;
end

NewVectorAdded = 0;
reduced_coeff = CoefMatrix(:, relevantDataIndices);
reduced_Data = Data (:, relevantDataIndices);
saveDebugDict = Dictionary(:,j);
saveDebugCoef = reduced_coeff(j,:);
%debug1 = sum(sum((reduced_Data - Dictionary*reduced_coeff).^2));
reduced_coeff (j, :) = 0;    % all but the j-th element

err_mat = reduced_Data - Dictionary * reduced_coeff;
[U S V flag] = svds((err_mat), 1);


% check for sign, flip U and V's sign if negative.
Idx_U = find(U<0);
Idx_V = find(V<0);

u1 = U; u1(Idx_U) = 0; v1 = V; v1(Idx_V) = 0;approx1 = norm(err_mat- u1*v1'*S);
u1 = zeros(size(U)); u1(Idx_U) = -U(Idx_U); v1 = zeros(size(V)); v1(Idx_V) = -V(Idx_V);approx2 = norm(err_mat- u1*v1'*S);
if (approx1<= approx2)
	betterDictionaryElement = U;
	betterDictionaryElement(Idx_U) = 0;
	coefs = V;
	coefs(Idx_V) = 0;
else
	betterDictionaryElement = zeros(size(U));
	betterDictionaryElement(Idx_U) = -U(Idx_U);
	coefs = zeros(size(V));
	coefs(Idx_V) = -V(Idx_V);
end

newAtomNorm = sqrt(betterDictionaryElement'*betterDictionaryElement);
betterDictionaryElement = betterDictionaryElement/newAtomNorm;
coefs = coefs * newAtomNorm;
% coefs(coefs<0) = 0;

newE = sum(sum(((reduced_Data - Dictionary(:,[1:j-1,j+1:end])*reduced_coeff([1:j-1,j+1:end],:))-betterDictionaryElement*coefs').^2));
oldE = sum(sum(((reduced_Data - Dictionary(:,[1:j-1,j+1:end])*reduced_coeff([1:j-1,j+1:end],:))-saveDebugDict*saveDebugCoef).^2));
if (newE>oldE)
	for iter = 1:30 % the number of iterations
		betterDictionaryElement = err_mat*coefs/(coefs'*coefs);
		betterDictionaryElement(betterDictionaryElement<0) = 0;
		coefs = err_mat'*betterDictionaryElement/(betterDictionaryElement'*betterDictionaryElement);
		coefs(coefs<0) = 0;
	end
	newAtomNorm = sqrt(betterDictionaryElement'*betterDictionaryElement);
	betterDictionaryElement = betterDictionaryElement/newAtomNorm;
	coefs = coefs * newAtomNorm;
	reduced_coeff(j,:) =  coefs;
end
reduced_coeff(j,:) =  coefs;
CoefMatrix (:, relevantDataIndices) =reduced_coeff;