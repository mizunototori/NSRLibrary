function Dictionary = cleardictionary(Dictionary,CoefMatrix,Data)
	[r,c]=size(Dictionary);
	[n,N]=size(CoefMatrix);
	norms = sqrt(sum(Dictionary.^2,1));
	Dictionary = Dictionary./(ones(r,1)*norms);
	CoefMatrix = CoefMatrix.*(norms'*ones(1,N));
	% Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));
	T2 = 0.999;%0.99;
	T1 = 3;
	K=size(Dictionary,2);

	Er=sum((Data-Dictionary*CoefMatrix).^2,1); % remove identical atoms

	G=Dictionary'*Dictionary; G = G-diag(diag(G));
	for jj=1:1:K,
	    if max(G(jj,:))>T2 || length(find(abs(CoefMatrix(jj,:))>1e-7))<=T1 ,       
	        [val,pos]=max(Er);
	        Er(pos(1))=0;
	        Dictionary(:,jj)=Data(:,pos(1))/norm(Data(:,pos(1)));
	        Dic(:,jj)=Data(:,pos(1));
	        G=Dictionary'*Dictionary; G = G-diag(diag(G));
    end;
end;

%{
function Dictionary = I_clearDictionary(Dictionary,CoefMatrix,Data)
T2 = 0.999;
T1 = 3;
K=size(Dictionary,2);
Er=sum((Data-Dictionary*CoefMatrix).^2,1); % remove identical atoms
G=Dictionary'*Dictionary; G = G-diag(diag(G));
for jj=1:1:K,
    if max(G(jj,:))>T2 | length(find(abs(CoefMatrix(jj,:))>1e-7))<=T1 ,
        [val,pos]=max(Er);
        Er(pos(1))=0;
        Dictionary(:,jj)=Data(:,pos(1))/norm(Data(:,pos(1)));
        G=Dictionary'*Dictionary; G = G-diag(diag(G));
    end;
end;	
%}