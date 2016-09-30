function [ratio,totalDistances] = calcrecovratio(true_dic,est_dic)
% Returns ratio of recoverd atoms and distance between a true dictionary and an estimated dictionary
%    Parameters:
%        true_dic : matrix
%           true dictionary
%        est_dic  : matrix

    num_recoverd_atoms =0;
    totalDistances = 0;
    num_atoms = size(true_dic,2);    

    for i = 1:num_atoms;
        
        % Aatom: all atom of the true dictionary
        Aatom = true_dic(:,i); 

        % Compare distance between an atom of true dictionary and atoms of the estimated dictionary.
        distances =mean((est_dic-repmat(Aatom,1,num_atoms)).^2);    
        
        % Get min distance and the atom of estimated dictionary
        [minValue,min_index] = min(distances); 

        % Eatom: an atom of the estimated dictionary which is closse to true dictionary
        Eatom = est_dic(:,min_index);

        % Measure error of Atom and Eatom
        errorOfElement = 1-abs(Eatom'*Aatom); % Assumed normarized?

        % Integrate error distance
        totalDistances = totalDistances+errorOfElement;

        % Integrate recoverd atoms
        num_recoverd_atoms = num_recoverd_atoms+(errorOfElement<0.01);%0.01
    end
    ratio = 100*(num_recoverd_atoms/num_atoms);
end