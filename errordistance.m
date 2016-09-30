function error = errordistance(A,B)   
    error = norm((A-B),'fro');
    %error = (sum(sum((A-B).^2)))/sum(sum((A.^2))); % ISSUE: is it correct?
end