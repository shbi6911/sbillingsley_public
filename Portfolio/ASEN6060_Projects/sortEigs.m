function [eig_trim,eig_sorted,indices,triv_ind] = sortEigs(eig_values,vrb)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     3-18-2025
%
%This function sorts an nx6 matrix of eigenvalues of the monodromy matrix for a
%family of periodic orbits in the CR3BP.  It trims the trivial eigenvalues 
% and makes sure the others are properly ordered. The order is maximum
% value first, then its reciprocal, then maximum of the remaining two, then 
% its reciprocal.  
%
% INPUTS:   eig_values  an nx6 matrix of eigenvalues of monodromy
%                       matrices from a periodic orbit family
%           vrb         optional verbosity flag
%
% OUTPUTS:  eig_trim    an nx4 matrix of eigenvalues with the trivial
%                       eigenvalues removed
%           eig_sorted  an nx4 sorted matrix of eigenvalues, without the
%                       trivial values, in reciprocal pairs
%           indices     an nx4 matrix of indices, each row of which
%                       corresponds to the reordering performed on the rows 
%                       of eig_trim, i.e. eig_trim(n,indices(n,:)) =
%                       eig_sorted(n,:)
%           triv_ind    an nx2 matrix of indices corresponding to trivial
%                       eigenvalues of the original matrix i.e.
%                       eig_values(n, indices(n,:)) = [] results in
%                       eig_sorted(n,:)
    arguments
        eig_values (:,6) {mustBeNumeric}
        vrb logical =0
    end
    sz = size(eig_values); sz = sz(1);
    eig_trim = zeros(sz,4);     %preallocate
    indices = zeros(sz,4);
    eig_sorted = zeros(sz,4);
    triv_ind = zeros(sz,2);

    for i = 1:sz                %remove trivials
        diffs = real(eig_values(i,:)) - 1;
        [~,I] = mink(abs(diffs),2);
        temp = eig_values(i,:); temp(I) = [];
        eig_trim(i,:) = temp;
        triv_ind(i,:) = I;
    end

    for i = 1:sz
        temp = eig_trim(i,:);
        [~,index] = max(temp,[],"ComparisonMethod",'real');  %find maximum
        indices(i,1) = index;
        diffs = temp - 1/(eig_trim(i,indices(i,1)));
        [~,index] = min(diffs,[],"ComparisonMethod",'abs'); %find reciprocal
        indices(i,2) = index;
        temp(indices(i,1:2)) = NaN;
        [~,index] = max(temp,[],"ComparisonMethod",'real');  %find max remaining
        indices(i,3) = index;
        diffs = temp - 1/(eig_trim(i,indices(i,3)));
         [~,index] = min(diffs,[],"ComparisonMethod",'abs');    %find reciprocal
        indices(i,4) = index;
        eig_sorted(i,:) = eig_trim(i,indices(i,:));
    end %i = 1:length(eig_sorted)
end %function