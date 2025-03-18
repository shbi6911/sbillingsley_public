function [pos_stab_mani,neg_stab_mani,pos_un_mani,neg_un_mani] = compManifold(state_matrix,N,d,vrb)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     3-18-2025
%
% this function takes in a matrix of 6-element state vectors corresponding
% to a corrected periodic orbit within the CR3BP.  It assumes each orbit
% has a stable and unstable manifold associated with it. It chooses N
% points along the orbit, finds their associated eigenvectors based on the
% monodromy matrix of the completed orbit, and takes a step d along the
% stable and unstable eigenvector at each point.  the function then outputs
% a matrix of these initial states for later integration.
%
% INPUTS    state_matrix    nx42 matrix of 42-element state vectors for a
%                           periodic orbit with State Transition Matrix
%           N               number of points along the orbit
%                           the function will use approximately this many
%                           points to compute manifolds
%           d               step distance along the eigenvector in
%                           NONDIMENSIONAL UNITS
%           vrb             optional verbosity flag
%
% OUTPUTS   pos_stab_mani   matrix of stable manifold initial conditions in
%                           a positive step direction
%           pos_un_mani     matrix of unstable manifold initial conditions
%                           in a positive step direction
%           neg_stab_mani   matrix of stable manifold initial conditions in
%                           a negative step direction
%           neg_un_mani     matrix of unstable manifold initial conditions 
%                           in a negative step direction
    arguments
        state_matrix (:,42) {mustBeNumeric}
        N   double
        d   double
        vrb logical =0
    end
    step = round(length(state_matrix)/N);
    indices = step:step:length(state_matrix);   %choose points
    %preallocate
    pos_stab_mani = zeros(length(indices),6);
    neg_stab_mani = zeros(length(indices),6);
    pos_un_mani = zeros(length(indices),6);
    neg_un_mani = zeros(length(indices),6);
    monodromy = reshape(state_matrix(end,7:42),6,6);    %extract monodromy matrix
    [mon_vec,mon_val] = eig(monodromy);             %get eigenvalues/vectors
    mon_val = diag(mon_val);
    [~,eig_sorted,sort_ind,triv_ind] = sortEigs(mon_val');   %sort eigenvalues
    mon_vec(:,triv_ind) = [];   %remove trivial eigenvectors
    mon_vec = mon_vec(:,sort_ind);     %sort eigenvectors
    %isolate stable and unstable eigenvectors of monodromy matrix
    unstab_vec = mon_vec(:,1);
    stab_vec = mon_vec(:,2);
    %loop through 
    for i = 1:length(indices)
        state = state_matrix(indices(i),1:6)'; %extract state variables
        phi_j = reshape(state_matrix(indices(i),7:42),6,6); %extract STM
        unstab_vec_pt = phi_j*unstab_vec;   %set eigenvectors for this point
        stab_vec_pt = phi_j*stab_vec;
        sf_unstab = 1/norm(unstab_vec_pt(1:3)); %scale eigenvectors by position
        sf_stab = 1/norm(stab_vec_pt(1:3));
        unstab_vec_pt = unstab_vec_pt*sf_unstab;
        stab_vec_pt = stab_vec_pt*sf_stab;
        pos_stab_mani(i,:) = (state + d*stab_vec_pt)';
        neg_stab_mani(i,:) = (state - d*stab_vec_pt)';
        pos_un_mani(i,:) = (state + d*unstab_vec_pt)';
        neg_un_mani(i,:) = (state - d*unstab_vec_pt)';
    end % for i = 1:length(indices)
end %function







