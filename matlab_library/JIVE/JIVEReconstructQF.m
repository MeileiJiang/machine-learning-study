function outstruct = JIVEReconstructQF(datablock, threshold, dataname, row_joint, ioutput)
%  function for JIVE Matrix re-construction
% Inputs:
%   datablock        - cells of data matrices {datablock 1, ..., datablock k}
%                    - Each data matrix is a d x n matrix that each row is
%   dataname         - a cell of strings: name of each data block; default
%                    - name is {'datablock1', ..., 'datablockk'}
%   row_joint        - orthonormal basis of estimated joint row space 
%                      output of 'JIVEJointSelectQF.m'
%    ioutput         - 0-1 indicator vector of output's structure 
%                      [cns, cns-loading, bss-joint, bss-loading-joint, matrix-joint,... 
%                       bss-indiv, bss-loading-indiv, matrix-indiv, matrix-res] 

%    Copyright (c) Qing Feng Jan Hannig & J. S. Marron 2016

k = length(datablock);
rjoint = size(row_joint, 1); % proposed joint rank

joint = {};
jointcns = {};
jointbss = {};
jointloading = {};
jointprojload = {};
individual = {};
indivdualcns = {};
indivdualloading = {};
individualbss = {};
res = {};
rankI = [];

% check and re-adjust the joint rank
rjointr = [];
for ib = 1:k;
    
    load_joint = datablock{ib}/row_joint;
    joint{ib} = load_joint*row_joint;
    jointprojload{ib} = normc(load_joint);

    %remove components below threshold
    rj = length(find(sqrt(sum(load_joint.^2,1))>threshold(ib)));
    rjointr = [rjointr; rj];

    % check whether any singular value of joint components droping below the threshold
    % if so, throw a note
    if rj ~= rjoint;  
        disp(['Note: Joint rank of ' dataname{ib} ' = ' num2str(rj) ' is not equal to Proposed Joint rank' ]);
    end;
end 

rjoint = min(rjointr);
row_joint = row_joint(1:rjoint,:);
disp(['Final Joint rank:' num2str(rjoint)]);

for ib = 1:k;
    
    load_joint = datablock{ib}/row_joint;
    joint{ib} = load_joint*row_joint;
    jointprojload{ib} = normc(load_joint);

    %remove components below threshold
    rj = length(find(sqrt(sum(load_joint.^2,1))>threshold(ib)));

    [t1,t2,t3] = svds(joint{ib},rj);
    joint{ib} = t1*t2*t3';
    jointloading{ib} = t1;
    jointbss{ib} = t2*t3';
    jointcns{ib} = t3';

    %  Individual reconstruction 
    % orthogonal basis od null space of joint
    p = null(joint{ib})'; 
    proj_joint_row = p'*p;

    indiv = datablock{ib}*proj_joint_row;
    s_indiv = svd(indiv);

    rI = length(find(s_indiv>threshold(ib)));
    [t1,t2,t3] = svds(indiv,rI);
    individual{ib} = t1*t2*t3';
    indivdualcns{ib} = t3';
    indivdualloading{ib} = t1;
    individualbss{ib} = t2*t3';
    row_indiv = t3;
    rankI(ib) = rI;
    disp(['Final individual ' dataname{ib} ' rank: ' num2str(rI)]);

    %residual
    row = [orth(joint{ib}')';orth(individual{ib}')'];
    p = null(row)'; %orthogonal basis od null space of joint
    proj = p'*p;
    res{ib} = datablock{ib}*proj;
end;

% return needed results based on ioutput
if ioutput(1) == 1; % output common normalized score
    outstruct.jointcns = jointcns;
else 
    outstruct.jointcns = [];
end

if ioutput(2) == 1; % output projection loadings of common normalized score
    outstruct.jointprojload = jointprojload;
else 
    outstruct.jointprojload = {};
end

if ioutput(3) == 1; % output block specific scores of each joint
    outstruct.jointbss = jointbss;
else 
    outstruct.jointbss = {};
end

if ioutput(4) == 1; % output the loading matrix of each joint block specific score
    outstruct.jointloading = jointloading;
else 
    outstruct.jointloading = {};
end

if ioutput(5) == 1; % output joint matrices
    outstruct.joint = joint;
else 
    outstruct.joint = {};
end

if ioutput(6) == 1; % output block specific scores of each individual
    outstruct.individualbss = individualbss;
else 
    outstruct.individualbss = {};
end

if ioutput(7) == 1; % output the loading matrix of each individual block specific score
    outstruct.indivdualloading = indivdualloading;
else 
    outstruct.indivdualloading = {};
end

if ioutput(8) == 1; % output individual matrices
    outstruct.individual = individual;
else 
    outstruct.individual = {};
end

if ioutput(9) == 1; % output residual matrices
    outstruct.res = res;
else 
    outstruct.res = {};
end


