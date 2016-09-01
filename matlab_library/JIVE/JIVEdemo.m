% script: JIVE demo

% Toy data example 
load('Dataexample/ToyData.mat');
% orgainze data matrices into cell struct
datablock{1} = X;
datablock{2} = Y; 
% plot scree plot 
JIVEPreVisualQF(datablock);
% try rank 
rank{1} = [1 2 10];
rank{2} = [2 3];
JIVEPreVisualQF(datablock,rank);
% try rank 
rank{1} = 2;
rank{2} = 2:12;
JIVEPreVisualQF(datablock, rank); 

% select rank as vecr = [2, 3]
vecr = [2, 3];
outstruct = JIVEMainQF(datablock, vecr);

% output matrices
paramstruct = struct('ioutput', [1, 1, 1, 1, 1, 1, 1, 1, 1]);
vecr = [2, 3];
outstruct = JIVEMainQF(datablock, vecr, paramstruct);

% visualize matrices
joint = outstruct.joint{1};
individual = outstruct.individual{1};
JIVEdecompVisualQF(joint, individual)