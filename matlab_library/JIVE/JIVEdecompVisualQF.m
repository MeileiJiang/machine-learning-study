function JIVEdecompVisualQF(joint, individual, dataname)
% heat-map view of JIVE decomposition 
% Input: estimates of joint, individual matrices 
%        dataname: name of the data 
% Output: graphics only 

%    Copyright (c) J. S. Marron, Jan Hannig & Qing Feng

[d, n] = size(joint);
[t1, t2] = size(individual);
if d~=t1 || n ~= t2;
    disp('Error Message: JIVEdecompVisualQF terminated as dimensions of joint and individual are not matched!')
    return;
end;

if nargin == 2;
    dataname = 'Data'; % set default value if no input
end

% set the color scheme 
b=0.99:-0.01:0;
bb=0:0.01:0.99;
blpha=[ones(length(b),1)*1 b' b'];bblpha=[bb' bb' ones(length(bb),1)*1 ];
map=[bblpha;[1 1 1]; blpha;];

vax = axisSM(joint+individual, (joint+individual)') ;
limup= max(max(vax), max(vax));
limd =  min(min(vax), min(vax));

figure;
ax = subplot(1, 3, 1);
imagesc(joint + individual);
colormap(map);
caxis([min(vax), max(vax)]);
set(gca, 'XTick', []) ; % remove X-axis
c=colorbar('southoutside', 'FontSize', 6); %add colorbar
%axpos = ax.Position;
axpos = get(ax,'Position') ;
    %  Modification of Qing's command to work with old Matlab
%cpos = c.Position;
cpos = get(c,'Position') ;
    %  Modification of Qing's command to work with old Matlab
cpos(2) = cpos(2) - 5*cpos(4);
cpos(4) = 0.5*cpos(4);
%c.Position = cpos;
set(c,'Position',cpos) ;
    %  Modification of Qing's command to work with old Matlab
%ax.Position = axpos;
set(ax,'Position',axpos) ;
    %  Modification of Qing's command to work with old Matlab
title([dataname ': JIVE Approx'], 'FontSize', 14) ;  % Add title for subplot


ax = subplot(1, 3, 2);
imagesc(joint);
colormap(map);
caxis([min(vax), max(vax)]);
set(gca, 'XTick', []) ; % remove X-axis
c=colorbar('southoutside', 'FontSize', 6); %add colorbar
%axpos = ax.Position;
axpos = get(ax,'Position') ;
    %  Modification of Qing's command to work with old Matlab
%cpos = c.Position;
cpos = get(c,'Position') ;
    %  Modification of Qing's command to work with old Matlab
cpos(2) = cpos(2) - 5*cpos(4);
cpos(4) = 0.5*cpos(4);
%c.Position = cpos;
set(c,'Position',cpos) ;
    %  Modification of Qing's command to work with old Matlab
%ax.Position = axpos;
set(ax,'Position',axpos) ;
    %  Modification of Qing's command to work with old Matlab
 % Add title for subplot
title('Joint Approx', 'FontSize', 14) ;
set(gca, 'YTick', []) ;

ax = subplot(1, 3, 3);
imagesc(individual);
colormap(map);
caxis([min(vax), max(vax)]);
set(gca, 'XTick', []) ; % remove X-axis
c=colorbar('southoutside', 'FontSize', 6); %add colorbar
%axpos = ax.Position;
axpos = get(ax,'Position') ;
    %  Modification of Qing's command to work with old Matlab
%cpos = c.Position;
cpos = get(c,'Position') ;
    %  Modification of Qing's command to work with old Matlab
cpos(2) = cpos(2) - 5*cpos(4);
cpos(4) = 0.5*cpos(4);
%c.Position = cpos;
set(c,'Position',cpos) ;
    %  Modification of Qing's command to work with old Matlab
%ax.Position = axpos;
set(ax,'Position',axpos) ;
    %  Modification of Qing's command to work with old Matlab
title('Individual Approx', 'FontSize', 14) ;  % Add title for subplot
set(gca, 'YTick', []) ;


print('-depsc', 'JIVEOutput');
print('-dpsc', 'JIVEOutput');
