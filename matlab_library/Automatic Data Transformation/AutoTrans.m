function [transformed_data, transformation] = AutoTrans(mdata, paramstruct)
% This program will automatically transform each row of inputted data set
% to make the distrubution of each row vector closer to standard normal
% distribution
%
% The transformation function belongs to a new parametrization of shifted
% logarithm family, in which the transformation parameter is choosen to 
% minimize the skewness/Andreson- Darling Statistics of the transformed data vecotr
%

% Inputs:
%         mdata - d x n matrix of data, 
%                 the n columns: each column corresponds with each sample
%                 the d rows: each row represnts a feature vector

%   paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%                    Version for easy copying and modification:
%    paramstruct = struct('',, ...
%                         '',, ...
%                         '',) ;
%
%    fields                 values
%    istat               Criterion for selecting the parameter of
%                           transformation funtion 
%                           1: Skewness 
%                           2. Andreson-Darling Test Statistics (default)
%
%   iplot                 [0 0 0] (default)
%                            indicator vector of whether to plot each data vector before and
%                            after transformation
%                            [kernel density estimate plot, Q-Q plot, sizer
%                            plot]
%
%   FeatureNames       Feature names of each data vector
%                                     default setting: 'Feature1' 
%                                     use strvcat to write each feature
%                                     name into each row 
%

% Assumes path can find personal functions:
% AutoTransPara.m
% ADStat.m
% kdeSM.m
% sizerSM.m


[d, n] = size(mdata) ;
transformed_data = [] ;
transformation = {} ;

%  First set all parameters to defaults
istat = 2 ;
iplot = [0 0 0] ;
FeatureNames = [] ;
for i=1:d;
    FeatureNames = strvcat(FeatureNames, ['Feature' num2str(i)]) ;
end;



if nargin>1;
    
    if isfield(paramstruct, 'istat');
        istat = paramstruct.istat;
    end;
    
    if isfield(paramstruct, 'iplot');
        iplot = paramstruct.iplot;
    end;
    
    if isfield(paramstruct, 'FeatureNames');
        FeatureNames = paramstruct.FeatureNames;
    end;
    
    %when the number of feature names is unmatached with number of features
    %needed to perform transformation, use default feature names 
    if d ~= size(FeatureNames, 1);
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        disp('!!! Warning: Number of Feature Names!!!!!');
        disp('!!! Unmatched with Number of Features!!!!');
        disp('!!! Use Default Set for Feature Names!!!!!!!');
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        FeatureNames=[];
        for i=1:d;
            FeatureNames = strvcat(FeatureNames, ['Feature' num2str(i)]) ;
        end;
    end;

end ;

for i = 1 : size(mdata,1) ;
    
    vari = mdata(i, :) ;
    
    % Standardize
    vari = (vari - mean(vari)) / std(vari) ;
      

%   if number of unique values is greater then 2    
    if length(unique(vari))<=2;
        disp ('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        disp ('!!!   Error from AutoOpTransQF.m:  !!!') ;
        disp ('!!!   Invalid Data Input                          !!!') ;
        disp ('!!!   Standardize only                           !!!') ;
        disp ('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        final_vari = vari;
        return ;
    else
        
        [final_vari text_k] = AutoTransPara(vari, istat);
    
    end;
    
transformed_data = [transformed_data;  final_vari] ;
transformation{i, 1} = [ FeatureNames(i, :) ': ' text_k] ;

%% Visualization of the transformation effect

% kernel density estimate plot before and after transformation
% use defualt density estimate: Sheather Jones Plug In

if  iplot(1) ;
    figure ;
    subplot(1,2,1);
    paramstructkde = struct ( 'ibigdot',0, ...
                                                    'isubpopkde',0, ...
                                                    'idatovlay', 1, ...
                                                    'datovlaymin', 0.4,...
                                                    'datovlaymax', 0.6,...
                                                    'titlestr', [ FeatureNames(i, :) ': KDE Before Transformation']) ;
    kdeSM(vari', paramstructkde);

    
    subplot(1,2,2);
    paramstructkde = struct ( 'ibigdot',0, ...
                                                    'isubpopkde',0, ...
                                                    'idatovlay', 1, ...
                                                    'datovlaymin', 0.4,...
                                                    'datovlaymax', 0.6,...
                                                    'titlestr', [ FeatureNames(i, :) ': KDE After Transformation']) ;
    kdeSM(final_vari', paramstructkde);
    savestr=['KDE_' FeatureNames(i, :)];
    print('-dpsc',savestr) ;
end;

%QQplot
if iplot(2);
    figure;
    paramstructQQ = struct ('titlestr', ['QQ Plot Contrast: ' FeatureNames(i, :) ]) ;
    QQPlotComp(vari, final_vari, paramstructQQ);
    savestr=['QQPlot' FeatureNames(i, :)];
    print('-dpsc', savestr) ;
end;


%Sizer
if iplot(3);
    figure ; 
    paramstructsizer = struct ('famoltitle',    [FeatureNames(i, :) ': Sizer Before Transformation'], ...
                                                     'imovie', 0, ...
                                                     'savestr', ['SizerBfTrans' FeatureNames(i, :)]);
    sizerSM(vari', paramstructsizer);
    

    figure ;
    paramstructsizer = struct ('famoltitle',    [FeatureNames(i, :) ': Sizer After Transformation'], ...
                                                     'imovie', 0, ...
                                                     'savestr', ['SizerAfTrans' FeatureNames(i, :)]);
    sizerSM(final_vari', paramstructsizer);

end;



end;


end

