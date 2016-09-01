disp('Running MATLAB script file AutoTranstest.m') ;
%
%    FOR DEVELOPMENT AND TESTING OF MATLAB FUNCTION AutoTrans,
%    General Purpose of automatical transformation of a given data vector

itest = 1 ;  % 1, 2, ..., 8

if itest == 1;      % check default setting 
    
    idat = exprnd(2, 20, 50) ;      % 20x50 data matrix
    
    [transformed_data, transformation] = AutoTrans(idat);
    
elseif itest ==2 ; % visual check of implementation
    
    idat = exprnd(2, 2, 50) ;      % 20x50 data matrix
    
    paramstruct = struct ( 'iplot', [1 0 0] ) ;
    [transformed_data, transformation] = AutoTrans(idat, paramstruct);
    
    paramstruct = struct ( 'iplot', [0 1 0] ) ;
    [transformed_data, transformation] = AutoTrans(idat, paramstruct);
    
    paramstruct = struct ( 'iplot', [0 0 1] ) ;
    [transformed_data, transformation] = AutoTrans(idat, paramstruct);
    
elseif itest == 3 ; % visual check of implementation of minimizing skewness
    
     idat = exprnd(3, 2, 100) ;      % 2x100 data matrix
    
     paramstruct = struct ( 'istat', 1) ;
     [transformed_data, transformation] = AutoTrans(idat, paramstruct);
    
elseif itest == 4;  % sample size more than 1000 (QQPlotcomp will plot 1000 quantiles)
    
    idat = exprnd (2, 5, 2000);         % 5x2000 data matrix
    
    paramstruct = struct ( 'iplot', [0 1 0]) ;
    [transformed_data, transformation] = AutoTrans(idat, paramstruct);
    
elseif itest == 5 ;  % check feature names
    
    idat = exprnd (2, 2, 200);         % 5x200 data matrix
    
    testname = strvcat('Testname1', 'Testname2') ;
    paramstruct = struct ( 'iplot', [1 1 1], ...
                                             'FeatureNames', testname) ;
    [transformed_data, transformation] = AutoTrans(idat, paramstruct);
    
elseif itest == 6; % check unmatched feature names 
    
     idat = exprnd (2, 1, 200);         % 5x200 data matrix
    
    testname = strvcat('Testname1', 'Testname2') ;
    paramstruct = struct ( 'iplot', [1 1 1], ...
                                             'FeatureNames', testname) ;
    [transformed_data, transformation] = AutoTrans(idat, paramstruct);
    
elseif itest == 7; % normal distributed data vector 
    
    idat = normrnd (1, 2, 1, 200);         % 1x200 data matrix
    
    paramstruct = struct ( 'iplot', [1 1 1]) ;
    [transformed_data, transformation] = AutoTrans(idat, paramstruct);
    
elseif itest == 8;  % data vector with binary values (  test error message )
    
     idat = binornd (1, 0.1, 1, 200)*10 + 5;       % 1x200 data matrix
     
     [transformed_data, transformation] = AutoTrans(idat);

end ;    %  of itest if-block



