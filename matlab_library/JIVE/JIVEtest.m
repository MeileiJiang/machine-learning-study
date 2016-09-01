disp('Running MATLAB script file JIVEtest.m') ;
%
%    FOR DEVELOPMENT AND TESTING OF JIVE MATLAB FUNCTIONS

itest = 24 ;     %  1,...,24,101,...,105


close all ;

load('Dataexample/ToyData.mat')
%Y = Y(1:500,:) ;


if itest == 1 ; 

  disp('Test Main Program, no ranks input') ;
  JIVEMainQF({X,Y})

elseif itest == 2 ;

  disp('Test Main Program, small ranks input, 2 & 1') ;
  outstruct = JIVEMainQF({X,Y},[2 1]) 

elseif  itest == 3  |  ...
        itest == 4  |  ...
        itest == 5  |  ...
        itest == 6  |  ...
        itest == 7  |  ...
        itest == 8  |  ...
        itest == 9  |  ...
        itest == 22  |  ...
        itest == 23  |  ...
        itest == 24  ;

  if itest == 3 ;
    r1 = 1 ;
    r2 = 2 ;
  elseif itest == 4 ;
    r1 = 2 ;
    r2 = 2 ;
  elseif itest == 5 ;
    r1 = 2 ;
    r2 = 3 ;
  elseif itest == 6 ;
    r1 = 3 ;
    r2 = 3 ;
  elseif itest == 7 ;
    r1 = 4 ;
    r2 = 5 ;
  elseif itest == 8 ;
    r1 = 10 ;
    r2 = 5 ;
  elseif itest == 9 ;
    r1 = 20 ;
    r2 = 30 ;
  elseif itest == 22 ;
    r1 = 30 ;
    r2 = 20 ;
  elseif itest == 23 ;
    r1 = 30 ;
    r2 = 5 ;
  elseif itest == 24 ;
    r1 = 99 ;
    r2 = 99 ;
  end ;
  disp(['Test Main Program, input ranks: ' num2str(r1) ' & ' num2str(r2)]) ;
  paramstruct = struct('iplot',[0 1], ...
                       'ioutput',[0 0 0 0 1 0 0 1 0]) ;
  outstruct = JIVEMainQF({X,Y},[r1; r2],paramstruct)
  Xj = outstruct.joint{1} ;
  Xi = outstruct.individual{1} ;
  Yj = outstruct.joint{2} ;
  Yi = outstruct.individual{2} ;
  JIVEdecompVisualQF(Xj,Xi,'X') ;
  JIVEdecompVisualQF(Yj,Yi,'Y') ;

  figure(2) ;
  vax = axis ;
  text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
       vax(3) + 0.6 * (vax(4) - vax(3)), ...
       ['Input ranks:  r1 = ' num2str(r1) ',  r2 = ' num2str(r2)], ...
       'FontSize',15) ;
  text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
       vax(3) + 0.5 * (vax(4) - vax(3)), ...
       ['Final Joint Rank = ' num2str(rank(Xj))], ...
       'FontSize',15) ;
  text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
       vax(3) + 0.4 * (vax(4) - vax(3)), ...
       ['Final Individual X Rank = ' num2str(rank(Xi))], ...
       'FontSize',15) ;
  text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
       vax(3) + 0.3 * (vax(4) - vax(3)), ...
       ['Final Individual Y Rank = ' num2str(rank(Yi))], ...
       'FontSize',15) ;
  orient landscape ;
  pstr = ['JIVEtestIT' num2str(itest) '-r' num2str(r1) '-r' num2str(r2)] ;
  print('-dpsc2',pstr) ;
  

elseif itest == 10 ;

  disp('Test Main Program, ranks 2,3, mean subtraction 0 0') ;
  paramstruct = struct('imean',[0 0],'ioutput',[0 0 0 0 1 0 0 1 0]) ;
  outstruct = JIVEMainQF({X,Y},[2; 3],paramstruct)
  Xj = outstruct.joint{1} ;
  Xi = outstruct.individual{1} ;
  Yj = outstruct.joint{2} ;
  Yi = outstruct.individual{2} ;
  JIVEdecompVisualQF(Xj,Xi,'X') ;
  JIVEdecompVisualQF(Yj,Yi,'Y') ;

elseif itest == 11 ;

  disp('Test Main Program, ranks 2,3, mean subtraction 1 1') ;
  paramstruct = struct('imean',[1; 1],'ioutput',[0 0 0 0 1 0 0 1 0]) ;
  outstruct = JIVEMainQF({X,Y},[2; 3],paramstruct)
  Xj = outstruct.joint{1} ;
  Xi = outstruct.individual{1} ;
  Yj = outstruct.joint{2} ;
  Yi = outstruct.individual{2} ;
  JIVEdecompVisualQF(Xj,Xi,'X') ;
  JIVEdecompVisualQF(Yj,Yi,'Y') ;

elseif itest == 12 ;

  disp('Test Main Program, ranks 2,3, mean subtraction 2 3') ;
  paramstruct = struct('imean',[2; 3],'ioutput',[0 0 0 0 1 0 0 1 0]) ;
  outstruct = JIVEMainQF({X,Y},[2; 3],paramstruct)
  Xj = outstruct.joint{1} ;
  Xi = outstruct.individual{1} ;
  Yj = outstruct.joint{2} ;
  Yi = outstruct.individual{2} ;
  JIVEdecompVisualQF(Xj,Xi,'X') ;
  JIVEdecompVisualQF(Yj,Yi,'Y') ;

elseif itest == 13 ;

  disp('Test Main Program, ranks 2,3, mean subtraction 2 2') ;
  paramstruct = struct('imean',[2 2],'ioutput',[0 0 0 0 1 0 0 1 0]) ;
  outstruct = JIVEMainQF({X,Y},[2; 3],paramstruct)
  Xj = outstruct.joint{1} ;
  Xi = outstruct.individual{1} ;
  Yj = outstruct.joint{2} ;
  Yi = outstruct.individual{2} ;
  JIVEdecompVisualQF(Xj,Xi,'X') ;
  JIVEdecompVisualQF(Yj,Yi,'Y') ;

elseif itest == 14 ;

  disp('Test Main Program, ranks 2,3, iplot 0 0') ;
  paramstruct = struct('iplot',[0 0]) ;
  JIVEMainQF({X,Y},[2; 3],paramstruct)

elseif itest == 15 ;

  disp('Test Main Program, ranks 2,1, iplot 1 1') ;
  paramstruct = struct('iplot',[1 1]) ;
  JIVEMainQF({X,Y},[2; 1],paramstruct)

elseif itest == 16 ;

  disp('Test Main Program, ranks 3,5, iplot 1 1') ;
  paramstruct = struct('iplot',[1 1]) ;
  JIVEMainQF({X,Y},[3; 5],paramstruct)

elseif itest == 17 ;

  disp('Test Main Program, ranks 3,5, iplot 0 1, boundp = 65 85') ;
  paramstruct = struct('iplot',[0 1], ...
                       'boundp',[65 85]) ;
  JIVEMainQF({X,Y},[3; 5],paramstruct)

elseif itest == 18 ;

  disp('Test Main Program, ranks 3,5, iplot 0 1, boundp = 20, threp = 20') ;
  paramstruct = struct('iplot',[0 1], ...
                       'boundp',20, ...
                       'threp',20) ;
  JIVEMainQF({X,Y},[3; 5],paramstruct)

elseif itest == 19 ;

  disp('Test Main Program, ranks 3,5, iplot 0 1, boundp = 20, threp = 80') ;
  paramstruct = struct('iplot',[0 1], ...
                       'boundp',20, ...
                       'threp',80) ;
  JIVEMainQF({X,Y},[3; 5],paramstruct)

elseif itest == 20 ;

  disp('Test Main Program, ranks 3,5, iplot 1 1') ;
  paramstruct = struct('iplot',[1 1], ...
                       'dataname', ...
                       {{'Name X' 'Name Y'}}) ;
  JIVEMainQF({X,Y},[3; 5],paramstruct)

elseif itest == 21 ;

  disp('Test Main Program, ranks 3,5') ;
  paramstruct = struct('ioutput',[1 0 0 0 1 1 0 1 0]) ;
  oustruct = JIVEMainQF({X,Y},[3; 5],paramstruct)
  disp('Check inner products of CNS outputs') ;
  outstruct.jointcns{1} * outstruct.jointcns{2}'
  disp('Check that is rotation matrix') ;
  det(outstruct.jointcns{1} * outstruct.jointcns{2}')

elseif itest == 22 ;

  disp('Test Main Program, big ranks input, 30 & 20') ;
  paramstruct = struct('ioutput',[0 0 0 0 1 0 0 1 0]) ;
  outstruct = JIVEMainQF({X,Y},[30; 20],paramstruct)
  Xj = outstruct.joint{1} ;
  Xi = outstruct.individual{1} ;
  Yj = outstruct.joint{2} ;
  Yi = outstruct.individual{2} ;
  JIVEdecompVisualQF(Xj,Xi,'X') ;
  JIVEdecompVisualQF(Yj,Yi,'Y') ;

elseif itest == 23 ;

  disp('Test Main Program, big ranks input, 30 & 5') ;
  paramstruct = struct('ioutput',[0 0 0 0 1 0 0 1 0]) ;
  outstruct = JIVEMainQF({X,Y},[30; 5],paramstruct)
  Xj = outstruct.joint{1} ;
  Xi = outstruct.individual{1} ;
  Yj = outstruct.joint{2} ;
  Yi = outstruct.individual{2} ;
  JIVEdecompVisualQF(Xj,Xi,'X') ;
  JIVEdecompVisualQF(Yj,Yi,'Y') ;

elseif itest == 24 ;

  disp('Test Main Program, big ranks input, 99 & 99') ;
  paramstruct = struct('ioutput',[0 0 0 0 1 0 0 1 0]) ;
  outstruct = JIVEMainQF({X,Y},[99; 99],paramstruct)
  Xj = outstruct.joint{1} ;
  Xi = outstruct.individual{1} ;
  Yj = outstruct.joint{2} ;
  Yi = outstruct.individual{2} ;
  JIVEdecompVisualQF(Xj,Xi,'X') ;
  JIVEdecompVisualQF(Yj,Yi,'Y') ;

elseif itest == 101 ;

  disp('Test Initial SVD Code') ;
  JIVEPreVisualQF({X,Y})

elseif itest == 102 ;

  disp('Test Initial SVD Code') ;
  JIVEPreVisualQF({X,Y},{[1 2], [2 3]})

elseif itest == 103 ;

  disp('Test Initial SVD Code') ;
  JIVEPreVisualQF({X,Y},{[1 2 5 10], [5 15 25 35 45 55 65 75 85]})

elseif itest == 104 ;

  disp('Test Initial SVD Code, numcompshow 20') ;
  JIVEPreVisualQF({X,Y},{[1 2], [2 3]},20)

elseif itest == 105 ;

  disp('Test Initial SVD Code, numcompshow 200') ;
  JIVEPreVisualQF({X,Y},{[1 2], [2 3]},200)


end ;


