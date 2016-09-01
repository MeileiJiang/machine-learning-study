function QQPlotComp(var1,var2, paramstruct)
%This function is used to create Q-Qplot for two vectors (var1, var2) compared with standard normal  
% Inputs:
%   data        - either n x 1 column vector of 1-d data
%                     or vector of bincts, when imptyp = -1
%
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
%    ndo                    number of plotted quantitles
%                               if the number of elements of min(var1,
%                               var2) < 1000, then use min (var1, var2)
%                               otherwise plot 1000 quantitles
%
%    dolcolor1         vector var1 overlay color
%                               string (any of 'r', 'g', 'k', etc.) for that single color
%                               default is 'b' (blue)
%                               1  Matlab 7 color default
%                               2  time series version (ordered spectrum of colors)
%                               nx3 color matrix:  a color label for each data point
%
%    dolcolor2         vector var2 overlay color
%                               string (any of 'r', 'g', 'k', etc.) for that single color
%                               default is 'g'
%                               1  Matlab 7 color default
%                               2  time series version (ordered spectrum of colors)
%                               nx3 color matrix:  a color label for each data point
%
%    dolmarkerstr1     Can be either a single string with symbol to use for marker,
%                                   e.g. 'o', '.'  default is '+' for var1
%                                   Or a character array (n x 1), of these symbols,
%                                   One for each data vector, created using:  strvcat
%
%    dolmarkerstr2     Can be either a single string with symbol to use for marker,
%                                   e.g. 'o', '.'  default is '*' for var2
%                                   Or a character array (n x 1), of these symbols,
%                                   One for each data vector, created using:  strvcat
%
%    dotsize          5  (default)  used for dot sizes
%                     
%    xlabelstr        string with x axis label
%                                    (only has effect when plot is made here)
%                           default is 'Theoretical Quantiles'
%
%    ylabelstr        string with y axis label
%                                    (only has effect when plot is made here)
%                           default is 'Quantiles of Input Sample'
%
%    titlestr         string with title (only has effect when plot is made here)
%                           default is empty string, '', for no title
%
%    linewidth       width of the 45 degree line indicating the standard
%                           normal distribution
%                           default is 2
%
%    linecolor        string with color for the 45 degree line
%                                    (only has effect when plot is made here)
%                           default is 'r' (red)
%                           use the empty string, '', for standard Matlab colors

if nargin == 1;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Error from QQPlotComp.m: !!!!!!!!!!') ;
    disp('!!!   Only one vector input  !!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Terminate Q-Q Plot Comparison  !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    return ;
end;

if length(var1)  ~= length(var2) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from QQPlotComp.m: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Two vectors have different numbers of elements !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
end ;


%  First set all parameters to defaults
if min(length(var1), length(var2)) < 1000;
    ndo = min(length(var1), length(var2)) ;
else
    ndo = 1000 ;
end;

dolcolor1 = 'b' ;
dolmarkerstr1 = '+' ;
dolcolor2 = 'g' ;
dolmarkerstr2 = '*' ;
linewidth = 2 ;
linecolor = 'r' ; 
dotsize = 5 ;
titlestr = '';
xlabelstr = 'Theoretical Quantiles' ;
ylabelstr = 'Quantiles of Input Sample' ;

%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 2 ;   %  then paramstruct has been added

  if isfield(paramstruct,'ndo') ;    %  then change to input value
    ndatovlay = getfield(paramstruct,'ndo') ; 
  end ;

  if isfield(paramstruct,'dolcolor1') ;    %  then change to input value
    dolcolor = getfield(paramstruct,'dolcolor1') ; 
  end ;
  
   if isfield(paramstruct,'dolcolor2') ;    %  then change to input value
    dolcolor = getfield(paramstruct,'dolcolor2') ; 
  end ;

  if isfield(paramstruct,'dolmarkerstr1') ;    %  then change to input value
    dolmarkerstr = getfield(paramstruct,'dolmarkerstr1') ; 
  end ;
  
  if isfield(paramstruct,'dolmarkerstr2') ;    %  then change to input value
    dolmarkerstr = getfield(paramstruct,'dolmarkerstr2') ; 
  end ;

  if isfield(paramstruct,'dotsize') ;    %  then change to input value
    ibigdot = getfield(paramstruct,'dotsize') ; 
  end ;

  if isfield(paramstruct,'linewidth') ;    %  then change to input value
    linewidth = getfield(paramstruct,'linewidth') ; 
  end ;

  if isfield(paramstruct,'linecolor') ;    %  then change to input value
    linecolor = getfield(paramstruct,'linecolor') ; 
  end ;

  if isfield(paramstruct,'titlestr') ;    %  then change to input value
    titlestr = getfield(paramstruct,'titlestr') ; 
  end ;
  
  if isfield(paramstruct,'xlabelstr') ;    %  then change to input value
    xlabelstr = getfield(paramstruct,'xlabelstr') ; 
  end ;

  if isfield(paramstruct,'ylabelstr') ;    %  then change to input value
    ylabelstr = getfield(paramstruct,'ylabelstr') ; 
  end ;

end ;  %  of resetting of input parameters

if nargin == 1;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Error from QQPlotComp.m: !!!!!!!!!!') ;
    disp('!!!   Only one vector input  !!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Terminate Q-Q Plot Comparison  !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    return ;
end;

if length(var1)  ~= length(var2) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from QQPlotComp.m: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Two vectors have different numbers of elements !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
end ;

if ndo > min(length(var1), length(var2)) ; 
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from QQPlotComp.m: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp( ['!!!   Numbers of elements is less then' num2str(ndo) '!!'] ) ;
    disp('!!!   Reset ndo to default !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    ndo = min(length(var1), length(var2)) ;
end;



vari=var1;
final_vari=var2;


p= (0+1/(ndo-1):1/(ndo-1):1-1/(ndo-1))';
y= icdf('normal',p,0,1);

qd_vari= quantile(vari, p);
qd_finalvari= quantile(final_vari,p);

plot(y,qd_vari, 'Marker', dolmarkerstr1,... 
                           'MarkerSize', dotsize, ...
                           'Linewidth',  linewidth-1, ...
                           'color', dolcolor1) ;
hold on;

plot(y,qd_finalvari, 'Marker', dolmarkerstr2,... 
                                   'MarkerSize', dotsize,...
                                   'Linewidth',  linewidth-1, ...
                                   'color', dolcolor2) ;

% plot the 45 degree line for comparison
plot([min(y) max(y)],[min(y) max(y)], 'color', linecolor, ...
                                                                'Linewidth',  linewidth, ...
                                                                'LineStyle', '--') ;
                                                            
xlabel( xlabelstr) ; 
ylabel( ylabelstr) ;
title( titlestr) ;

hold off;

end

