function [final_vari text_k] = AutoTransPara(vari, istat)
% This function will go through the series of transformations and will
% choose one that has the smallest absolute skewness or Anderson-Darling
% Statistics
% Input:
%   - vari              data vector to be transformed 

%   - istat              Criterion for selecting the parameter of
%                           transformation funtion 
%                           1: Skewness 
%                           2. Andreson-Darling Test Statistics (default)


% Output:
%   - final_vari        -   data vector after tranformation
%   - text_k              -   transformatio information

% Assumes path can find personal functions:
% ADStat.m


% The first  input arguments is required; the other is set to the default
% value if no input
if nargin < 2;
    istat = 2;    
end;


% Find the transformation parameter minimizing 
%  1. absoulte skewness (istat=1)
%  2. Andreson Darling Test Statistics (istat=2)

% get the interquartile range of the values in the data vector
iqrange_vari = iqr(vari);

% use the interquartile range as the data range if not zero 
% otherwise use the full data range
if iqrange_vari == 0;
    range_vari = range(vari);
else
    range_vari = iqrange_vari;
end;

% display the information of the data vector before transformation
if istat == 1; 
    val_vari = skewness(vari);
    disp ( 'Transformation Criterion: Minimize Skewness' ) ;
    disp (['Skewness Before Transformation: ' num2str(val_vari)]) ;
else
    if istat == 2;
     val_vari = ADStat(vari) ;
     disp ( 'Transformation Criterion: Minimize Log (Andreson_Darling Test Statistic) (Standard Normal)' ) ;
     disp ( ['Log A-D Stat Before Transformation: ' num2str(log(val_vari))] ) ;
    end;
end;

final_vari = vari;  % initial the transformed data vector
final_val = val_vari;   % value of istat given current transformation parameter
text_k = 'no tranformation';    % tranformation information

val_full = [ ] ;
para_full = [ ] ;
lestr = [ ] ;
[sortvari, index] = sort(vari) ;

% choosing transformation with smallest absoulte skewness/andrewson darling statistics

for i=[-9:0.01:-0.01 0.01:0.01:9];
    
     %change the scale of beta
     beta = sign(i)*(exp(abs(i))-1);
     
     % reparameterization 
     alpha = abs(1/(beta));
     
     % decided by sign of beta       
     if beta>0;
         vark1 = log(vari-min(vari)+alpha*range_vari);
     else
         vark1 = -log(max(vari)-vari+alpha*range_vari);
     end;    
         
     % standardize the transformed data vector 
     vark1 = (vark1 - mean(vark1)) / std(vark1) ;
     
     % truncate the outliers in the tails
     [d,n] = size(vari) ;
     p95 = evinv(0.95,0,1) ;
     ak = (2*log(n*d))^(-0.5) ;
     bk = (2*log(n*d))^0.5 -(log(log(n*d))+log(4*pi))/((2*log(n*d))^0.5) ;
     TH = p95*ak+bk;
 
     vark1(vark1>TH) = repmat(TH, 1, sum(vark1>TH)) ;
     vark1(vark1<-TH) = repmat(-TH, 1, sum(vark1<-TH)) ;
     
     % Restandardize 
     vark1 = (vark1 - mean(vark1)) / std(vark1) ;
     
     if istat == 1 ;
         val_vark1 = skewness(vark1) ;
     elseif istat == 2;
         val_vark1 = ADStat(vark1) ;
     end;
    
     val_full = [val_full val_vark1] ;
     para_full = [para_full; alpha i beta] ;
    
     % replace the previous parameter if the current istat is smaller 
     if abs(val_vark1) < abs(final_val) ;
        final_vari = vark1 ;
        final_val = val_vark1 ;
        if istat == 1;
            if skewness(vari) > 0 ;
               text_k =  [' transformation Parameter Beta =' num2str(beta)];
            else
               text_k =  [' transformation Parameter Beta =' num2str(-beta)];
            end;
        else
               text_k =  [' transformation Parameter Beta =' num2str(beta)];
        end;
     end; 

end;


% display the statistics of the data vector after transformation
disp( ['Selected Transformation: ' text_k] );

disp( ['Skewness After Transformation: ' num2str(skewness(final_vari))] );
           
disp( ['Log A-D Stat After Transformation: ' num2str(log(ADStat(final_vari)))] );
       

end

