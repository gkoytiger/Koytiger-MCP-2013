function [cf gof] = fit_Kd(Xaxis,Yaxis, Kd_threshold)
%function fits Kd to series of concentrations and fluorescence intensities,
%using a robust fit


%fo_ = fitoptions('method','NonlinearLeastSquares','Robust','On','Lower',[0 0],'Upper',[Inf  Kd_threshold]); 
fo_ = fitoptions('method','NonlinearLeastSquares','Robust','On','Lower',[0 -1],'Upper',[Inf  Kd_threshold]); 

%st_ = [max(Yaxis) 1 ]; %statpoint for Kd is 1uM and Fmax is Maximum value of Fluorescence
st_ = [max(Yaxis) Kd_threshold];

set(fo_,'Startpoint',st_);
ft_ = fittype('Fmax * X / (Kd + X)',...
    'dependent',{'y'},'independent',{'X'},...
    'coefficients',{'Fmax', 'Kd'});

% Fit this model using new data
[cf gof] = fit(Xaxis,Yaxis,ft_,fo_); %returns fit object and goodness of fit object
%beta = coeffvalues(cf); %returns the fitted coefficient values (Fmax, Kd)
%Fmax = beta(1);
%Kd = beta(2);
%R2 = gof.rsquare;