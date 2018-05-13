%%%%%%%%%%%%%%%%%%%%
%%%% Calculate the initial n and Rs for the fit of the error equations 

 n0=1.5;
 RS0=1E-8;
 RSH0=500000;
 Initial values
b0=[n0, abs(RS0), abs(RSH0)]; 
%%%%%%%%%%%%%%%%%%%%


%%%%%%%% parameters for the fit
% options.Algorithm = 'levenberg-marquardt';
% options = optimset('DiffMinChange',1e0,'DiffMaxChange',1);
options = optimset('Display','none',... %change 'none' to 'iter' to see details
'TolX',1e-18,...
'TolFun',1e-18,...
'Algorithm','levenberg-marquardt');

% these the are data point from V=0 to V=Voc, +1
Is=I(ijsc:ivoc+1,1);
Vs=V(ijsc:ivoc+1,1);

% err=@(c0,Vs) -1*(Vs/c0(2)-c0(3)*(c0(2)*c0(4)+c0(2)*c0(5)+Vs)/(c0(2)*(c0(3)+c0(2)))+(c0(1)*k*T)/(c0(2))...
%     *lambertw((c0(2)*c0(5)*c0(3))/((c0(3)+c0(2))*(c0(1)*k*T))*exp(c0(3)*(c0(2)*c0(4)+c0(2)*c0(5)+Vs)/((c0(1)*k*T)*(c0(2)+c0(3))))))-Is;

errl=@(b0,Vs) -1*(Vs/b0(2)-b0(3)*(b0(2)*(Jsc+(b0(2)*Jsc-Voc)/b0(3))/(1-exp((b0(2)*Jsc-Voc)/(b0(1)*k*T)))+b0(2)*Voc/b0(3)+Vs)/(b0(2)*(b0(2)+b0(3)))+...
    (b0(1)*k*T)/(b0(2))*lambertw((b0(2))/(b0(1)*k*T)*((Jsc-(Voc)/(b0(2)+b0(3)))*exp(-Voc/(b0(1)*k*T)))/(1-exp((b0(2)*Jsc-Voc)/(b0(1)*k*T)))...
    *exp(b0(3)*(b0(2)*(Jsc+(b0(2)*Jsc-Voc)/b0(3))/(1-exp((b0(2)*Jsc-Voc)/(b0(1)*k*T)))+b0(2)*Voc/b0(3)+Vs)/(b0(1)*k*T*(b0(2)+b0(3))))))-Is;





%%%%%%%%%%%%%%%%%%%%
%%% Perform fit
% bl=lsqnonlin(@(b0) errl(b0,Vs),b0,[],[],options);

% With easily adjustable constraints
bl=lsqnonlin(@(b0) errl(b0,Vs),b0,[0,0,0],[2,1E10,1E12],options);