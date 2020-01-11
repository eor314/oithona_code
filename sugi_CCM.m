function [ SugiCorr , predY , predX , origY , origX ] = sugi_CCM( xx , yy , tau , E  )

% Convergent Cross Mapping from 
% Sugihara et al., "Detecting Causality in Complex Ecosystems." Science
% 2012
% Inputs:
%   xx = first data series
%   yy = second data series
%   tau = time lag 
%   E = embedding dimension
% Outputs:
%   SugiCorr = matrix of corrrelation between the perdicted and actual data
%   predX = predited time series of X from CCM of Y
%   predY = predicted time series of Y from CCM of X
%   origY = orignal Y starting from T = 1 + (E - 1)*tau
%   origX = orignal X starting from T = 1 + (E - 1)*tau


L=length(xx);
T=1+(E-1)*tau;
Xm=zeros((L-T+1),E);
Ym=zeros((L-T+1),E);
SugiN=E+1;
N = L-T+1;

%% RECONTRUCTIONS OF ORIGINAL SYSTEMS

for t=1:(L-T+1)
    Xm(t,:)=xx((T+t-1):-tau:(T+t-1-(E-1)*tau));
    Ym(t,:)=yy((T+t-1):-tau:(T+t-1-(E-1)*tau));
end
%%

predX=zeros(N,1);
predY=zeros(N,1);

origY=yy(T:end);
origX=xx(T:end);

parfor j=1:N
    % neighborhood search 

    [n1,d1]=knnsearch(Xm,Xm(j,:),'k',E+2);
    [n2,d2]=knnsearch(Ym,Ym(j,:),'k',E+2);
    susY=origY(n1(2:end)); 
    susX=origX(n2(2:end));

    % CMM

    SugsusY=susY(1:SugiN);
    SugsusX=susX(1:SugiN);
    Sugid1=d1(:,2:SugiN+1);
    Sugid2=d2(:,2:SugiN+1);
    u1=exp(-Sugid1./(Sugid1(:,1)*ones(1,SugiN)));
    u2=exp(-Sugid2./(Sugid2(:,1)*ones(1,SugiN)));
    w1=u1./(sum(u1,2)*ones(1,SugiN));
    w2=u2./(sum(u2,2)*ones(1,SugiN));
    predY(j)= w1*SugsusY;
    predX(j)= w2*SugsusX;

end

SugiCorr1=corrcoef(origY,predY);
SugiCorr(2,1)=SugiCorr1(1,2);

SugiCorr2=corrcoef(origX,predX);
SugiCorr(1,1)=SugiCorr2(1,2);

end