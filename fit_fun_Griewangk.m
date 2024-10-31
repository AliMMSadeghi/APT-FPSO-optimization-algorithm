function y=fit_fun_Griewangk(x)

%--------------------------------------------------
%                       Griewangk method 4  [-600 600]
N=length(x(1,:));
nPt=length(x(:,1));
ii=sqrt(ones(nPt,1)*[1:1:N]);
y=1+sum(x.^2,2)/4000-prod(cos(x./ii),2);