clear all
clear classes
fun = @(x) x^2;
x0 = -1;
f0 = fun(x0);
df0 = (fun(x0+1e-3)-f0)/1e-3;
d = -df0;
lsObj = lineSearch(fun,x0,d,'golden',[],[],'directionSwitch',true,'c1',0.2,'c2',1e-9);
[xout,fval,stepSize,relStepSize,nFevalOut,nGradEvalOut,exitflag,message] = lsObj.solve