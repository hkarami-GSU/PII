function composite_val=R0s(params0)

%beta0/gamma
composite_val=(params0(:,4)*params0(:,1) + (1 - params0(:,4))*params0(:,2)) / params0(:,5);




