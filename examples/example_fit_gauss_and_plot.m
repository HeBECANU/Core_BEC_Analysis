% voigt fit play




%%

x_data=linspace(0,1,1e3);
cof_names={'sigma','mu','amp','offset'};
gauss_fun=@(b,x) b(3)*exp(-(1/2)*((x-b(2))./b(1)).^2)+b(4);
beta0=[0.1,0.4,2,0.1];
y_data=gauss_fun(beta0,x_data);
% add a touch of noise
y_data=y_data+randn(size(y_data))*1e-1;

predictor=x_data;
response=y_data;


opt = statset('TolFun',1e-10,'TolX',1e-10,...
    'MaxIter',1e4,... %1e4
    'UseParallel',1);
% weights=ones(size(predictor));
% %'Weights',weights

fitobj=fitnlm(predictor,response,gauss_fun,beta0,...
    'options',opt,...
    'CoefficientNames',cof_names);

%
xplotvalues=linspace(min(predictor),max(predictor),1e4);
xplotvalues=col_vec(xplotvalues);
size(xplotvalues)
amp_pred=fitobj.predict(xplotvalues); %'Alpha',1-erf(1/sqrt(2)),'Prediction','observation'
size(amp_pred)

size(xplotvalues)
[amp_pred,ci]=predict(fitobj,xplotvalues,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation'); %'Alpha',1-erf(1/sqrt(2)),'Prediction','observation'
size(amp_pred)

shaded_ci_lines=true;
color_shaded=[0.9,1,1];

figure(1)
clf

hold on
if shaded_ci_lines
    patch([xplotvalues', fliplr(xplotvalues')], [ci(:,1)', fliplr(ci(:,2)')], color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
else
    plot(xplotvalues,ci(:,1),'-','LineWidth',1.5)
    plot(xplotvalues,ci(:,2),'-','LineWidth',1.5)
end  

plot(predictor,response,'k.')

plot(xplotvalues,amp_pred,'-','LineWidth',1.0)

hold off



