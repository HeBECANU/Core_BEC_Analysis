function [sigma,out_st]=prop_error_cov(f,varlist,vals,errs,cov)
% propagate the error in the output of a function given errors in the input
% include the effects of variable covariance


% THIS FUNCTION COULD DO WITH FURTHER TESTS

% readings
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty
% https://www.physi.uni-heidelberg.de/~menzemer/Stat0708/statistik_vorlesung_6a.pdf


% based on code by  Brad Ridder 2007. Feel free to use this under the BSD guidelines. If
%you wish to add to this program, just leave my name and add yours to it.

if nargin<5
    cov=[];
end

out_st=[];
fsym=str2sym(f);
varlist_sym=cellfun(@str2sym,varlist);

% find the partial derivative of the expression y=f(x1,...xi) with respect to each x
iimax = numel(varlist);
pdy_pdxi = vpa(ones(1,iimax));
for ii = 1:iimax
    pdy_pdxi(ii) = diff(fsym,varlist_sym(ii),1);
end

normal_var=sum((subs(pdy_pdxi,varlist_sym,vals).^2).*(errs.^2));

if ~isempty(cov) && ~isdiag(cov)
    pd_prod_matrix= repmat(pdy_pdxi,[iimax,1]);
    pd_prod_matrix=pd_prod_matrix.*transpose(pd_prod_matrix);
    % set the diagnoal to zero
    pd_prod_matrix = pd_prod_matrix - diag(diag(pd_prod_matrix));
    % we can then find the conponent of the variance from the covariance terms
    cov_var=pd_prod_matrix.*cov;
    cov_var=subs(cov_var,varlist_sym,vals);
    cov_var=sum(cov_var(:));
else
    cov_var=0;
end

error1 =sqrt(normal_var+cov_var);
sigma = double(error1);
value=subs(fsym,varlist,vals);
value=double(value);

out_st.var_cov=double(cov_var);
out_st.var_normal=double(normal_var);
out_st.std=sigma;
out_st.std_wout_cov=double(sqrt(normal_var));
out_st.value=value;

     
end
     