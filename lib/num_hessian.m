function H=num_hessian(fun,x_list,delta)


%here I find the numerical hessian of the the passed scalar function potential
%it requires(1+2*3+4*3) 19 evaluations of the potential
%better than the 9*4=36 for naive simple implementation
%for the xx,yy,zz derivatives only need 1+2*3=7
%evaluations each (instead of 4*3=12)
%df/dxdx=f(x+2h)-2f(x)+f(x-2h)/h^2
%delta=delta*(0.5+rand(1)); %noise up the step to prevent false minima

num_pts=size(x_list,1);
dim=size(x_list,2);
fx=fun(x_list);
H=zeros(dim,dim,num_pts);
for i=1:dim
    for j=1:dim
        if i~=j && j>i
            %next we calculate the dj term at plus i
            dj_plus_i=(fun(x_list...
                                    +repmat([zeros(1,j-1),delta,zeros(1,dim-j)],num_pts,1)...
                                    +repmat([zeros(1,i-1),delta,zeros(1,dim-i)],num_pts,1))...
                        -fun(x_list...
                                    -repmat([zeros(1,j-1),delta,zeros(1,dim-j)],num_pts,1)...
                                    +repmat([zeros(1,i-1),delta,zeros(1,dim-i)],num_pts,1)))...
                        /(2*delta);
            
            %next we calculate the dj term at minus i
            dj_minus_i=(fun(x_list...
                                    +repmat([zeros(1,j-1),delta,zeros(1,dim-j)],num_pts,1)...
                                    -repmat([zeros(1,i-1),delta,zeros(1,dim-i)],num_pts,1))...
                        -fun(x_list...
                                    -repmat([zeros(1,j-1),delta,zeros(1,dim-j)],num_pts,1)...
                                    -repmat([zeros(1,i-1),delta,zeros(1,dim-i)],num_pts,1)))...
                        /(2*delta);
            H(i,j,:)=(dj_plus_i- dj_minus_i)/(2*delta);
        elseif j<i
            H(i,j,:)=H(j,i,:);  
        elseif j==i
            H(i,j,:)=(fun(x_list+repmat([zeros(1,j-1),2*delta,zeros(1,dim-j)],num_pts,1))...
            -2*fx...
            +fun(x_list-repmat([zeros(1,j-1),2*delta,zeros(1,dim-j)],num_pts,1)))/(4*delta^2);
        end
    end
end



end