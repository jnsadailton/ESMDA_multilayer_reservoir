function [alpha] = get_alpha(Na,Ne,Cd,d_i)
    if Na == 1
        alpha = 1;
    else
        D=(d_i - mean(d_i,2))/(sqrt(Ne-1));
        Gi_d=sqrt(Cd)\D;
        N=min(length(d_i),Ne);
        [~,S,~] = svd(Gi_d);
        diagS=diag(S(:,1:N));
        
        lamb=mean(diagS);
        alpha_1=max(lamb^2,Na);
    
        f=@(x) 1 - alpha_1*(1-1/x)- x^(-Na);
    
        beta = bissec(f,0.00001,1);

        for i=1:Na
            alpha(i) = beta.^(i-1)*alpha_1;
        end
    end
    
    function [sol] = bissec(f, a, b)
        err=1;
        cont=0;
        fa  = f(a);
        while err>1e-4 
            sol = (a+b)/2;
            fsol = f(sol);
            if fsol==0 || (b-a)/2 < 1e-4
                break
            end
            if fa*fsol>0
                a=sol;
                fa=fsol;
            else
                b=sol;
            end
        end
    end
end

