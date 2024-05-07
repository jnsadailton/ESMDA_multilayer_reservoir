function [deltapl] = p_derivative(t,deltaP)
     n= length(deltaP);
    deltapl(n)=0;
    deltapl(1)=0;
    for i=2:(n-1)
        deltapl(i)=(deltaP(i+1)-deltaP(i))/log(t(i+1)/t(i))*log(t(i)/t(i-1))/log(t(i+1)/t(i-1))+(deltaP(i)-deltaP(i-1))/log(t(i)/t(i-1))*log(t(i+1)/t(i))/log(t(i+1)/t(i-1));
    end
end

