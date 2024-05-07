function [pwf] = valuate_pwf(t,rf,rw,rs,lambchapos,lambchapws,lambchapor,lambchapwr,etaws,etaos,etawr,etaor,qinj,h)
    function[val]=pl1(u)
        A=defineM(u,rf,rw,rs,lambchapos,lambchapws,lambchapor,lambchapwr,etaws,etaos,etawr,etaor,h);
        b=zeros(length(A),1);        
        b(1)=qinj/u;
        x=A\b;
        A11=x(1);
        B11=x(2);
        if rs(1)==0
            val=A11*besselk(0,rw.*sqrt(u./etawr(1)))+B11*besseli(0,rw.*sqrt(u./etawr(1)));
        else
            val=A11*besselk(0,rw.*sqrt(u./etaws(1)))+B11*besseli(0,rw.*sqrt(u./etaws(1)));
        end
    end

    s1=0;
    N=12;
    for j=1:N
        s2=0;
        for k=floor((j+1)/2):min(j,N/2)
            s2=s2+(k^(1+N/2)*factorial(2*k)/(factorial(N/2-k)*(factorial(k))^2*factorial(j-k)*factorial(2*k-j)));  
        end
        vj=s2*(-1)^(j+N/2);
        s1=s1+vj*pl1(log(2)*j./t);
    end
    pwf=s1.*log(2)./t;
end