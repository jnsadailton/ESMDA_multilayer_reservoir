function [q] = valuate_q(t,rf,rw,rs,lambchapos,lambchapws,lambchapor,lambchapwr,etaws,etaos,etawr,etaor,qinj,h,layer)
    function[val]=pl1(u)
        A=defineM(u,rf,rw,rs,lambchapos,lambchapws,lambchapor,lambchapwr,etaws,etaos,etawr,etaor,h);
        b=zeros(length(A),1);        
        b(1)=qinj/u;
        x=A\b;
        A11=x(5*(layer-1)+1);
        B11=x(5*(layer-1)+2);
        if rs(layer)==0
            val=-A11*besselk(1,rw.*sqrt(u./etawr(layer))).*sqrt(u./etawr(layer))+B11*besseli(1,rw.*sqrt(u./etawr(layer))).*sqrt(u./etawr(layer));
        else
            val=-A11*besselk(1,rw.*sqrt(u./etaws(layer))).*sqrt(u./etaws(layer))+B11*besseli(1,rw.*sqrt(u./etaws(layer))).*sqrt(u./etaws(layer));
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
    q=s1.*log(2)./t;
end