function [M] = defineM(u,rf,rw,rs,lambchapos,lambchapws,lambchapor,lambchapwr,etaws,etaos,etawr,etaor,h)
%% Using all equations to define the M Matrix
    tam = length(rs);
    number_equations_acoupling_between_regions = zeros(tam,1);
    ideal_size = 5*tam;
    for i = 1:tam
        if rs(i) == 0
            number_equations_acoupling_between_regions(i)=0;
            etaws(i)=etawr(i);
            lambchapws(i)=lambchapwr(i);
            etaos(i)=etaor(i);
            lambchapos(i)=lambchapor(i);
            rs(i)=rw+eps;
        else
            number_equations_acoupling_between_regions(i)=2;
        end
    end
    M=[];
    M_cell = cell(1,tam);
    for i = 1:tam
        M_cell{i} = zeros(ideal_size,ideal_size/tam);
        if i==1 && tam~=1
            eq2=[besselk(0,rw.*sqrt(u./etaws(i))) besseli(0,rw.*sqrt(u./etaws(i)))];
            M_cell{i}(ideal_size/tam+1,1:length(eq2))=eq2;
        elseif i==tam && tam~=1
            eq2=[-besselk(0,rw.*sqrt(u./etaws(i))) -besseli(0,rw.*sqrt(u./etaws(i)))];
            M_cell{i}(ideal_size/tam*(tam-1)+1,1:length(eq2))=eq2;
        else
            if tam~=1
                eq21=[-besselk(0,rw.*sqrt(u./etaws(i))) -besseli(0,rw.*sqrt(u./etaws(i)))];
                eq22=[besselk(0,rw.*sqrt(u./etaws(i))) besseli(0,rw.*sqrt(u./etaws(i)))];
                M_cell{i}(ideal_size/tam*(i-1)+1,1:length(eq21))=eq21;
                M_cell{i}(ideal_size/tam*(i)+1,1:length(eq22))=eq22;
            end
        end
        eq1=[rw*besselk(1,rw.*sqrt(u./etaws(i))).*sqrt(u./etaws(i))*h(i)*lambchapws(i) -rw*besseli(1,rw.*sqrt(u./etaws(i))).*sqrt(u./etaws(i))*h(i)*lambchapws(i)];
        M_cell{i}(1,1:length(eq1))=eq1;

        if rf(i) < rs(i)
            eq_rf=[besselk(0,rf(i).*sqrt(u./etaws(i))) besseli(0,rf(i).*sqrt(u./etaws(i))) -besselk(0,rf(i).*sqrt(u./etaos(i))) -besseli(0,rf(i).*sqrt(u./etaos(i)));
                        besselk(1,rf(i).*sqrt(u./etaws(i))) -besseli(1,rf(i).*sqrt(u./etaws(i))) -lambchapos(i)/lambchapws(i).*sqrt(etaws(i)/etaos(i)).*besselk(1,rf(i).*sqrt(u./etaos(i))) lambchapos(i)/lambchapws(i).*sqrt(etaws(i)/etaos(i)).*besseli(1,rf(i).*sqrt(u./etaos(i)))];
            eq_rs=[0 0 besselk(0,rs(i).*sqrt(u./etaos(i))) besseli(0,rs(i).*sqrt(u./etaos(i))) -besselk(0,rs(i).*sqrt(u./etaor(i)));
                        0 0 besselk(1,rs(i).*sqrt(u./etaos(i))) -besseli(1,rs(i).*sqrt(u./etaos(i))) -lambchapor(i)/lambchapos(i)*sqrt(etaos(i)/etaor(i))*besselk(1,rs(i).*sqrt(u./etaor(i)))];  
            if i~=1
                M_cell{i}(2+ideal_size/tam*(i-1):3+ideal_size/tam*(i-1),1:length(eq_rf))=eq_rf;
                M_cell{i}(4+ideal_size/tam*(i-1):5+ideal_size/tam*(i-1),1:length(eq_rs))=eq_rs;
            else
                M_cell{i}(2:3,1:length(eq_rf))=eq_rf;
                M_cell{i}(4:5,1:length(eq_rs))=eq_rs;
            end
        else
            eq_rs=[besselk(0,rs(i).*sqrt(u./etaws(i))) besseli(0,rs(i).*sqrt(u./etaws(i))) -besselk(0,rs(i).*sqrt(u./etawr(i))) -besseli(0,rs(i).*sqrt(u./etawr(i))) ;
                    besselk(1,rs(i).*sqrt(u./etaws(i))) -besseli(1,rs(i).*sqrt(u./etaws(i))) -lambchapwr(i)/lambchapws(i).*sqrt(etaws(i)/etawr(i)).*besselk(1,rs(i).*sqrt(u./etawr(i))) lambchapwr(i)/lambchapws(i).*sqrt(etaws(i)/etawr(i)).*besseli(1,rs(i).*sqrt(u./etawr(i)))];
            eq_rf=[0 0 besselk(0,rf(i).*sqrt(u./etawr(i))) besseli(0,rf(i).*sqrt(u./etawr(i))) -besselk(0,rf(i).*sqrt(u./etaor(i)));
                    0 0 besselk(1,rf(i).*sqrt(u./etawr(i))) -besseli(1,rf(i).*sqrt(u./etawr(i))) -lambchapor(i)/lambchapwr(i)*sqrt(etawr(i)/etaor(i))*besselk(1,rf(i).*sqrt(u./etaor(i)))];
                if i~=1
                    M_cell{i}(2+ideal_size/tam*(i-1):3+ideal_size/tam*(i-1),1:length(eq_rs))=eq_rs;
                    M_cell{i}(4+ideal_size/tam*(i-1):5+ideal_size/tam*(i-1),1:length(eq_rf))=eq_rf;
                else
                    M_cell{i}(2:3,1:length(eq_rs))=eq_rs;
                    M_cell{i}(4:5,1:length(eq_rf))=eq_rf;
                end
        end
        M=[M M_cell{i}];
    end
end

