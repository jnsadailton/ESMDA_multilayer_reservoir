function [q_vec] = get_q(t,rf,rw,rs,ks,k,qinj,h,muo,muw,kro,krw,phi,ct,alphat,layer)
    size_rf=length(rf);
    if size_rf == 1
        q_vec = qinj*ones(length(t),1);
    else
        for i=1:length(rs)
            lambchapos(i)=ks(i)*kro(1)/muo;
            lambchapor(i)=k(i)*kro(1)/muo;
            lambchapws(i)=ks(i)*krw(2)/muw;
            lambchapwr(i)=k(i)*krw(2)/muw;
            etaws(i)=ks(i)*krw(2)*alphat/(muw*phi*ct);
            etawr(i)=k(i)*krw(2)*alphat/(muw*phi*ct);
            etaos(i)=ks(i)*kro(1)*alphat/(muo*phi*ct);
            etaor(i)=k(i)*kro(1)*alphat/(muo*phi*ct);
        end
        for i=1:length(t)
                for j = 1:size_rf
                    rf_vec(j) = rf{j}(i);
                end
            q_vec(i,1)=valuate_q(t(i),rf_vec,rw,rs,lambchapos,lambchapws,lambchapor,lambchapwr,etaws,etaos,etawr,etaor,qinj,h,layer);
        end        
    end

end

