clear
close all
clc
tic
rw=0.108;
n = 1; % n>1 if you want to compare diffents compiles
for k = 1:n
%% Cases
% % Case A


qinj=500;   %flow-rate
muo=5.1;    %oil viscosity
phi=0.32;   %porosity
kj=[1000 1000]; %permeability per layer
h=[10 15];  %thickness per layer
ksj=[500 100];  %skin zone permeability per layer
rsj=[0.3 0.5];  %skin zone radii

% % Case B
% qinj=200;
% phi=0.12;
% muo=1.0;
% kj=[600 600];
% h=[10 10];
% ksj=[240 0];
% rsj=[0.25 0];


% % Case C
% qinj = 400;
% phi = 0.15;
% muo = 1.5;
% kj = [1000 1200 800];
% h = [10 8 7];
% ksj = [250 400 300];
% rsj = [0.4 0.5 0.2];
% % rsj = [0.4 0.5 0.3];


%%Case D
% qinj=100;
% phi=0.25;
% muo=2.3;
% kj=[600 800];
% h=[15 10];
% ksj=[150 200];
% rsj=[0.2 0.2];

%% General parameters

muw=0.5; %Water viscosity
cr=8.0E-5;
cw=4.04E-5;
co=1.14E-4;
swi=0;
sor=0;
ct=cr+swi*cw+(1-sor)*co;

kro=[0.54 0];
krw=[0 0.17];


alphap=19.03;
alphat=0.0003484;
qinj=alphap*qinj;



%% Defining terms
%Definings lambchap and etas
for i=1:length(rsj)
    lambchapos(i)=ksj(i)*kro(1)/muo;
    lambchapor(i)=kj(i)*kro(1)/muo;
    lambchapws(i)=ksj(i)*krw(2)/muw;
    lambchapwr(i)=kj(i)*krw(2)/muw;
    etaws(i)=ksj(i)*krw(2)*alphat/(muw*phi*ct);
    etawr(i)=kj(i)*krw(2)*alphat/(muw*phi*ct);
    etaos(i)=ksj(i)*kro(1)*alphat/(muo*phi*ct);
    etaor(i)=kj(i)*kro(1)*alphat/(muo*phi*ct);
end

%time
t(1)=1e-4;
tFluxo=10;
j=1;

%Defining the water front radii on each layer
rf=cell(1,length(rsj));
for i = 1:length(rsj)
    while t(j) <= tFluxo
        rf{i}(j)=sqrt(rw*rw+((qinj/alphap)*t(j)/(24*phi*pi*h(i)*(1-swi-sor)))); %raio da frente de avanço 1 
        j=j+1;
        t(j)=t(j-1)*1.2;
    end
    j=1;
end
t=t(1:end-1);



%% Calculating pwf, defining ESMDA parameters
    n_layers = length(rsj);
    pwf = get_pwf(t,rf,rw,rsj,ksj,kj,qinj,h,muo,muw,kro,krw,phi,ct,alphat); 
    
    include_q = 1; %zero in case you don't want to include flow-rate
    if include_q
        root_flow = [2,16,32,48,63]; %flow-rate indices to consider
        q_i = [];
        for i = 1:(n_layers-1)
            dpdr = get_q(t,rf,rw,rsj,ksj,kj,qinj,h,muo,muw,kro,krw,phi,ct,alphat,i);
            q = -lambchapws(i)*h(i)*dpdr*rw;
            q_i = [q_i; q(root_flow)];
        end
    end
    root_rand = randn(length(pwf),1);
    pwf_no_error = pwf;
    pwf = pwf + 0.015*root_rand.*pwf; %include error
    dpwf = p_derivative(t,pwf);  %calculate dpwf

    p_normalize = 0;
    q_normalize = 0;
    q_norm = @(x) x/qinj;
    p_norm = @(x) x/max(pwf);
    
    if p_normalize
        maxpwf=max(pwf);
        pwf = p_norm(pwf);
    end
    
    if q_normalize && include_q
        q_i = q_norm(q_i);
    end
    
    rs_true = rsj;
    ks_true = ksj;
    kj_true = kj;
    
    Ne = 1000;
    Na = 4;
    
    d_obs = pwf;
    if n_layers ~=1 && include_q
        d_obs = [pwf;q_i];
    end
    
    Nd = length(d_obs);
    d_i = zeros(Nd, Ne);

    Cd = 5*1e-2*eye(Nd);
    Cd(1:64,:) = Cd(1:64,:)/5;  
    
    %% Limits of parameters
    cte_rsj = 0; % if is cte =1, if is not cte = 0
    cte_ksj = 0;
    cte_kj = 0;
    
    rj_min = rw;
    rj_max = 1.5;

    ksj_min = 10;
    ksj_max = 1000;

    kj_min = 50;
    kj_max = 5000; 

    A = 3;
    B = 3;
    if cte_rsj == 0
         for i=1:n_layers
             xmin=rj_min;
             xmax=rj_max;
             rsj1(i,:) = xmin + (xmax - xmin)*betarnd(A,B,1,Ne);
         end
        rsj = rsj1';
        rsj_ini=rsj;
    else
        rsj_ini=rsj;
        for j=1:n_layers
            rsj(j,1:Ne)=rsj_ini(i);
        end
        rsj=rsj';
    end
    rsj(rsj<=rw)=rw+eps;
    rsj(rsj>=2)=2;
    if cte_ksj == 0
        for i=1:n_layers      
             xmin=ksj_min;
             xmax=ksj_max;
             ksj1(i,:) = xmin + (xmax - xmin)*betarnd(A,B,1,Ne);


        end
        ksj= ksj1';
        ksj_ini=ksj;
    else
        ksj_ini=ksj;
        for j=1:n_layers
            ksj(j,1:Ne)= ksj_ini(j);
        end
        ksj=ksj';
    end
    
    if cte_kj == 0
        for i=1:n_layers
             xmin=kj_min;
             xmax=kj_max;
             kj1(i,:) = xmin + (xmax - xmin)*betarnd(A,B,1,Ne);
        end
        kj= kj1';
        
        kj_ini=kj;
    else
        kj_ini=kj;
        for j=1:n_layers
            kj(j,1:Ne)= kj_ini(j);
        end
        kj=kj';
    end
    
     %% Applying ESMDA
    qi_not_normalized = zeros(5*(n_layers-1),Ne);
    for i = 1:Ne
            p = get_pwf(t,rf,rw,rsj(i,:),ksj(i,:),kj(i,:),qinj,h,muo,muw,kro,krw,phi,ct,alphat);
            p_not_normalized(:,i) = p;
            if p_normalize
                p = p_norm(p);
            end
            if n_layers ~= 1 && include_q
                q_tot = [];
                for j = 1:(n_layers-1)
                    dpdr = get_q(t,rf,rw,rsj(i,:),ksj(i,:),kj(i,:),qinj,h,muo,muw,kro,krw,phi,ct,alphat,j);
                    q = -lambchapws(j)*h(j)*dpdr*rw;
                    qi_not_normalized((1+length(root_flow)*(j-1)):length(root_flow)*j,i) =  q(root_flow);                
                    if q_normalize
                        q = q_norm(q);   
                    end
                    q_tot = [q_tot; q(root_flow)];
                end
            else
                q_tot=[];
            end
            d_i(:, i) = [p;q_tot];
    end
    d_k = d_i;
    [alphal]=get_alpha(Na,Ne,Cd,d_i);    
    for ii = 1:Na
        if ii>1
            for i = 1:Ne
                p = get_pwf(t,rf,rw,rsj(i,:),ksj(i,:),kj(i,:),qinj,h,muo,muw,kro,krw,phi,ct,alphat);
                if p_normalize
                    p = p_norm(p);
                end
                if n_layers ~= 1 && include_q
                    q_tot = [];
                    for j = 1:(n_layers-1)
                        dpdr1 = get_q(t,rf,rw,rsj(i,:),ksj(i,:),kj(i,:),qinj,h,muo,muw,kro,krw,phi,ct,alphat,j);
                        q = -lambchapws(j)*h(j)*dpdr1*rw;
                        if q_normalize
                            q= q_norm(q);
                        end
                        q_tot = [q_tot; q(root_flow)];
                    end
                else
                    q_tot=[];
                end
                d_k(:, i) = [p;q_tot];
            end
        end
        
        tam = length(d_obs);
        d_ucjk = zeros(size(d_k));
    
        for i = 1:Ne
              d_ucjk(:,i) = d_obs + sqrt(alphal(ii))*sqrt(Cd)*randn(size(d_obs));
        end
        m=[];
        if cte_rsj == 0
            m = [rsj'];
        end
        if cte_ksj == 0 
            m = [m; log(ksj)'];
        end
        if cte_kj == 0
            m = [m; log(kj)'];
        end

        deltaDk = (d_k - mean(d_k,2))/sqrt(Ne-1);
        deltaMk = (m - mean(m,2))/sqrt(Ne-1);
        pseudoMk = pinv(deltaMk);
        Gkd = sqrt(Cd)\deltaDk;
        K = (Gkd')/(Gkd*Gkd'+alphal(ii)*eye(Nd));
        m = m + pseudoMk\(K*( sqrt(Cd)\(d_ucjk - d_k)  ));
    
        if cte_rsj == 0
            rsj = m(1:n_layers, :)';
            rsj(rsj<=rj_min)=rj_min+eps;
        end
    
        if cte_rsj == 0
            rsj = m(1:n_layers, :)';
            rsj(rsj<=rj_min)=rj_min+eps;
        end

        if cte_ksj == 0
            ksj = exp(m(n_layers+1:2*n_layers, :)');
            ksj(ksj>1000) = 1000;
        end
        if cte_kj == 0
            kj = exp(m(2*n_layers+1:3*n_layers, :)');
            kj(kj>9000) = 9000;
        end
    
    end
%final update of data
    for i = 1:Ne
        p = get_pwf(t,rf,rw,rsj(i,:),ksj(i,:),kj(i,:),qinj,h,muo,muw,kro,krw,phi,ct,alphat);
        if n_layers ~=1 && include_q
            q_tot = [];
            for j = 1:(n_layers-1)
                dpdr = get_q(t,rf,rw,rsj(i,:),ksj(i,:),kj(i,:),qinj,h,muo,muw,kro,krw,phi,ct,alphat,j);
                q = -lambchapws(j)*h(j)*dpdr*rw;
                q_tot = [q_tot; q(root_flow)];
            end
        else
            q_tot=[];
        end
        d_k(:, i) = [p;q_tot];    
    end
    rs{k}=rsj;
    perm{k}=kj;
    ks{k}=ksj;
end
somar=zeros(Ne,n_layers);
somak=somar;
somaks=somar;
for j = 1:n
    somar = somar+ rs{j};
    somak = somak + perm{j};
    somaks = somaks + ks{j};
end
rsj = somar/n;
kj = somak/n;
ksj = somaks/n;

%% Plots
sizepwf = length(pwf);

d_i(1:sizepwf,:)=p_not_normalized;

if p_normalize
    pwf = pwf*maxpwf;
end

%Pwf e dpwf
figure
% initial ensemble
loglog(t,d_i(1:sizepwf,1), 'Color', [.5 .5 .5])
hold on
% final ensemble
loglog(t, d_k(1:sizepwf,1), 'y','LineWidth',2)
% observed data
loglog(t,pwf,'o', 'MarkerSize', 5, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [0 0 0])
%Derivative
loglog(t,dpwf,'^','MarkerEdgeColor','g','MarkerFaceColor','g')
% Final ensemble mean
loglog(t,mean(d_k(1:sizepwf,:),2),'b','LineWidth',2)


loglog(t,d_i(1:sizepwf,2:Ne), 'Color', [.5 .5 .5])
loglog(t, d_k(1:sizepwf,2:Ne), 'y','LineWidth',2)

loglog(t,mean(d_k(1:sizepwf,:),2),'b','LineWidth',2)

loglog(t,pwf,'o', 'MarkerSize', 5, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [0 0 0])
loglog(t,dpwf,'^','MarkerEdgeColor','g','MarkerFaceColor','g')



legend('Initial Ensemble','Final Ensemble','Observed Data(\DeltaP)','Derivative','Mean','Location','southeast')


grid on
xlabel('time (h)')
ylabel('\DeltaP (kgf/cm^2)')



p = gcf;
nomepressao=sprintf('caso_%s_pressao_.eps',caso);
exportgraphics(p,nomepressao,'ContentType','vector')

%flow-rate graphic
if n_layers~=1 && include_q
    d_i(length(pwf)+1:Nd,:) = qi_not_normalized;
    if q_normalize
        q_i = q_i*qinj;
    end
    for j = 1:(n_layers - 1)
        figure
        % Initial Ensemble
        loglog(t(root_flow),(d_i(length(pwf)+1+length(root_flow)*(j-1):(Nd-length(root_flow)*(n_layers-1-j)),1))/alphap, 'Color', [.5 .5 .5])
        hold on
        % Final Ensemble
        loglog(t(root_flow), (d_k(length(pwf)+1+length(root_flow)*(j-1):(Nd-length(root_flow)*(n_layers-1-j)),1))/alphap, 'y','LineWidth',2)
        % observed data
        loglog(t(root_flow), (q_i(1+length(root_flow)*(j-1):length(root_flow)*j))/alphap,'o', 'MarkerSize', 5, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [0 0 0])
        
         % mean of final ensemble
        loglog(t(root_flow),mean((d_k(length(pwf)+1+length(root_flow)*(j-1):(Nd-length(root_flow)*(n_layers-1-j)),:))/alphap,2),'b','LineWidth',2)
 
        loglog(t(root_flow), (d_i(length(pwf)+1+length(root_flow)*(j-1):(Nd-length(root_flow)*(n_layers-1-j)),2:Ne))/alphap, 'Color', [.5 .5 .5])
        loglog(t(root_flow), (d_k(length(pwf)+1+length(root_flow)*(j-1):(Nd-length(root_flow)*(n_layers-1-j)),2:Ne))/alphap, 'y','LineWidth',2)     
        loglog(t(root_flow),mean((d_k(length(pwf)+1+length(root_flow)*(j-1):(Nd-length(root_flow)*(n_layers-1-j)),:))/alphap,2),'b','LineWidth',2)
 
        loglog(t(root_flow),(q_i(1+length(root_flow)*(j-1):length(root_flow)*j))/alphap,'o', 'MarkerSize', 5, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [0 0 0])        
        
        grid on
        xlabel('time (h)')
        ylabel('Q(m^3/d)')
        legend('Initial Ensemble','Final Ensemble','Observed Data(Q)','Mean')
        q = gcf;
        nameofplotq=sprintf('caso_%s_vazao_layer%i.eps',caso,j);
        exportgraphics(q,nameofplotq,'ContentType','vector')
    end
end

%Create histograms
for i=1:n_layers
    tam=3-(cte_kj+cte_ksj+cte_rsj);
    figure
    if cte_rsj == 0
        histogram(rsj_ini(:,i),(0.05:0.05:(max(rsj_ini(:,i))+0.1)),'FaceColor','g')
        hold on
        histogram(rsj(:,i),(0.05:0.05:(max(rsj(:,i))+0.1)),'FaceColor','b')
        xline(rs_true(i),'-',{'r_s true'},'LineWidth',2,'Color','r')
        legend('Initial Ensemble','Final Ensemble')
    end
    xlabel('Skin zone radius')
    ylabel('Frequency')
    rsj_fig = gcf;
        nameofplotrsj=sprintf('caso_%s_rs_layer%i.eps',caso,i);
        exportgraphics(rsj_fig,nameofplotrsj,'ContentType','vector')
    figure
    if cte_ksj == 0 
        histogram(ksj_ini(:,i),20,'FaceColor','g')
        hold on
        histogram(ksj(:,i),'FaceColor','b')
        if rs_true(i)==0
            ks_true(i)=kj_true(i);
        end
        xline(ks_true(i),'-',{'k_s true'},'LineWidth',2,'Color','r')
        legend('Initial Ensemble','Final Ensemble')
        set(gcf,'Position',[400 75 600 600])
    end
    xlabel('Skin zone Permeability')
    ylabel('Frequency')
    ksj_fig = gcf;
        nameofplotksj=sprintf('caso_%s_kslayer%i.eps',caso,i);
        exportgraphics(ksj_fig,nameofplotksj,'ContentType','vector')
    figure
    if cte_kj == 0
        histogram(kj_ini(:,i),20,'FaceColor','g')
        hold on
        histogram(kj(:,i),'FaceColor','b')
        xline(kj_true(i),'-',{'k true'},'LineWidth',2,'Color','r')
        legend('Initial Ensemble','Final Ensemble')
    end
    xlabel('Permeability')
    ylabel('Frequency')
    kj_fig = gcf;
    nameofplotkj=sprintf('caso_%s_klayer%i.eps',caso,i);
        exportgraphics(kj_fig,nameofplotkj,'ContentType','vector')
end

rsj_median = median(rsj);
ksj_median = median(ksj);
kj_median = median(kj);
pwf_median = get_pwf(t,rf,rw,rsj_median,ksj_median,kj_median,qinj,h,muo,muw,kro,krw,phi,ct,alphat);
pwf_mean = get_pwf(t,rf,rw,mean(rsj),mean(ksj),mean(kj),qinj,h,muo,muw,kro,krw,phi,ct,alphat);
dpwf_median = p_derivative(t,pwf_median);
dpwf_mean = p_derivative(t,pwf_mean);
figure
loglog(t,pwf_no_error,'o', 'MarkerSize', 5, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [0 0 0])
hold on
loglog(t,pwf_median,'o', 'MarkerSize', 5, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
hold on
loglog(t,pwf_mean,'o', 'MarkerSize', 5, 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 0])
hold on
loglog(t,p_derivative(t,pwf_no_error),'^','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0])
hold on
loglog(t,dpwf_median,'^','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
hold on
loglog(t,dpwf_mean,'^','MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1])

legend('Pwf','Pwf_{median}','Pwf_{mean}','dPwf','dPwf_{median}','dPwf_{mean}')
grid on
xlabel('t')
ylabel('P')
toc