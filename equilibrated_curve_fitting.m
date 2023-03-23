function[fx_p,fx]= equilibrated_curve_fitting()
% Curve fitting through equilibrated error in each element
%********* Input **********************
% f : f(x) function to be fitted 
% p : Order of approximation for curve fitting
% n : Number of segments required to fit the curve, which will change adaptively so that error comes under the desired tolerance 'tol'
% x1st : first point (or the left boundarary) of the domain
% xlast : Last point (or the right boundarary) of the domain
% domain :[x1st,xlast];
x=sym('x');
f=x^5-x^4+x^3-x^2+x-1; %sin(7*x);
x1st=0;
xlast=2*pi;
p=2;%Order of approximation
n=2;% no of segments(Elements)
nitr=8;
tol=10^-5;% tolerence while equilibrating error  tol=0.001%=10^-5
Jzero=10^-5;
lw=1.5;%line width in plot
%*********Initial Calculations ***********
df=diff(f,x);
L=xlast-x1st; % length in x-dirn
J_f=int(((f)^2+(df)^2),x,x1st,xlast);
Jt_d=double(J_f*tol);
errtol=Jt_d*10^-2;% diff in previous Jtotal and Jtotal 1% of desired Jtotal
%********Initialization *****************
x_p=sym(zeros(p+1,1));
hd=(zeros(n,1));
K=zeros(p+1,p+1);R=zeros(p+1,1);
Y_i=zeros(p-1,1,n);
K_io=zeros(p-1,2,n);
fx_app=zeros(n,p);%approximated functional values in each segment
fx_points=zeros(n,p+1);
signal=0;cal=0;Itr=1;le=2;
Jtotal=Jt_d+1;
%********** Ploting function ************
hold on 
fp=ezplot(f,[x1st,xlast]);
set(fp,'linewidth',lw);
set(fp,'color',1/255*[0,150,100]);
legendInfo{1}='fn';
%% $$$$$ Main program $$$$$$$$$$$$$$
% plotStyle = {'-.m+',':bx','--r^','-.co','rd','mp','bh','k<'};
plotStyle = {'m+','bx','r^','co','rd','mp','bh','k<'};
e1=sym('x1');e2=sym('x2'); 
[Ns,dNs]=localBasisfn();
%=======================================Main programe =======================================
while (Jtotal>Jt_d)
    xp=linspace(x1st,xlast,n+1);% nodes in x-dirn
    yp=subs(f,xp);xp1=xp;yp1=yp;
    h0=ones(n,1)*(xlast-x1st)/n;
    elsize(Itr,1)=log(h0(1));
    fprintf('n = %d',n); fprintf('\n');
    for itr=1:nitr    
        Kg=zeros(n+1,n+1);
        Rg=zeros(n+1,1);
        for e=1:n % loop over no of elements/segments
            xp_seg=[xp(e),xp(e+1)];
            phi=subs(Ns,[e1,e2],xp_seg);% basis fun for segment with end points xp_seg
            dphi=subs(dNs,[e1,e2],xp_seg);
            [Klocal,Rlocal,Yi_b,Kio_b]=segmentmatrices(phi,dphi,xp_seg);% using static condensation ***********
            % ******** Assembly  **************
            Kg(e,e)=Kg(e,e)+Klocal(1,1);           %        Connctivity like 1st order  
            Kg(e,e+1)=Kg(e,e+1)+Klocal (1,2); %         L(element,localnode)=(element-1)*p+(localnode) 
            Kg(e+1,e)=Kg(e+1,e)+Klocal (2,1); % -->> L(e,1)=e; L(e,2)=e+1;
            Kg(e+1,e+1)=Kg(e+1,e+1)+Klocal (2,2);
            Rg(e)=Rg(e)+Rlocal(1);
            Rg(e+1)=Rg(e+1)+Rlocal(2);
            Y_i(:,:,e)=Yi_b;
            K_io(:,:,e)=Kio_b;
        end
        yap= approx_solution(Rg,Kg,yp); % Value of function at nodal locations
        %********************************
        Jtotal=0;
        for el=1:n % loop over no of elements/segments
            xp_seg=[xp(el),xp(el+1)];
            Yo=[yap(el),yap(el+1)]';
            if p>1
                Yi=Y_i(:,:,el)-K_io(:,:,el)*Yo; %back substitution  
            else
                Yi=[];
            end
            Y=[Yo(1);Yi;Yo(2)];
            phi=subs(Ns,[e1,e2],xp_seg);% basis fun for segment with end points xp_seg; 
            fbar=Y'*phi;
            dfbar=diff(fbar,x);
            J=double(int(((f-fbar)^2+(df-dfbar)^2),x,xp(el),xp(el+1)));
            Jtotal=Jtotal+J;
            if J<Jzero
                skl=1.0;
            else
                skl=((n*J/tol)^(1/(2*p+1)));%scale
                %skl=(J/Jt_d)^(1/(2*p+1));%scale
            end
            hd(el)=h0(el)/skl;
            fx_app(el,:)=[Yo(1);Yi];
            fx_points(el,:)=(linspace(xp(el),xp(el+1),p+1));
        end       
        Lbar=sum(hd);
        hd_s=(hd*L/Lbar);%scaled h desired
        err=Jtotal;
        %********************************
        if err<=Jt_d
            cal=1;% mark the itration when Jtotal statrs <=Jtotal desired
            le=le+1;
            if le==3%Plot interpolation fit when Jtotl <= Jtotal desired for the 1st time in itr
                plot(xp1,yp1,'-.gv','LineWidth',lw,'MarkerSize',10)
                legendInfo{2}='itr1 or lnitrfit ';
            end
        end
        if itr==1&&cal==1
            signal=1;
        elseif itr >=2
            differr=errold-err;
            if abs(differr)<errtol||(itr>2&&abs(differr)~=differr)% error < tolrence OR diff in error changes sign
                signal=2;
            end
        end
        if (cal==1&&signal~=2)||(le==3)%first appearance of Jtotal<=Jt_d
            fx=double([reshape(fx_app',[n*p,1]);yp(end)]);
            fx_p=unique(reshape(fx_points',[n*(p+1),1]));
            if itr~=1% already ploted for itr=1
                plot(xp,yap,plotStyle{itr},'LineWidth',lw,'MarkerSize',10);
                legendInfo{le}=strcat('itr',num2str(itr));
            end
        end
        if signal==2
            break
        end
        fprintf('Itr = %d',Itr);fprintf('; itr = %d',itr);
        fprintf('; Jtotal = %d',err);fprintf('\n');
        if signal==1||(itr>=2&&cal==1)
            break
        end
        %********************************
        for i=1:n
            xp(i+1)=xp(i)+hd_s(i);
        end
        h0=hd_s;
        errold=err;
    end
    %******* required no of segments *******
    Err(Itr,1)=log(err);
    if Jtotal>Jt_d 
        nr=ceil(n*(Jtotal/Jt_d)^(1/(2*p+1)));%!!! depends on previous itrations
        n=nr;
    end
    signal=0;
    Itr=Itr+1;
    fprintf('\n');
end
%================================ Curve and error plot ========================================
le=length(legendInfo)+1;
plot(fx_p,fx,':k*','LineWidth',1.45,'MarkerSize',9.5);% plotting fitted curve
legendInfo{le}='approx curve ';
legend(legendInfo);legend('location','northeast');legend('boxoff');
hold off
figure
plot(-elsize,Err,'-b*')
title('Improvement of error norm with increasing element number')
xlabel('log of inverse of element size') % x-axis label
ylabel('log of total error norm') % y-axis label
fprintf('Last Jtotal = %d',err);fprintf('; Desired Jtotal = %d',Jt_d);
%% $$$$$$$$$ Lobal basis fn(N) $$$$$$$$  
    function [N,dN]=localBasisfn()
        N=sym(zeros(p+1,1));dN=sym(zeros(p+1,1));
        Lel=e2-e1;% Length of physical element
        x_p(1)= e1; % the first physical node 
        for ip=1:p
            x_p(ip+1)=x_p(1)+ip*Lel/p; % Ponits corresponding to internal nodes in the physical element based on 'p'
        end
        for r=1:p+1 
            N(r)=1;
            for s=1:p+1        
                if(s~=r)  % This will be executed only for j not equal to i (as per lagrangian polynominals definition)
                    N(r)=N(r)*(x-x_p(s))/(x_p(r)-x_p(s));            
                end        
            end
            dN(r)=diff(N(r),x);  
        end
    end
%% $$$$$$$$ Local Matrices $$$$$$$$$$$
    function [Ksb,Rsb,Yi_b,Kio_b]=segmentmatrices(N,dN,xne)
        for  kl=1:p+1
            for jl=1:p+1
                K(kl,jl)=int((dN(kl)*dN(jl)+N(kl)*N(jl)),x,xne(1),xne(2));
            end
            R(kl)=int((df*dN(kl)+f*N(kl)),x,xne(1),xne(2));
        end
        [Ksb,Rsb,Yi_b,Kio_b]=localmatrices(K,R);
        % Matrices for Static Condensation in the segment
        function [Koo_bar,Ro_bar,Yi_bar,Kio_bar]=localmatrices(K,R)
            Koo=[K(1,1),K(1,p+1);K(p+1,1),K(p+1,p+1)];
            Ro=[R(1);R(p+1)];
            Koi=K([1,p+1],2:p);Kio=K(2:p,[1,p+1]);
            Kii=K(2:p,2:p);
            Ri=R(2:p);
            if isempty(Kii)
                Kii=1;
            end
            Yi_bar=(Kii^-1)*Ri;
            Kio_bar=(Kii^-1)*Kio;
            Koo_bar=Koo-Koi*Kio_bar;
            Ro_bar=Ro-Koi*Yi_bar;
        end
    end
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    function ya= approx_solution(F,KofB,yp)
        %*****applying end condition*********
        for ii=1:n+1
            F(ii)=F(ii)-KofB(ii,1)*yp(1)-KofB(ii,n+1)*yp(n+1);
        end
        F(1)=yp(1);
        F(n+1)=yp(n+1);
        KofB(1,:)=0;KofB(n+1,:)=0;
        KofB(:,1)=0;KofB(:,n+1)=0;
        KofB(1,1)=1;KofB(n+1,n+1)=1;
        %******solution********************
        ya=linsolve(KofB,F);% y approximate
    end
end