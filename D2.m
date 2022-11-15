%%
%无并行情况的电导计算，全空间k点求和，无能量范围选择（修正后的仅取能量范围在0.05eV内的k点）
%测试1
clear;clc;
tic
knum1=21;
knum2=21;
bnum=2;
kx=pi*linspace(-1,1,knum1);
ky=pi*linspace(-1,1,knum2);
d2k=(kx(2)-kx(1))*(ky(2)-ky(1));
sigx=[0,1;1,0];sigy=[0,-i;i,0];sigz=[1,0;0,-1];
% eigw=zeros(bnum,bnum,knum2,knum1);
% eigvl=zeros(knum2,knum1,bnum);
eigvl=zeros(2,1);
eigw=zeros(2,2);
% Berry=zeros(knum2,knum1);
kt=0.0026*3;%0.0026电子伏特，相当于30K
tau=1;%弛豫时间
%可调参数，用于向ode函数传递参数
ele=0.1;
% B=1;
%电导计算中的两个因子函数，g1,g2
g1=zeros(1,2);
ltime=100;points=10001;
tspan=linspace(0,-ltime,points);%离散化轨道历史时刻
dt=ltime/(points-1);
Bpoints=10;
sigmaB=zeros(2,2,Bpoints);
fermi=2.5;%费米面高度
for bi=1:Bpoints
    B=0.2*bi;
    sigma=zeros(2,2);%电导矩阵
for i=1:knum1
    for j=1:knum2
        k0=[kx(i);ky(j)];%初始时刻的k点
        [t,p] =ode45(@(t,p) korbit2(t,p,ele,B),tspan,k0);%输出轨道历史时刻动量p
        g2=zeros(1,2);
        for l=1:points
        %求解哈密顿量本征方程
            H=(2*cos(p(l,1))+2*cos(p(l,2))-3)*sigz+sin(p(l,1))*sigx-sin(p(l,2))*sigy;
            [w,e]=eig(H);
            eigvl(:)=diag(e);
            eigw(:,:)=w;

            pH=zeros(bnum,bnum,2);   %对哈密顿量的1阶导，后一个指标代表x，y
            ppH=zeros(bnum,bnum,2,2);%对哈密顿量的2阶导，后两个指标代表x，y
            pH(:,:,1)=-2*sin(p(l,1))*sigz+cos(p(l,1))*sigx;
            pH(:,:,2)=-2*sin(p(l,2))*sigz-cos(p(l,2))*sigy;
            pxxH=-2*cos(p(l,1))*sigz-sin(p(l,1))*sigx;
            pyyH=-2*cos(p(l,2))*sigz+sin(p(l,2))*sigy;  
            ppH(:,:,1,1)=pxxH;
            ppH(:,:,2,2)=pyyH;

            ber=-2*imag(eigw(:,2)'*pH(:,:,1)*eigw(:,1)*eigw(:,1)'*pH(:,:,2)*eigw(:,2))./(eigvl(2)-eigvl(1))^2;
            invD=1/(1+ele*B*ber); %压缩系数修正

            eps=[0,1;-1,0];%2维时3阶全反对称张量变为2阶张量
            v=zeros(2,bnum);%2维速度，第一个指标代表x，y

            %以下计算D，V，T
            D=zeros(bnum,bnum,2);     
            Ddag=zeros(bnum,bnum,2);
            V=zeros(bnum,bnum,2);
            Vdag=zeros(bnum,bnum,2);
            T=zeros(bnum,bnum,2,2);%后两个指标代表x,y
            Tdag=zeros(bnum,bnum,2,2);


            for m=1:bnum
                for n=1:bnum
                    if m~=n
                        for index=1:2
                        V(m,n,index)=eigw(:,m)'*pH(:,:,index)*eigw(:,n);
                        D(m,n,index)=V(m,n,index)/(eigvl(n)-eigvl(m));
                        end
                    end
                end
            end

            for index=1:2
                Vdag(:,:,index)=V(:,:,index)';
                Ddag(:,:,index)=D(:,:,index)';
            end

            %计算T
            for a=1:2
                for b=1:2
                    if a==b
                        for m=1:bnum
                           for n=1:bnum
                                if m==n
                                    T(m,m,a,a)=D(m,:,a)*D(:,m,a);
                                else
                                    T(m,n,a,a)=-2/(eigvl(m)-eigvl(n))*(1/2*eigw(:,m)'*ppH(:,:,a,a)*eigw(:,n)+...
                                    V(m,:,a)*D(:,n,a)-V(n,n,a)*D(m,n,a));
                                end                                                                    
                           end
                        end

                    else
                        for m=1:bnum
                            for n=1:bnum
                                if m==n
                                    T(m,m,a,b)=1/2*(D(m,:,a)*D(:,m,b)+D(m,:,b)*D(:,m,a));
                                else
                                    T(m,n,a,b)=-1/(eigvl(m)-eigvl(n))*(eigw(:,m)'*ppH(:,:,a,b)*eigw(:,n)+...
                                    V(m,:,a)*D(:,n,b)+V(m,:,b)*D(:,n,a)-V(n,n,a)*D(m,n,b)-V(n,n,b)*D(m,n,a));
                                end
                            end
                        end        
                    end
                    Tdag(:,:,a,b)=T(:,:,a,b)';
                end
            end

            %计算无磁矩时的速度
            for bn=1:bnum
                for index=1:2
                    v(index,bn)=real(eigw(:,bn)'*pH(:,:,index)*eigw(:,bn));
                end
            end

            %计算磁矩导数，用于修正速度。每次计算时都需要对pm归0化
            pm=zeros(2,bnum);%z方向磁矩偏导计算
            for bn=1:bnum %计算bn带的磁矩
                for mu=1:2 %磁矩导数的两个分量指标
                    for a=1:2
                        for b=1:2
                            pm(mu,bn)=pm(mu,bn)-1i*ele/2*eps(a,b)*(-Tdag(bn,:,mu,a)*V(:,bn,b)+Ddag(bn,:,a)*V(:,:,mu)*D(:,bn,b)-...
                            V(bn,bn,mu)*Ddag(bn,:,a)*D(:,bn,b)-Vdag(bn,:,a)*T(:,bn,mu,b));
                        end
                    end
                    v(mu,bn)=v(mu,bn)-B*pm(mu,bn);%至此计算出修正速度
                end
            end

            if l==1
              %计算初始点（轨道终点）的磁矩，用以计算能量。 
              M=-ele*imag(eigw(:,2)'*pH(:,:,1)*eigw(:,1)*eigw(:,1)'*pH(:,:,2)*eigw(:,2))./(eigvl(2)-eigvl(1));
              energy=eigvl(2)-B*M-fermi;
              g1(:)=v(:,2); %由于磁场垂直于二维面，这里的第一个因子函数等于\tilde{v} 
              s=1-1/(1+exp(energy/kt));
              pef=-1/kt*s/(1+exp(energy/kt));%费米分布函数求导，energy
            end
            %接下来计算g2，需要乘上衰减因子后对时间积分（离散求和），只计算上能带，即能带2
            g2(:)=g2(:)+dt/tau*exp(t(l)/tau)*v(:,2)*invD;
        end
        
        for a=1:2
           for b=1:2
             sigma(a,b)=sigma(a,b)+ele^2*tau/(2*pi)^2*d2k*pef*g1(a)*g2(b); 
           end
        end
        sigma
        
    end
end
sigmaB(:,:,bi)=sigma;
end
toc
%%
%拓扑非平凡，有反带，修正前后的能带计算以及贝里曲率和磁矩绘制
clear;clc;
bnum=2;%2能级系统
mu=2;%费米能
tau=1; %弛豫时间
B=20;%磁场强度
ele=0.1;%电子电量
%在2维第一布里渊区内均匀取点
knum1=101;
knum2=101;
kx=pi*linspace(-1,1,knum1);
ky=pi*linspace(-1,1,knum2);

sigx=[0,1;1,0];sigy=[0,-1i;1i,0];sigz=[1,0;0,-1];
eigw=zeros(bnum,bnum,knum2,knum1);
eigvl=zeros(knum2,knum1,bnum);
Berry=zeros(knum2,knum1); %暂时只考虑上能带贝里曲率
pH=zeros(bnum,bnum,2);
etidle=zeros(knum2,knum1);
M=zeros(knum2,knum1);
for i=1:knum1  
    for j=1:knum2
        H=(2*cos(kx(i))+2*cos(ky(j))-3)*sigz+sin(kx(i))*sigx-sin(ky(j))*sigy;
        [w,e]=eig(H);
        eigvl(j,i,:)=diag(e);
        eigw(:,:,j,i)=w; 
        
        pH(:,:,1)=-2*sin(kx(i))*sigz+cos(kx(i))*sigx;
        pH(:,:,2)=-2*sin(ky(j))*sigz-cos(ky(j))*sigy;
        Berry(j,i)=-2*imag(eigw(:,2,j,i)'*pH(:,:,1)*eigw(:,1,j,i).*eigw(:,1,j,i)'*pH(:,:,2)*eigw(:,2,j,i))./(eigvl(j,i,1)-eigvl(j,i,2))^2;      
        M(j,i)=-ele*imag(eigw(:,2,j,i)'*pH(:,:,1)*eigw(:,1,j,i)*eigw(:,1,j,i)'*pH(:,:,2)*eigw(:,2,j,i))/(eigvl(j,i,2)-eigvl(j,i,1));
        etidle(j,i)=eigvl(j,i,2)-B*M(j,i);
    end 
end
figure(1)
surf(kx,ky,eigvl(:,:,2));
hold on;
surf(kx,ky,1.38*ones(knum2,knum1),'FaceAlpha',0.9);
figure(2)
surf(kx,ky,etidle(:,:));
hold on;
surf(kx,ky,1.38*ones(knum2,knum1),'FaceAlpha',0.9);
% surf(kx,ky,0.496*ones(knum2,knum1),'FaceAlpha',0.9);
figure(3)
surf(kx,ky,Berry);
figure(4)
surf(kx,ky,M)
%%
%数值解解出的磁矩偏导数
clear;clc;
tic
knum1=101;
knum2=101;
ele=0.1;
bnum=2;
kx=pi*linspace(-1,1,knum1);
ky=pi*linspace(-1,1,knum2);
sigx=[0,1;1,0];sigy=[0,-i;i,0];sigz=[1,0;0,-1];
m=zeros(knum2,knum1);
pm=zeros(knum2,knum1,2,bnum);
eigvl=zeros(2,1);
eigw=zeros(2,2);
eps=[0,1;-1,0];%2维时3阶全反对称张量变为2阶张量
for i=1:knum1
    for j=1:knum2
            H=(2*cos(kx(i))+2*cos(ky(j))-3)*sigz+sin(kx(i))*sigx-sin(ky(j))*sigy;
            [w,e]=eig(H);
            eigvl(:)=diag(e);
            eigw(:,:)=w*exp(4i*(kx(i)+ky(j)))*sigz;

            pH=zeros(bnum,bnum,2);   %对哈密顿量的1阶导，后一个指标代表x，y
            ppH=zeros(bnum,bnum,2,2);%对哈密顿量的2阶导，后两个指标代表x，y
            pH(:,:,1)=-2*sin(kx(i))*sigz+cos(kx(i))*sigx;
            pH(:,:,2)=-2*sin(ky(j))*sigz-cos(ky(j))*sigy;
            pxxH=-2*cos(kx(i))*sigz-sin(kx(i))*sigx;
            pyyH=-2*cos(ky(j))*sigz+sin(ky(j))*sigy;  
            ppH(:,:,1,1)=pxxH;
            ppH(:,:,2,2)=pyyH;
            
            v=zeros(2,bnum);%2维速度，第一个指标代表x，y

            %以下计算D，V，T
            D=zeros(bnum,bnum,2);     
            Ddag=zeros(bnum,bnum,2);
            V=zeros(bnum,bnum,2);
            Vdag=zeros(bnum,bnum,2);
            T=zeros(bnum,bnum,2,2);%后两个指标代表x,y
            Tdag=zeros(bnum,bnum,2,2);
            VD=zeros(bnum,bnum,2);
            VDdag=zeros(bnum,bnum,2);


            for m=1:bnum
                for n=1:bnum
                    if m~=n
                        for index=1:2
                        V(m,n,index)=eigw(:,m)'*pH(:,:,index)*eigw(:,n);
                        D(m,n,index)=V(m,n,index)/(eigvl(n)-eigvl(m));
                        end
                    else 
                        for index=1:2
                          V(m,m,index)=eigw(:,m)'*pH(:,:,index)*eigw(:,m);
                        end  
                    end
                end
            end

            for index=1:2
                Vdag(:,:,index)=V(:,:,index)';
                Ddag(:,:,index)=D(:,:,index)';
            end

            %计算T
            for a=1:2
                for b=1:2
                    if a==b
                        for m=1:bnum
                           for n=1:bnum
                                if m==n
%                                     T(m,m,a,b)=1/2*(D(m,:,a)*D(:,m,b)+D(m,:,b)*D(:,m,a));
                                    T(m,m,a,a)=D(m,:,a)*D(:,m,a);
                                else
                                    T(m,n,a,a)=-2/(eigvl(m)-eigvl(n))*(1/2*eigw(:,m)'*ppH(:,:,a,a)*eigw(:,n)+...
                                    V(m,:,a)*D(:,n,a)-V(n,n,a)*D(m,n,a));
                                end                                                                    
                           end
                        end

                    else
                        for m=1:bnum
                            for n=1:bnum
                                if m==n
                                    T(m,m,a,b)=1/2*(D(m,:,a)*D(:,m,b)+D(m,:,b)*D(:,m,a));
                                else
                                    T(m,n,a,b)=-1/(eigvl(m)-eigvl(n))*(eigw(:,m)'*ppH(:,:,a,b)*eigw(:,n)+...
                                    V(m,:,a)*D(:,n,b)+V(m,:,b)*D(:,n,a)-V(n,n,a)*D(m,n,b)-V(n,n,b)*D(m,n,a));
                                end
                            end
                        end        
                    end
                    Tdag(:,:,a,b)=T(:,:,a,b)';
                end
            end

            %计算无磁矩时的速度
            for bn=1:bnum
                for index=1:2
                    v(index,bn)=real(eigw(:,bn)'*pH(:,:,index)*eigw(:,bn));
                end
            end

            %计算磁矩导数，用于修正速度。每次计算时都需要对pm归0化
            
            for index=1:2
                for m=1:2
                    for n=1:2
                     VD(m,n,index)=(eigvl(m)-eigvl(n))*D(m,n,index);  
                    end
                end
            end
        
        for index=1:2
            VDdag(:,:,index)=VD(:,:,index)';
        end
            
               bn=2;
                for mu=1:2 %磁矩导数的两个分量指标
                    for a=1:2
                        for b=1:2
                            pm(j,i,mu,bn)=pm(j,i,mu,bn)-1i*ele/2*eps(a,b)*(-Tdag(bn,:,mu,a)*V(:,bn,b)+Ddag(bn,:,a)*V(:,:,mu)*D(:,bn,b)-...
                            V(bn,bn,mu)*Ddag(bn,:,a)*D(:,bn,b)-Vdag(bn,:,a)*T(:,bn,mu,b));
%                             pm(j,i,mu,bn)=pm(j,i,mu,bn)-1i*ele/2*eps(a,b)*(Tdag(bn,:,mu,a)*VD(:,bn,b)+Ddag(bn,:,a)*V(:,:,mu)*D(:,bn,b)-...
%                             V(bn,bn,mu)*Ddag(bn,:,a)*D(:,bn,b)+VDdag(bn,:,a)*T(:,bn,mu,b));

                        end
                    end
                end
    end
     
end
figure(1)
surf(kx,ky,real(pm(:,:,1,2)));
figure(2)
surf(kx,ky,real(pm(:,:,2,2)));
toc
%%
%解析解计算的磁矩以及磁矩偏导数
clear;clc;
B=10;%磁场强度
ele=0.1;%电子电量
%在2维第一布里渊区内均匀取点
knum1=101;
knum2=101;
m=zeros(knum2,knum1);
pxm=zeros(knum2,knum1);
pym=zeros(knum2,knum1);
kx=pi*linspace(-1,1,knum1);
ky=pi*linspace(-1,1,knum2);
for i=1:knum1
   for j=1:knum2
       m(j,i)=-ele*imag(Moment(kx(i),ky(j)));
       pxm(j,i)=-ele*imag(pxM(kx(i),ky(j)));
       pym(j,i)=-ele*imag(pyM(kx(i),ky(j)));
       
   end
end
figure(6)
surf(kx,ky,pxm);
figure(7)
surf(kx,ky,pym);
figure(8)
surf(kx,ky,m)

%%
%拓扑平凡，没有反带
clear;clc;
bnum=2;%2能级系统
mu=2;%费米能
tau=1; %弛豫时间
B=10;%磁场强度
ele=0.1;%电子电量
%在2维第一布里渊区内均匀取点
knum1=101;
knum2=101;
kx=pi*linspace(-1,1,knum1);
ky=pi*linspace(-1,1,knum2);

sigx=[0,1;1,0];sigy=[0,-1i;1i,0];sigz=[1,0;0,-1];
eigw=zeros(bnum,bnum,knum2,knum1);
eigvl=zeros(knum2,knum1,bnum);
Berry=zeros(knum2,knum1); %暂时只考虑上能带贝里曲率
% D=zeros(bnum,bnum,2);     %压缩系数
% Ddag=zeros(bnum,bnum,2);
% V=zeros(bnum,bnum,2);
% Vdag=zeros(bnum,bnum,2);
% T=zeros(bnum,bnum,2,2);%后两个指标代表x,y
% Tdag=zeros(bnum,bnum,2,2);
% ppH=zeros(bnum,bnum,2,2);%对哈密顿量的2阶导，后两个指标代表x，y
% v=zeros(2,bnum);%2维速度，第一个指标代表x，y
pH=zeros(bnum,bnum,2);
eps=[0,1;-1,0];%2维时3阶全反对称张量变为2阶张量
etidle=zeros(knum2,knum1);
M=zeros(knum2,knum1);
for i=1:knum1  
    for j=1:knum2
        H=(2*cos(kx(i))+2*cos(ky(j))-3)*sigz+sin(kx(i))*sigx-sin(ky(j))*sigy;
        [w,e]=eig(H);
        eigvl(j,i,:)=diag(e);
        eigw(:,:,j,i)=w; 
        
        pH(:,:,1)=-2*sin(kx(i))*sigz+cos(kx(i))*sigx;
        pH(:,:,2)=-2*sin(ky(j))*sigz-cos(ky(j))*sigy;
        Berry(j,i)=-2*imag(eigw(:,2,j,i)'*pH(:,:,1)*eigw(:,1,j,i).*eigw(:,1,j,i)'*pH(:,:,2)*eigw(:,2,j,i))./(eigvl(j,i,1)-eigvl(j,i,2))^2;      
        M(j,i)=-ele*imag(eigw(:,2,j,i)'*pH(:,:,1)*eigw(:,1,j,i)*eigw(:,1,j,i)'*pH(:,:,2)*eigw(:,2,j,i))/(eigvl(j,i,2)-eigvl(j,i,1));
        etidle(j,i)=eigvl(j,i,2)-B*M(j,i);
    end 
end
figure(1)
surf(kx,ky,eigvl(:,:,2));
figure(2)
surf(kx,ky,etidle(:,:));
hold on;
surf(kx,ky,0.496*ones(knum2,knum1),'FaceAlpha',0.9);
figure(3)
surf(kx,ky,Berry);
figure(4)
surf(kx,ky,M*B)

%%
% hold on;
% surf(kx,ky,eigvl(:,:,2));
% surf(kx,ky,0.818*ones(knum2,knum1),'FaceAlpha',0.8);
% surf(kx,ky,1.38*ones(knum2,knum1),'FaceAlpha',1);
%%
k=find(eigvl(:,:,2)<0.818)
mat=eigvl(:,:,2);
mat(k)
%%
%寻找落在基准能量上的k点散点，观察轨道形状
error=0.01;
[x,y]=find(abs(etidle(:,:)-0.496)<error);
l=length(x);
% surf(kx(x),ky(y),etidle(x,y));
ener=zeros(1,l);
for i=1:l
    ener(i)=etidle(x(i),y(i));
end
figure(10)
scatter3(kx(x),ky(y),ener);
%%
%拓扑非平凡，有反带，在导带底部2次展开，进行能带计算以及贝里曲率和磁矩绘制
clear;clc;
bnum=2;%2能级系统
mu=2;%费米能
tau=1; %弛豫时间
B=1;%磁场强度
ele=0.1;%电子电量
%在2维第一布里渊区内均匀取点
knum1=101;
knum2=101;
kx=pi*linspace(-1,1,knum1);
ky=pi*linspace(-1,1,knum2);

sigx=[0,1;1,0];sigy=[0,-1i;1i,0];sigz=[1,0;0,-1];
eigw=zeros(bnum,bnum,knum2,knum1);
eigvl=zeros(knum2,knum1,bnum);
Berry=zeros(knum2,knum1); %暂时只考虑上能带贝里曲率
pH=zeros(bnum,bnum,2);
etidle=zeros(knum2,knum1);
M=zeros(knum2,knum1);
for i=1:knum1  
    for j=1:knum2
%         H=(2*cos(kx(i))+2*cos(ky(j))-3)*sigz+sin(kx(i))*sigx-sin(ky(j))*sigy;
        H=(1-kx(i)^2-ky(j)^2)*sigz+kx(i)*sigx-ky(j)*sigy;
        [w,e]=eig(H);
        eigvl(j,i,:)=diag(e);
        eigw(:,:,j,i)=w; 
        
%         pH(:,:,1)=-2*sin(kx(i))*sigz+cos(kx(i))*sigx;
%         pH(:,:,2)=-2*sin(ky(j))*sigz-cos(ky(j))*sigy;
        pH(:,:,1)=-2*kx(i)*sigz+1*sigx;
        pH(:,:,2)=-2*ky(j)*sigz-1*sigy;
        Berry(j,i)=-2*imag(eigw(:,2,j,i)'*pH(:,:,1)*eigw(:,1,j,i).*eigw(:,1,j,i)'*pH(:,:,2)*eigw(:,2,j,i))./(eigvl(j,i,1)-eigvl(j,i,2))^2;      
        M(j,i)=-ele*imag(eigw(:,2,j,i)'*pH(:,:,1)*eigw(:,1,j,i)*eigw(:,1,j,i)'*pH(:,:,2)*eigw(:,2,j,i))/(eigvl(j,i,2)-eigvl(j,i,1));
        etidle(j,i)=eigvl(j,i,2)-B*M(j,i);
    end 
end
figure(1)
surf(kx,ky,eigvl(:,:,2));
hold on;
surf(kx,ky,1.38*ones(knum2,knum1),'FaceAlpha',0.9);
figure(2)
surf(kx,ky,etidle(:,:));
hold on;
% surf(kx,ky,0.496*ones(knum2,knum1),'FaceAlpha',0.9);
figure(3)
surf(kx,ky,Berry);
figure(4)
surf(kx,ky,M)
%%
%电子空穴补偿模型
clear;clc;
bnum=2;%2能级系统
mu=2;%费米能
tau=1; %弛豫时间
B=10;%磁场强度
ele=0.1;%电子电量
%在2维第一布里渊区内均匀取点
knum1=101;
knum2=101;
kx=pi*linspace(-0.5,1.5,knum1);
ky=pi*linspace(-1,1,knum2);

sigx=[0,1;1,0];sigy=[0,-1i;1i,0];sigz=[1,0;0,-1];
eigw=zeros(bnum,bnum,knum2,knum1);
eigvl=zeros(knum2,knum1,bnum);
Berry=zeros(knum2,knum1); %暂时只考虑上能带贝里曲率
pH=zeros(bnum,bnum,2);
etidle=zeros(knum2,knum1);
M=zeros(knum2,knum1);
for i=1:knum1  
    for j=1:knum2
        H=-2*cos(kx(i))*eye(2)-(2*cos(ky(j))-3)*sigz;
        [w,e]=eig(H);
        eigvl(j,i,:)=diag(e);
        eigw(:,:,j,i)=w; 
        
%         pH(:,:,1)=-2*sin(kx(i))*sigz+cos(kx(i))*sigx;
%         pH(:,:,2)=-2*sin(ky(j))*sigz-cos(ky(j))*sigy;
%         Berry(j,i)=-2*imag(eigw(:,2,j,i)'*pH(:,:,1)*eigw(:,1,j,i).*eigw(:,1,j,i)'*pH(:,:,2)*eigw(:,2,j,i))./(eigvl(j,i,1)-eigvl(j,i,2))^2;      
%         M(j,i)=-ele*imag(eigw(:,2,j,i)'*pH(:,:,1)*eigw(:,1,j,i)*eigw(:,1,j,i)'*pH(:,:,2)*eigw(:,2,j,i))/(eigvl(j,i,2)-eigvl(j,i,1));
%         etidle(j,i)=eigvl(j,i,2)-B*M(j,i);
    end 
end
figure(1)
surf(kx,ky,eigvl(:,:,2));
hold on;
surf(kx,ky,eigvl(:,:,1));
surf(kx,ky,0*ones(knum2,knum1),'FaceAlpha',0.9);
xlabel('kx');
ylabel('ky');
% figure(2)
% surf(kx,ky,etidle(:,:));
% hold on;
% surf(kx,ky,1.38*ones(knum2,knum1),'FaceAlpha',0.9);
% % surf(kx,ky,0.496*ones(knum2,knum1),'FaceAlpha',0.9);
% figure(3)
% surf(kx,ky,Berry);
% figure(4)
% surf(kx,ky,M)