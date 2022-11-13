function dkdt = korbit2(t,k,ele,B)
        bnum=2;
        sigx=[0,1;1,0];sigy=[0,-1i;1i,0];sigz=[1,0;0,-1];   
        eigvl=zeros(1,2);
        eigw=zeros(2,2);
        H=(2*cos(k(1))+2*cos(k(2))-3)*sigz+sin(k(1))*sigx-sin(k(2))*sigy;
        [w,e]=eig(H);
        eigvl(:)=diag(e);
        eigw(:,:)=w; 
%         eigw(:,:)=w*exp(1i*(k(1)+k(2)))*sigz; %这里对上下能带的波函数做规范变换，k轨道结果不变，说明不依赖于规范
        
        pH=zeros(bnum,bnum,2);   %对哈密顿量的1阶导，后一个指标代表x，y
        ppH=zeros(bnum,bnum,2,2);%对哈密顿量的2阶导，后两个指标代表x，y
        pH(:,:,1)=-2*sin(k(1))*sigz+cos(k(1))*sigx;
        pH(:,:,2)=-2*sin(k(2))*sigz-cos(k(2))*sigy;
        pxxH=-2*cos(k(1))*sigz-sin(k(1))*sigx;
        pyyH=-2*cos(k(2))*sigz+sin(k(2))*sigy;  
        ppH(:,:,1,1)=pxxH;
        ppH(:,:,2,2)=pyyH;

        ber=-2*imag(eigw(:,2)'*pH(:,:,1)*eigw(:,1)*eigw(:,1)'*pH(:,:,2)*eigw(:,2))./(eigvl(1)-eigvl(2))^2;       
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
        VD=zeros(bnum,bnum,2);
        VDdag=zeros(bnum,bnum,2);

        
        for m=1:bnum
            for n=1:bnum
                if n~=m
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
                                T(m,m,a,a)=D(m,:,a)*D(:,m,a);
%                                 T(m,m,a,a)=0;
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
%                                 T(m,m,a,b)=D(m,:,b)*D(:,m,a); 
%                                 T(m,m,a,b)=0;
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
        
        for index=1:2
           for m=1:2
               for n=1:2
                VD(m,n,index)=(eigvl(m)-eigvl(n))*D(m,n,index);  
               end
           end
           VDdag(:,:,index)=VD(:,:,index)';
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
%                         pm(mu,bn)=pm(mu,bn)-1i*ele/2*eps(a,b)*(-Tdag(bn,:,mu,a)*V(:,bn,b)+Ddag(bn,:,a)*V(:,:,mu)*D(:,bn,b)-...
%                         V(bn,bn,mu)*Ddag(bn,:,a)*D(:,bn,b)-Vdag(bn,:,a)*T(:,bn,mu,b));
                        pm(mu,bn)=pm(mu,bn)-1i*ele/2*eps(a,b)*(Tdag(bn,:,mu,a)*VD(:,bn,b)+Ddag(bn,:,a)*V(:,:,mu)*D(:,bn,b)-...
                        V(bn,bn,mu)*Ddag(bn,:,a)*D(:,bn,b)+VDdag(bn,:,a)*T(:,bn,mu,b));
                    end
                end
                v(mu,bn)=v(mu,bn)-B*real(pm(mu,bn));%这里需要取实部，虚部可能来自于计算误差
%                 v(mu,bn)=v(mu,bn)-B*pm(mu,bn);
            end
        end
        %只计算上能带
        f=-invD*ele*v(2,2)*B;
        g= invD*ele*v(1,2)*B;
dkdt = [f;g];