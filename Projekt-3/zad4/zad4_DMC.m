clear;
Tp=0.5;

kk=600;
Upp=0;
Ypp=0;

y(1:kk)=Ypp;
u(1:kk)=Upp;
yzad(1:kk)=-0.3;
yzad(150:kk)=10;
yzad(300:kk)=0.5;
yzad(450:kk)=1.5;

Umax = 1;
Umin = -1;
E = 0;

% Ustawienia (wykresy u i y po kolei tak jak podane poniżej parametry)
%N          600,300,300,150,100
%Nu         300,150,50,25,10
%lambda     10,10,5,1,1
D = 600; 
N = 11; 
Nu = 2; 
lambda = 10;
u(20:kk) = Upp+0.01; %%%%%%%%
y(1:kk)=Ypp;
for k=7:kk
    y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2)); 
end
S=(y(21:kk)-Ypp)/(u(end)-u(1));
S(kk-20:kk)=S(end);


M=zeros(N,Nu);
for i=1:N
    for j=1:Nu
        if (i>=j)
            if(i-j+1<=D)
                M(i,j)=S(i-j+1);
            else
                M(i,j)=S(D);
            end
        end
    end
end

Mp = zeros(N, D-1);
for i=1:N
    for j=1:D-1
        if i + j <= D
            Mp(i, j) = S(i + j) - S(j);
        else
            Mp(i, j) = S(D) - S(j);
        end
    end
end

L = lambda * eye(Nu);
K = (M'*M + L)^(-1) * M';
dUp(1:D-1)=0;
dUp=dUp';
u(1:kk) = Upp;
y(1:kk) = Ypp;
yzad(1:kk)=Ypp;
yzad(150:kk)=Ypp+0.01;

for k=7:kk

    y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2));

    for i=1:D-1
        dUp(i)=0;
        if k-i-1>0
            dUp(i)=u(k-i)-u(k-1-i);
        end
    end
     
    Y = y(k) * ones(N,1);
    Yzad = yzad(k) * ones(N, 1);
    dU = K * (Yzad - Y - Mp * dUp);
    
    u(k) = u(k-1) + dU(1);

    u(k) = min(u(k), Umax);
    u(k) = max(u(k), Umin);
    
    E = E + (yzad(k)- y(k))^2;
end

figure;
stairs((1:kk),u);
title("u, E="+ num2str(E)); 
xlabel('k'); 
ylabel('u');
print('pust_dmc_u_5_regulatory_nr_5.png','-dpng','-r400');

figure;
hold on
stairs((1:kk),y);
stairs((1:kk),yzad, ':');
title("yzad, y, E="+ num2str(E)); 
xlabel('k'); 
ylabel('y');
legend('y','yzad');
title("y, E="+ num2str(E))
print('pust_dmc_y_5_regulatory_nr_5.png','-dpng','-r400');