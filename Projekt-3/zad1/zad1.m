kk=200;

Upp=0;
Ypp=0;

u(1:kk)=Upp;
y(1:kk)=Ypp;

Umin=-1;
Umax=1;

for k=7:kk
    y(k)=symulacja_obiektu15y_p3(u(k-5),u(k-6),y(k-1),y(k-2));
end

figure
stairs(1:kk,u)
title('Sterowanie w punkcie pracy:')
xlabel('k')
ylabel('u')
print('zad1_u.png','-dpng','-r400')

figure
stairs(1:kk,y)
title('Wyjście w punkcie pracy:')
xlabel('k')
ylabel('y')
print('zad1_y.png','-dpng','-r400')