function y=HermiteSplineInterpolation(x,tx,ty)

%coefficineti di Hermite
h00=@(t) (1+2*t)*(1-t)^2;
h10=@(t) t*(1-t)^2;
h01=@(t) t^2*(3-2*t);
h11=@(t) t^2*(t-1);

%tempo mappato per intervalli qualsiasi
tm=@(x,x0,x1) (x-x0)/(x1-x0);

%polinomio interpolante
p=@(t,x0,x1,y0,y1,m0,m1) h00(tm(t,x0,x1))*y0+h10(tm(t,x0,x1))*(x1-x0)*m0+h01(tm(t,x0,x1))*y1+h11(tm(t,x0,x1))*(x1-x0)*m1;

%derivata sui campioni, a differenze finite
m=@(y_1,y0,y1,t_1,t0,t1) 0.5*((y1-y0)/(t1-t0)+(y0-y_1)/(t0-t_1));

y=zeros(1,length(ty));

j=1;
mi0=m(0,0,x(1),tx(1)-2*tx(2),tx(1)-tx(2),tx(1)); %approssimazione
mi1=m(0,x(1),x(2),tx(1)-tx(2),tx(1),tx(2)); %approssimazione: t_-1-t_0=t_0-t_1 e x_-1=0
while(ty(j)<tx(1))
    y(j)=p(ty(j),tx(1)-tx(2),tx(1),0,x(1),mi0,mi1);
    j=j+1;
end

mi0=mi1;
mi1=m(x(1),x(2),x(3),tx(1),tx(2),tx(3));
while(ty(j)<tx(2))
    y(j)=p(ty(j),tx(1),tx(2),x(1),x(2),mi0,mi1);
    j=j+1;
end

mi0=m(x(1),x(2),x(3),tx(1),tx(2),tx(3));
for i=2:length(tx)-2
    mi1=m(x(i),x(i+1),x(i+2),tx(i),tx(i+1),tx(i+2));
    while(ty(j)<tx(i+1))
        y(j)=p(ty(j),tx(i),tx(i+1),x(i),x(i+1),mi0,mi1);
        j=j+1;
    end
    mi0=mi1;    
end

i=length(tx)-1;
dt=(tx(i+1)-tx(i)); %stima dt da ultimi campioni
mi1=m(x(i),x(i+1),0,tx(i),tx(i+1),tx(i+1)+dt); %approssimazione
while(ty(j)<tx(i+1))
    y(j)=p(ty(j),tx(i),tx(i+1),x(i),x(i+1),mi0,mi1);
    j=j+1;
end

i=i+1;
mi0=mi1;
mi1=m(x(i),0,0,tx(i),tx(i)+dt,tx(i)+2*dt); %approssimazione
while(j<=length(ty))
    y(j)=p(ty(j),tx(i),tx(i)+dt,x(i),0,mi0,mi1);
    j=j+1;
end

end