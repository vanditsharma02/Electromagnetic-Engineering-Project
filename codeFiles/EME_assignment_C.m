%fprintf("2");

A = 1;
B = 1;

NA = 50;
NB = 50;

gaussquad = 3;
H = 1;

eps = 8.85*(10^(-12));
pie = 3.14;
plates = NA*NB;

P = zeros(NB,NA);

fun = @(x,y) 1./(sqrt(x.^2 + y.^2));
q = integral2(fun,-A/(NA*2),A/(NA*2),-B/(NB*2),B/(NB*2));




%disp(q);
%temp = 8*log(1+(2^0.5))*(A/(2*NA));

X = [];
Y = [];

[u,v] = lgwt(gaussquad,-(A/(2*NA)),(A/(2*NA)));
disp(u);
disp(v);

for i = 1:plates
	for j = 1:plates
		
        if j == i
            P(i,j) =(1/((4*pie*eps)))*q;

        else
            temp = 0;
            for k = 1:gaussquad
                for l = 1:gaussquad
                    x1 = [((mod(j-1,NA)+1/2)*(B/NB)+u(k)) ((floor((j-1)/NA)+1/2)*(A/NA)+u(l))];
                    x2 = [((mod(i-1,NA)+1/2)*(B/NB)) ((floor((i-1)/NA)+1/2)*(A/NA))];
                    dist = (((x1(1)-x2(1))^2)+((x1(2)-x2(2))^2))^0.5;
                    temp = temp+ v(k)*v(l)*(1/((4*pie*eps)*(dist)));
                end
            end
            P(i,j) = temp;
        end
    end
    X = [X,(mod(i-1,NA)+1/2)*(B/NB)];
    Y = [Y,(floor((i-1)/NA)+1/2)*(A/NA)];
end

%disp(P);
lin1=[];
lin2=[];
for j = 1:plates
    lin1(j,1)=1/2;
    lin2(j,1)=-1/2;
end
%disp(lin);
P = P./((B/NB)*(A/NA));
Q1 = pinv(P)*lin1;
Q2 = pinv(P)*lin2;

surfQ1=[,];
for i = 1:NA
	for j = 1:NB
        surfQ1(i,j) = Q1(((i-1)*NA)+j);
    end
end       

surfQ2=[,];
for i = 1:NA
	for j = 1:NB
        surfQ2(i,j) = Q2(((i-1)*NA)+j);
    end
end     
        
disp(Q1);
plot3(X,Y,Q1);

figure;
surf(surfQ1);

hold on;
surf(surfQ2);
hold off;
%%%%%%%%%%%%%%%%%%%%%% PART 2-FINDING ELECTRIC FIELD %%%%%%%%%%%%%%%%%%%%%%

[c,d,e] = meshgrid(-1:0.2:2,-1:0.2:2,H/10:H/10:9*(H/10));
u=0;
for i = 1:plates
    u = u +(1/(4*pie*eps)).*(c-X(i)).*Q1(i).*((c-X(i)).^2+(d-Y(i)).^2+(e).^2).^(-1.5);
    u = u +(1/(4*pie*eps)).*(c-X(i)).*Q2(i).*((c-X(i)).^2+(d-Y(i)).^2+(e-H).^2).^(-1.5);
end
v=0;
for i = 1:plates
    v = v +(1/(4*pie*eps)).*(d-Y(i)).*Q1(i).*((c-X(i)).^2+(d-Y(i)).^2+(e).^2).^(-1.5);
    v = v +(1/(4*pie*eps)).*(d-Y(i)).*Q2(i).*((c-X(i)).^2+(d-Y(i)).^2+(e-H).^2).^(-1.5);
end
w=0;
for i = 1:plates
    w = w +(1/(4*pie*eps)).*(e).*Q1(i).*((c-X(i)).^2+(d-Y(i)).^2+(e).^2).^(-1.5);
    w = w +(1/(4*pie*eps)).*(e-H).*Q2(i).*((c-X(i)).^2+(d-Y(i)).^2+(e-H).^2).^(-1.5); 
end


figure;
quiver3(c,d,e,u,v,w,'linewidth',1);


capacitance = 0;
for i = 1:plates
    capacitance = capacitance+Q1(i);
end
capacitance = capacitance;
disp(capacitance);

