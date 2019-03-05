x = [];
y = [];

for i=1:2:41
     y = [y,cap(i,1,1,3)];
     x= [x,i];
end    
plot(x,y);
    

x = [];
y = [];

for i=1:2:11
     y = [y,cap(30,1,1,i)];
     x= [x,i];
end
plot(x,y);
    