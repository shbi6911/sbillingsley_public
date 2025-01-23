metal = "Alum25";
y=zeros(1,8);
tss=zeros(1,8);

for ii=1:8
    [tss(ii),y(ii)]=steadystate(const.data.(metal).time,const.data.(metal).("th"+string(ii)),25,0.00001);
end

figure(4);
hold on;
for ll=1:8
    plot(const.data.(metal).time,const.data.(metal).("th"+string(ll)));
end

yline(y);
xline(tss);

title("Temperature Progression of Metal");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
hold off;
legend (const.legend,'location','eastoutside');