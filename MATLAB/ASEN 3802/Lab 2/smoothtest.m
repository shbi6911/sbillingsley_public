%x=linspace(0,const.data.steel.time(end),((length(const.data.steel.time)-1)*10)+1);
%interpolated=interp1(const.data.steel.time,const.data.steel.th1,x);

smoothed=smoothdata(const.data.steel.th1,'gaussian',100);

figure(1);
hold on;
%plot(const.data.steel.time,const.data.steel.th1);
plot(const.data.steel.time,smoothed);



