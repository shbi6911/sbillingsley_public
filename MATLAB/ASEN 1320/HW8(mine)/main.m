%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

xData = -5:0.175:30;
yData = GenerateData (xData);
M = zeros(1,4);
B = zeros(1,4);

templog = xData < 10;
tempindex = find(templog);
xDatatemp = xData(tempindex);
yDatatemp = yData(tempindex);
[M(1),B(1)] = LeastSquares(xDatatemp,yDatatemp);

templog = (xData >= 10) & (xData < 15);
tempindex = find(templog);
xDatatemp = xData(tempindex);
yDatatemp = yData(tempindex);
[M(2),B(2)] = LeastSquares(xDatatemp,yDatatemp);

templog = (xData >= 15) & (xData < 20);
tempindex = find(templog);
xDatatemp = xData(tempindex);
yDatatemp = yData(tempindex);
[M(3),B(3)] = LeastSquares(xDatatemp,yDatatemp);

templog = xData >= 20;
tempindex = find(templog);
xDatatemp = xData(tempindex);
yDatatemp = yData(tempindex);
[M(4),B(4)] = LeastSquares(xDatatemp,yDatatemp);

xFit = -5.0:0.007:30;
yFit = PiecewiseLeastFit (M, B, xFit);

makePlot (xData, yData, xFit, yFit);
