function v = deflection(w0,w1,L,E,I,x)
C1 = (w1-w0)/(120*L);
C2 = w0/24;
C3 = ((w0+w1)*L)/12;
C4 = (L^2)*(w0+(2*w1))/12;
C5 = 1/(E*I);
v = C5*(-(C1*x^5)-(C2*x^4)+(C3*x^3)-(C5*x^2));
end