function f = loading(w0,w1,L)
f = @(x) ((w1-w0)/L)*x + w0;
end