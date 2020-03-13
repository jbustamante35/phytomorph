term.closeDB
term.openDB
a = term(1,1);
b = term(1,1);
z = term(1,1);
c = a*b*z;
term.query
e = c + a;
d = c - a;
term.query
tmp = e.getTerm(1);
%%
as = simplex(rand(1,2),rand(1,3,4),1,1);
bs = simplex(rand(1,2),rand(1,3,4),1,1);
as2 = simplex(rand(1,2),rand(1,3,4),1,1);
bs2 = simplex(rand(1,2),rand(1,3,4),1,1);
cs = as + bs;
ds = as*bs;
ds2 = as2*bs2;
delta = ds - ds2;

