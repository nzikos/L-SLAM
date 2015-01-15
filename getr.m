function r=getr(d,f)
s1=0.2;
s2=0.05;

s11 = ((cos(f)*s1)^2 + (d*sin(f)*s2)^2);
s22 = ((sin(f)*s1)^2 + (d*cos(f)*s2)^2);
s12 = cos(f)*sin(f) * (s1^2 - d^2*s2^2);

r=[s11 s12; s12 s22];

