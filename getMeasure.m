function [d,th,idx]=getMeasure(map,p,phi)

T=[cos(phi) -sin(phi) p(1);sin(phi) cos(phi) p(2);0 0 1];

map2=T\[map; ones(1,size(map,2))];

map2=map2(1:2,:);

d=sqrt(sum(map2.^2));
th=atan2(map2(2,:),map2(1,:));

mask=d>5 & d<50 & ((th<pi/2 & th>pi/6) | (th>-pi/2 & th<-pi/6));

d(~mask)=1000;

[~,idx]=sort(d,'ascend');
d=d(idx(1));
th=th(idx(1));
idx=idx(1);
