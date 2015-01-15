clear all;clc;
% stream = RandStream.getGlobalStream;
% reset(stream);
% rand(1,19);
for j=1:1
    map=rands(2,30)*30;
    rp=[0;-15]; p=rp; phi=0; fs=[]; P=zeros(2,2); Q=eye(2)*0.2^2;
    pos=[]; pos.p=p; pos.pi=p; pos.P=P; pos.Pfp=P; pos.Pi=P; pos.rp=rp; pos.Q=Q; mov=[]; ismov=1;
    for i=1:5
        prt(i).pos=pos; prt(i).fs=fs;
        prt(i).w=1; prt(i).phi=phi;
    end
    
    for i=2:70        
        u=2; phi=phi+2*pi/50;
        rp=updatePos(rp,u,phi); pos(i).rp=rp;
        [d,th,idx(i-1)]=getMeasure(map,rp,phi);
        
        prt=updateSLAM(prt,d+randn*0.1,th+randn*0.03,u+randn*0.12,i,phi);
        [~,midx]=max([prt.w]);
        mov=visualize(pos,prt,midx,map,mov,ismov,1);
    end
    err=[pos.rp]-[prt(midx).pos.p];
    res(j)=mean(sqrt(sum(err.^2)));
end
mean(res)
% mean(res2)
if ismov movie2avi(mov,'lslam.avi');end