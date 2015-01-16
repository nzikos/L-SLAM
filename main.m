clear all;clc;
addpath('tools')

noOfParticles=5;
map=rands(2,30)*30;

realPosition=[0;-15]; p=realPosition; phi=0; fs=[]; P=zeros(2,2); Q=eye(2)*0.2^2;
pos=[]; pos.p=p; pos.pi=p; pos.P=P; pos.Pfp=P; pos.Pi=P; pos.rp=realPosition; pos.Q=Q; 
mov=[]; saveAvi=1;


for i=1:noOfParticles
    prt(i).pos=pos; prt(i).fs=fs;
    prt(i).w=1; prt(i).phi=phi;
end

for i=2:70
    u=2; phi=phi+2*pi/50;
    realPosition=updatePos(realPosition,u,phi); pos(i).rp=realPosition;
    [d,th,~]=getMeasure(map,realPosition,phi);
    
    prt=updateSLAM(prt,d+randn*0.1,th+randn*0.03,u+randn*0.12,i,phi);
    [~,midx]=max([prt.w]);
    mov=visualize(pos,prt,midx,map,mov,saveAvi,1);
end
err=[pos.rp]-[prt(midx).pos.p];
res=mean(sqrt(sum(err.^2)));
mean(res)

if saveAvi movie2avi(mov,'lslam.avi');end