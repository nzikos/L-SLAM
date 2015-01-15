function p=updateSLAM(p,d,th,u,t,phi)

for i=1:length(p)
    p(i).phi=phi+randn*0.03;
    [p(i).fs,p(i).pos,w]=updateParticle(p(i).fs,d,th,u,t,p(i).phi,p(i).pos);
    p(i).w=w;
end
sumw=sum([p.w]);
for i=1:length(p)
    p(i).w=p(i).w/sumw;
end
p=resample(p);



function [fs,pos,mw]=updateParticle(fs,d,th,u,t,phi,pos)
Q=getTransitionCov();
p=updatePos(pos(t-1).p,u,phi);
P=pos(t-1).P+Q;
Pfp=pos(t-1).Pfp;
pos(t-1).pii=p;
[idd,mw]=getAssociation(fs,d,th,phi,p);

if isempty(idd)  % New feature
    fs=[fs addFeature(d,th,t,p,P,phi)];
    pos(t).p=p; pos(t).pi=p; pos(t).P=P; pos(t).Pi=P; pos(t).Q=Q; pos(t).Pfp=Pfp;
else             % Update feature
    [fs(idd),K,y]=updateFeature(fs(idd),d,th,t,p,P,phi);
    pos(t).p=fs(idd).p; pos(t).pi=fs(idd).p;
    pos(t).P=fs(idd).Pp; pos(t).Pi=fs(idd).Pp; pos(t).Q=Q; pos(t).Pfp=fs(idd).Pfp;
    [fs,pos]=smooth(fs,pos,t);
end

function [fs,pos]=smooth(fs,pos,st)
tsir=[fs.t];
H=-[-1 0 0 0;0 -1 0 0 ];
Q=getTransitionCov();
for t=st-1:-1:2
    
    K=pos(t).Pi/(pos(t).Pi+Q)*(pos(t+1).P-pos(t).Pi-Q)/(pos(t).Pi+Q);
    y=(pos(t).Pi+Q)/(pos(t+1).P-pos(t).Pi-Q)*(pos(t+1).p-pos(t).pii);
    
    if sqrt(abs(det(K)))<0.015 break; end
    
    pos(t).p=pos(t).pi+K*y;
    temp=(eye(2)+K)*pos(t).Pi;
    
    pos(t).P=temp;
    idx=find(tsir==t);
    if isempty(idx) continue; end
    f=fs(idx(1));
    
    Ktemp=[K; f.Pfp'*K'/f.Pp];
    x=[f.p;f.m]+Ktemp*y;
    Ptemp=(eye(4)+Ktemp*H)*[f.Pp f.Pfp;f.Pfp' f.Pf];
    
    f.Pp_s  = Ptemp(1:2,1:2);
    f.Pf_s  = Ptemp(3:4,3:4);
    f.Pfp_s = Ptemp(1:2,3:4);
    f.m_s   = x(3:4);
    f.p_s   = x(1:2);
    fs(idx(1))=f;
end

function f=addFeature(d,th,t,p,P,phi)
R=getr(d,th+phi);
f.t=t;
f.p=p;
f.Pp=P;
f.Pfp=P;
f.Pf=P+R;
f.m=p+d*[cos(th+phi);sin(th+phi)];

f.p_s=p;
f.Pp_s=P;
f.Pfp_s=P;
f.Pf_s=P+R;
f.m_s=p+d*[cos(th+phi);sin(th+phi)];

function [f,K,y]=updateFeature(f,d,th,t,p,P,phi)
Q=getTransitionCov();
R=getr(d,th+phi);

pfp=f.Pfp_s;
P=[f.Pp_s+Q pfp; pfp' f.Pf_s];

f.t=t;

z=d*[cos(th+phi);sin(th+phi)];
x=[p;f.m_s];
H=[-1 0 1 0;0 -1 0 1];
y=z-H*x;

S=H*P*H'+R;
K=P*H'/S;
x=x+K*y;
P=(eye(4)-K*H)*P;
f.Pp=P(1:2,1:2);
f.Pf=P(3:4,3:4);
f.Pfp=P(1:2,3:4);
f.m=x(3:4);
f.p=x(1:2);
K=K(1:2,1:2);

function a=checkCov(a)
a=(a+a')/2;
b=sqrt(diag(a)*diag(a)');
mask=abs(a)>b;
a(mask)=a(mask)./abs(a(mask)).*b(mask)*0.9;

function [idx,mw]=getAssociation(fs,d,th,phi,p)
idx=[];
p0=10^-5;
mw=p0;
if isempty(fs) return; end
R=getr(d,th+phi);
H=[-1 0 1 0;0 -1 0 1];
z=d*[cos(th+phi);sin(th+phi)];
H2=-[-1 0 0 0;0 -1 0 0 ];
Q=getTransitionCov();
for i=1:length(fs)
    f=fs(i);

    Ptemp=[f.Pp_s+Q f.Pfp_s;f.Pfp_s' f.Pf_s];
    x=[p;f.m_s];
    y=z-H*x;
    S=H*Ptemp*H'+R;
    w(i)=exp(-0.5*y'/S*y)*(abs(det(S))^-0.5) /(2*pi);
end

[mw,idx]=max(w);
if mw<p0
    idx=[];
    mw=p0;
end




