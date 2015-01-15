function p=resample(p)

N= length(p);
w= [p.w];

[keep, Neff] = stratified_resample(w);

if Neff < N*0.75
    p= p(keep);
    for i=1:N, p(i).ws= 1/N; end
end


function [keep, Neff] = stratified_resample(w)

w= w / sum(w); % normalise
Neff= 1 / sum(w .^ 2); 

len= length(w);
keep= zeros(1,len);
select = stratified_random(len); 
w= cumsum(w); 

ctr=1; 
for i=1:len
   while ctr<=len & select(ctr)<w(i)
       keep(ctr)= i;
       ctr=ctr+1; 
   end
end



function s = stratified_random(N)

k = 1/N;
di = (k/2):k:(1-k/2); % deterministic intervals
s = di + rand(1,N) * k - (k/2); % dither within interval