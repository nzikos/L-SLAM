function mov=visualize(rp,prt,midx,map,mov,ismov,isplot)

pos=prt(midx).pos;
fs=prt(midx).fs;

if ~isplot
    return
end
temp=[rp.rp];
plot(temp(1,:),temp(2,:),'--rs')
axis([-30 30 -30 30]);axis square;hold on;grid on
for j=1:length(pos)
    plot_gaussian_ellipsoid(pos(j).p, pos(j).P,2);
    plot(pos(j).p(1),pos(j).p(2),'r*')
end

for j=1:length(fs)
    plot_gaussian_ellipsoid(fs(j).m_s, fs(j).Pf_s,2);
    plot(map(1,:),map(2,:),'g*')
end

for i=1:length(prt)
    plot(prt(i).pos(end).p(1),prt(i).pos(end).p(2),'b*')
end

hold off

if ismov
    mov=[mov getframe];
end
pause(0.01)