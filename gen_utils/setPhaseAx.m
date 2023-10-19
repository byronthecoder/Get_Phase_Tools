function setPhaseAx(axIn)
    set(axIn,'ytick',[0,pi,2*pi],'yticklabel',{'0','\pi','2\pi'})
    set(axIn,'ylim', [-0.01,2*pi+0.01])
