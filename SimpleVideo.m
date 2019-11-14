%aviobj = avifile(['comparison'],'fps',15);

wr = VideoWriter('comparison_Re.avi');
wr.FrameRate = 15;
open(wr);

fig = figure('visible','off','color','white','Units','Normalized','Position',[0,0,0.5,0.5],'renderer','zbuffer');
s1 = axes('NextPlot','replacechildren','FontSize',12,'FontName','sansserif','DataAspectRatio',[1 1 1],'XGrid','on','box','on');

colorbar
lx100 = out.velocity.ux(2:2:end,2:2:end,2:2:end);
lx50 = out2.velocity.ux(1:end,1:end,2:end);

blub = permute(lx100-lx50,[2 1 3]);



for i=2:600,
    
    



    
    
    h1 = imagesc(blub(:,:,i),[-0.05 0.05]);
    
    %wr = addframe(wr,fig);
    writeVideo(wr,getframe(fig))
    %print(fig,'-dpng',['plots','/shore-',filename,'T', num2str(round(out.time_y(i))),'.png'],'-r80')
    delete(h1)
    
    i
   
end

%saveas(fig,['plots','/Depo-',filename,'T' num2str(round(out.time_y(t(end)))),'.fig'])
close(wr)
%aviobj = close(aviobj);
close(fig)