x = 0:1/360:1;

figure(1)

filename = 'avaModClusters.gif';

for n = 1:1000
    n
      
      scatter3(maptsne(:,1), maptsne(:,2), maptsne(:,3), 30, rAvaSizes);
      view(x(n)*360, 5);
      drawnow

      frame = getframe(1);

      im = frame2im(frame);

      [imind,cm] = rgb2ind(im,256);

      if n == 1;

          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);

      else

          imwrite(imind,cm,filename,'gif','WriteMode','append');

      end

end