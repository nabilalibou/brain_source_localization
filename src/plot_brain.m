function plot_brain(S, path_png, mesh, t)

         figure();

         set(gcf, 'PaperUnits', 'inches');
         x_width=4 ;y_width=4;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]);trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3), S);
        
        
                        
        title(t,'FontSize',16); 
        axis off;
        az=-158.7;el= 40.4;
        view (az, el)
        saveas(gcf,path_png) 