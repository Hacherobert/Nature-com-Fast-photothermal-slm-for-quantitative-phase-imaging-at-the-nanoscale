%Generation d'un profile HSD en forme de disque de rayon "rayon", Nx,Ny le
%nombre de pixel selon x,y de l'image, pix_size taille du pixel en mètre,
%Q_Ir chaleur emise par la source (on suppose distribution de source
%homogène) en W/m^2
%HSD= profil HSD
% Q_tot= chaleur totale emis par la source en W

function [HSD,Q_tot]=HSD_generation(Nx,Ny,radius,z_source,pix_size,Q_Ir)

    power_per_pix=Q_Ir*pix_size^2;              %Puissance émise par pixel
    z_s=(round(z_source/pix_size))*pix_size;

    Q_tot=0;

    HSD=zeros(Nx,Ny);
        
    for i=1:Nx

        for j=1:Ny

                if (sqrt((i-Nx/2)^2+(j-Ny/2)^2)<radius/pix_size) 

                    HSD(i,j)=power_per_pix;
                    Q_tot=Q_tot+power_per_pix;

                end

            

        end

    end
          
    
end