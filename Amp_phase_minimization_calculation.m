
%%% Estimate the amplitude (amp), phase(phase), offset (offsett) of a curve (im)from a fitting
%%% algorithment based and least square minimization. error is the
%%% correlation between a cos function and the curve. im can be 3D (same image with different phase-shift)or
%%% 4D (stack of image with different phase-shift), amp,phase,offset,error will be respectively 1D,2D and 3D.
%%% To avoid local minimum, the amplitude (amp), period (per), offset
%%% (offset) and phase (phase) can be set. If not, they have to be set to
%%% 0. 
%%% fixamp=1 does the minimization with a fix value of amp otherwise it has to be set to 0,
%%% fixper=1 does the minimization with a fix value of per otherwise it has to be set to 0,

function [amp,phase,offset,error]=Amp_phase_minimization_calculation(im,amp,fixamp,per,fixper,offset,phase,fact_iter)

dim_im = length(size(im));

    if dim_im == 4
    
        length_im_1  = size(im,1);
        length_im_2  = size(im,2);
        Nb_zscan     = size(im,3);
        Nb_sampling  = size(im,4);
        amp          = zeros(length_im_1,length_im_1,Nb_zscan);
        phase        = zeros(length_im_1,length_im_1,Nb_zscan);
        offset       = zeros(length_im_1,length_im_1,Nb_zscan);
        error        = zeros(length_im_1,length_im_1,Nb_zscan);
        
        for i=1:Nb_zscan 
            i
            im = reshape(im(:,:,i,:),[length_im_1,length_im_2,Nb_sampling]);           
            
            for k=1:length_im_1  

                for l=1:length_im_2
                    
                    y=reshape(im(k,l,:),[Nb_sampling,1]);
                    [s,~,rsq]     = fit_cos(y,amp,fixamp,per,fixper,offset,phase,fact_iter); 
                 
                    if (rsq(2,1)<0.2)&&((k>1)&&(l>1))
                        
                        s(1)=(amp(k-1,l,i)+amp(k,l-1,i)+amp(k-1,l-1,i))/3;
                        s(2)=(phase(k-1,l,i)+phase(k,l-1,i)+phase(k-1,l-1,i))/3;
                        
                    end
                        
                    
                     amp(k,l,i)   = s(1);  
                     phase(k,l,i) = s(2);  
                     error(k,l,i) = rsq(2,1);
                     offset(k,l,i)= s(3);

                end

            end

        end
        
    elseif dim_im==3
        
        length_im_1=size(im,1);    
        length_im_2=size(im,2);        
        Nb_sampling=size(im,3); 
        
        
        amp      = zeros(length_im_1,length_im_1);
        phase    = zeros(length_im_1,length_im_1);
        error    = zeros(length_im_1,length_im_1);

        
        for k=1:length_im_1  

            for l=1:length_im_2
                
                y=reshape(im(k,l,:),[Nb_sampling,1]);
                [s,~,rsq]  = fit_cos(y,amp,fixamp,per,fixper,offset,phase,fact_iter);                  
                amp(k,l)   = s(1);  
                phase(k,l) = s(2);  
                error(k,l) = rsq(2,1);         

            end

        end
      
    

       
    end



end

