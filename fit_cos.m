%Fit a curve (y) with cosine function
%s0 is the vector containing the initial values of the cosine function
%(amplitude, phase, period and offset)
% s is the vector containing the parameter of the fitted cosine function
% rsq is the correclation between the curve y and the fitted cosine function

function [s,s0,rsq] = fit_cos(y,amp,fixamp,per,fixper,offset,phase,fact_iter)
    N_sampling=length(y);
    x      = 1:N_sampling;                                               
    maxy   = max(y);
    miny   = min(y);

    %% estimation of the initial parameters
    % estimation of the amplitude
    if amp==0
        
        amp    = maxy-miny;                                             
        
    end
    
    % Estimation of the offset
    if offset==0
        
        offset = amp/2+miny;                                                 
        
    end
    
    yoff  = y-offset; 
    
    % Find zero-crossings
    zero_x = x(xor(yoff(:) <= 0,circshift(yoff(:),[1 0]) <= 0));       
    
    
   
    
   % Period estimation
    if per==0
        
        if length(zero_x)==1                                                

            ydiff=diff(yoff);
            ydiff=[ydiff(1);ydiff];

            zero_externa_length = max(abs(zero_x-x(logical((yoff == max(yoff)) + (yoff == min(yoff))))));        
            per = 4*zero_externa_length; 

        else

            per = 2*mean(diff(zero_x));

        end
    end
    
    % Phase estimation
    if phase==0

        nb_phase_ech=round(per/20);


        if isempty(zero_x)

            if mean(yoff)<0

                phase=3*pi/2;
            else

                phase=pi/2;
            end

        else
            if zero_x(1) == 1

                zero_x = zero_x(2:end);

                if zero_x(1) > per/2

                    zero_x(1)=per/2;
                end

            end

            if (zero_x(1)+nb_phase_ech)>N_sampling

               y(zero_x(1):N_sampling)=mean(y(1:zero_x(1))); 

               zero_x(1)= round(N_sampling/2);
               phase = pi/2-zero_x(1)*2*pi/per;

            else      




                if mean(diff(yoff(zero_x(1):zero_x(1)+nb_phase_ech)))>0

                    phase = pi*3/2-zero_x(1)*2*pi/per;

                else

                    phase = pi/2-zero_x(1)*2*pi/per;

                end
            end

        end
    end
%% Least square mean minimization

    if fixamp==1
        
        fit = @(b,x)  amp/2.*(cos(2*pi/per*x + b(1))) + b(2);                                           % Function to fit
        fcn = @(b) sum((fit(b,x) - y').^2);                                                             % Least-Squares cost function
        s0  = [ phase;  offset];                                                                        % Initial parameters
        options = optimset('MaxFunEvals',200*fact_iter*length(s0),'MaxIter',200*fact_iter*length(s0));
        s   = fminsearch(fcn, s0,options);                                                              % Minimise Least-Squares 
        rsq = corrcoef(y, fit(s,x)).^2;   
        s   = [amp per s(1) s(2)];
    else
    
        if fixper==0
         
            fit = @(b,x)  b(1).*(cos(2*pi/b(2)*x + b(3))) + b(4);                                               % Function to fit
            fcn = @(b) sum((fit(b,x) - y').^2);                                                                 % Least-Squares cost function
            s0  = [amp/2;  per;  phase;  offset];                                                               % Initial parameters
            options = optimset('MaxFunEvals',200*fact_iter*length(s0),'MaxIter',200*fact_iter*length(s0));
            s   = fminsearch(fcn, s0,options);                                                                  % Minimise Least-Squares    
            rsq = corrcoef(y, fit(s,x)).^2;
   
        else

            fit = @(b,x)  b(1).*(cos(2*pi/per*x + b(2))) + b(3);                                               % Function to fit
            fcn = @(b) sum((fit(b,x) - y').^2);                                                                 % Least-Squares cost function
            s0  = [amp/2; phase;  offset];                                                                      % Initial parameters
            options = optimset('MaxFunEvals',200*fact_iter*length(s0),'MaxIter',200*fact_iter*length(s0));
            s   = fminsearch(fcn, s0,options);                                                                 % Minimise Least-Squares 
            rsq = corrcoef(y, fit(s,x)).^2;   
            s   = [s(1) per s(2) s(3)];


        end
    end

end
