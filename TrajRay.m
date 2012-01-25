%
%
%   Function:TragRay
%	Description: Returns Rayleigh Tragectory 
%
%

function Y = TrajRay(Tsim, Trep, P_velo, P_dist)

    x = 0:Trep:(Tsim-Trep);
    num=Tsim/Trep;
    
    p = raylpdf(x,1);
    p=p*P_velo;
    p=p+P_velo;
    
    
    subplot(1,2,1);
    plot(x,p);
    Title('Input Velocity Trajectory');
    xlabel('Time(Sec.)')
    ylabel('velocity(m/s)');
    
     Dist(1)=P_dist;
    
    for i=2:num
        Dist(i)=Dist(i-1)-p(i-1)*Trep;
    end
    
    subplot(1,2,2);
    plot(x,Dist);
    Title('Input Distance Trajectory');
    xlabel('Time(Sec.)')
    ylabel('distance(m)');
    
    Y=p;
    
end

