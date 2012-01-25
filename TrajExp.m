%
%
%   Function:TragExp
%	Description: Returns Exponential Tragectory
%
%

function Y = TrajExp(Tsim, Trep, P_velo, P_dist)

    x=0:Trep:(Tsim-Trep);
    num=Tsim/Trep;
       
    Traj=P_velo*exp(-0.4*x);
    Traj=Traj+P_velo;
    subplot(1,2,1);
    plot(x,Traj);
    Title('Input Velocity Trajectory');
    xlabel('Time(Sec.)')
    ylabel('velocity(m/s)');
    
    Dist(1)=P_dist;
    
    for i=2:num
        Dist(i)=Dist(i-1)-Traj(i-1)*Trep;
    end
    
    subplot(1,2,2);
    plot(x,Dist);
    Title('Input Distance Trajectory');
    xlabel('Time(Sec.)')
    ylabel('distance(m)');
    
    Y=Traj;
    

end

