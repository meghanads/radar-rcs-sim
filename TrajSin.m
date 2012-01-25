%
%
%   Function:TragConst
%	Description: Returns constant Tragectory(i.e. 
%
%

function Y = TrajSin(Tsim, Trep, P_velo, P_dist)

    x=0:Trep:(Tsim-Trep);
    num=Tsim/Trep;
    
    step=(2*pi)/num;
    m=-pi:step:pi-step;
       
    Traj=P_velo*cos(2*m)*3/2;

    Traj=Traj + 2*P_velo;
    size(x)
    size(Traj)
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

