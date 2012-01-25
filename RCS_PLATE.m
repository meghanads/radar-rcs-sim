%
%
%   Function: RCS_PLATE
%	Description: calculates RCS of Scatterers.
%
%

function Y = RCS_PLATE(len, angle, width, rot_l, rot_h, rot_step, len_dipole, k, p_dist)

    n=1;

    for a= rot_l:rot_step:rot_h
        
        phi_angle=angle-a;    % angle made by plate with  y-axis
        
        proj=len*cos(pi*(phi_angle)/180);
        
        num_dipole=ceil(proj/len_dipole);
        
        if(phi_angle ==0)
            
            phi_dipole=0;
            
        elseif(mod(phi_angle, 90)==0)
            
            phi_dipole=0;
            
        else
            
            phi_dipole=len_dipole*tan(pi*phi_angle/180);
            
        end
        
        e_dipole(n)=0;
        
        for m=1:num_dipole
            
            e_dipole(n)=e_dipole(n) + exp(-2*j*(phi_dipole*(m-1)+p_dist)*k);
            
        end
        
        e_dipole(n)=e_dipole(n)*width*len_dipole;
        
        n=n+1;
        
    end
    
    Y=e_dipole;
    
        
        
    
