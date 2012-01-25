%
%
%   Function: RCS_SCATT
%
%       Description: calculates RCS for Scatterers
%                       
% Should:
%           abs(ang_xaxis)<=180
%
%

function Y = RCS_SCATT(rcs, ref_dist, dist, ang_xaxis, rot_l, rot_h, rot_step, k, len_dipole)

        n=1;
        for a=rot_l:rot_step:rot_h
            
            phi_angle=(ang_xaxis-a);
            
            if(phi_angle==0)
                
                e_scat(n)=exp(-2*j*(ref_dist - dist)*k);
                
                            
            elseif(mod(phi_angle,180)==0)
                
                e_scat(n)=exp(-2*j*(ref_dist+dist)*k);
                
                
            elseif(mod(phi_angle,90)==0)
                
                e_scat(n)=exp(-2*j*ref_dist*k);
                
                
            elseif(abs(phi_angle) < 90)
                
                e_scat(n)=exp(-2*j*(ref_dist-dist*cos(pi*phi_angle/180))*k);
                
            else
                
                e_scat(n)=exp(-2*j*(ref_dist+dist*cos(pi*phi_angle/180))*k);
            end
            
            
            e_scat(n)=e_scat(n)*(rcs);
            %_scat(n)=e_scat(n)*len_dipole;
            
            n=n+1;
            
        end
        
        Y=e_scat;
end

                
            
            