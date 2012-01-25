% randomize
% FUNCTION: randomize
% generate number 1 or 0 randomely


function Y=RANDOM()
if((rand(1)-0.5)*2>0)
    Y=1;
else
    Y=-1;
    
end