function [ g ] = aboveMeanG( x, meandB)
if iscell(x)
   g= GCell(db(GVector(x))-meandB,101, 10,30,1);
else
    g=db(x)-meandB;
end
end

