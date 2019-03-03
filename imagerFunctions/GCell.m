%==============================================================================%
% GCell()                                                                      %
%==============================================================================%
% [Tx][Rx][f],[Sx] => {Tx,Rx,Sx}. Note, this is NOT the order of radio [Rx][Tx]

function c=GCell(g,fCount,txCount,rxCount,sxCount)
    c=cell(txCount,rxCount,sxCount);
    for i=1:size(c,1) % Tx
        for j=1:size(c,2) % Rx
            for k=1:size(c,3) % Sx
                c{i,j,k}=g(fCount*(rxCount*(i-1)+(j-1))+(1:fCount),k);
            end
        end
    end
end