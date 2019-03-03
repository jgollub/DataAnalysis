%==============================================================================%
% GVector()                                                                    %
%==============================================================================%
% {Tx,Rx} => [Tx][Rx][f].  Note, this is NOT the order of radio [Rx][Tx]

function v=GVector(g)
    v=zeros(numel(g{1,1})*size(g,1)*size(g,2),size(g,3));
    for i=1:size(g,1) % Tx
        for j=1:size(g,2) % Rx
            for k=1:size(g,3) % Sx
                v(numel(g{1,1})*(size(g,2)*(i-1)+(j-1))+(1:numel(g{1,1})),k)=g{i,j,k};
            end
        end
    end
end