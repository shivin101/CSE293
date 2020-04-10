function [matchNN,matchDist,vis]=matchNeighbor(B,I,matchNN,matchDist,vis,idx)
    matchId = I(1);
    if vis(I(1))==1
        for i=1:size(matchNN,1)
            if(matchNN(i))==I(1)
                if(matchDist(i)>B(1))
                    matchDist(i)=9999;
                    matchNN(i)=0;
                    matchDist(idx)=B(1);
                    matchNN(idx)=I(1);
%                 else
                    
                end
            end
        end
    else
        vis(I(1))=1;
        matchNN(idx)=I(1);
        matchDist(idx)=B(1);
    end
end