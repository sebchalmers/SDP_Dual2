function [ out ] = HP( P, IP, M )
%HPM transformation for the NT method
    out = 0.5*(P*M*IP + IP.'*M.'*P.');


end

