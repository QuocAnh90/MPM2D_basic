%% Function test_new_bspline
% Input: Knot vector, degree and global positions
% Output: Value of material points in basis functions

function [N_vec , B_vec] = test_new_bspline(Xi,deg,pos_p_glob)

Xp_repmat = repmat(pos_p_glob, [1, length(Xi)-1]);
Xi_repmat = repmat(Xi', [length(pos_p_glob), 1]);

% Create N_i,0 for all material points
N_vec = Xi_repmat(:,1:end-1) <= Xp_repmat &...
        Xp_repmat < Xi_repmat(:,2:end);
B_vec = zeros(length(pos_p_glob),length(Xi)-1);

    for p=1:deg,
        N_L = (   Xp_repmat(:,1:end-p)-Xi_repmat(:,1:end-p-1) )./...
            ( Xi_repmat(:,1+p:end-1)-Xi_repmat(:,1:end-p-1) ).*N_vec(:,1:end-1);
        N_L(isnan(N_L)) = 0;
        N_R = ( Xi_repmat(:,1+p+1:end)-Xp_repmat(:,1:end-p) )./...
            ( Xi_repmat(:,1+p+1:end)-Xi_repmat(:,1+1:end-p) ).*N_vec(:,2:end);
        N_R(isnan(N_R))=0;
        
        
        B_L = ( p*ones(length(pos_p_glob),length(Xi)-p-1))./...
            ( Xi_repmat(:,1+p:end-1)-Xi_repmat(:,1:end-p-1) ).*N_vec(:,1:end-1);
        B_L(isnan(B_L)) = 0;
        B_R = ( p*ones(length(pos_p_glob),length(Xi)-p-1) )./...
              ( Xi_repmat(:,1+p+1:end)-Xi_repmat(:,1+1:end-p) ).*N_vec(:,2:end);
        B_R(isnan(B_R))=0;
        
        N_vec = N_L + N_R;
        B_vec = B_L - B_R;
    end
%    N_vec_new = sparse(N_vec);
%    B_vec_new = sparse(B_vec);
%    clear N_L N_R B_L B_R N_vec B_vec
end
